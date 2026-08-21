/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.platform

import java.util.regex.Pattern

import groovy.transform.CompileStatic

/**
 * Describe how resource label keys and values must be normalised to be accepted
 * by the API of a given executor backend.
 *
 * One policy is provided for each backend consuming the {@code resourceLabels}
 * directive, and it is meant to be applied to the labels derived from the workflow
 * metadata only — see {@link nextflow.processor.TaskConfig#getResourceLabels(ResourceLabelPolicy)}.
 * User-declared labels are never rewritten: an invalid value there is the user's own
 * responsibility, and silently changing it would be worse than the API error.
 *
 * Two rules are common to all policies:
 * <li>An entry whose key sanitises to an empty string is dropped, since there is no
 *   key left to apply it under.
 * <li>An entry whose *value* sanitises to an empty string is dropped only by the
 *   policies rejecting empty values in the position where the labels are applied,
 *   namely {@link #GOOGLE} (Batch resource labels) and {@link #K8S} (pod labels).
 *   AWS and Azure accept an empty value, therefore the entry is retained.
 *
 * Truncation is applied last and never leaves an invalid trailing character.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
abstract class ResourceLabelPolicy {

    /** No-op policy, for a backend accepting the labels verbatim e.g. the Seqera scheduler */
    static final ResourceLabelPolicy IDENTITY = new IdentityPolicy()

    /** AWS Batch job tags: permissive charset, key up to 128 chars, value up to 256 */
    static final ResourceLabelPolicy AWS = new AwsPolicy()

    /** Google Batch labels: lowercase {@code [a-z0-9_-]}, key starting with a letter, up to 63 chars */
    static final ResourceLabelPolicy GOOGLE = new GooglePolicy()

    /** Kubernetes pod labels: optional DNS subdomain prefix, alphanumeric bounded name and value */
    static final ResourceLabelPolicy K8S = new K8sPolicy()

    /** Azure Batch pool metadata: near identity, minus the reserved {@code microsoft} name prefix */
    static final ResourceLabelPolicy AZURE = new AzurePolicy()

    private final String name

    protected ResourceLabelPolicy(String name) {
        this.name = name
    }

    String getName() { name }

    /**
     * @param key The label key to be normalised
     * @return The key accepted by the target API, or an empty string when nothing is left of it
     */
    abstract String sanitizeKey(String key)

    /**
     * @param value The label value to be normalised
     * @return The value accepted by the target API, possibly empty
     */
    abstract String sanitizeValue(String value)

    /**
     * @return {@code true} when the target API rejects an empty label value
     */
    protected boolean dropEmptyValues() { false }

    /**
     * Normalise the given labels, dropping the entries that cannot be represented.
     *
     * @param labels The labels to be normalised
     * @return A new map holding the normalised entries, in the same order as the source map
     */
    Map<String,String> sanitize(Map<String,String> labels) {
        if( labels == null )
            return Collections.<String,String>emptyMap()
        if( labels.isEmpty() )
            return labels
        final result = new LinkedHashMap<String,String>(labels.size())
        for( Map.Entry<String,String> entry : labels.entrySet() ) {
            final key = sanitizeKey(entry.key)
            // no key, no label
            if( !key )
                continue
            final value = sanitizeValue(entry.value)
            if( !value && dropEmptyValues() )
                continue
            result.put(key, value)
        }
        return result
    }

    @Override
    String toString() {
        return "ResourceLabelPolicy[$name]"
    }

    /**
     * Truncate the given string to the specified maximum length
     */
    protected static String truncate(String value, int max) {
        return value.length() > max ? value.substring(0, max) : value
    }

    /**
     * Strip the leading and trailing characters that are not alphanumeric, as required
     * by the Kubernetes label syntax
     */
    protected static String trimNonAlphanumeric(String value) {
        int start = 0
        int end = value.length()
        while( start < end && !Character.isLetterOrDigit(value.charAt(start)) )
            start++
        while( end > start && !Character.isLetterOrDigit(value.charAt(end-1)) )
            end--
        return value.substring(start, end)
    }

    /**
     * AWS Batch tag syntax: letters, digits, spaces and {@code + - = . _ : / @},
     * key up to 128 characters, value up to 256.
     */
    @CompileStatic
    private static class AwsPolicy extends ResourceLabelPolicy {

        private static final Pattern INVALID = ~/[^A-Za-z0-9 +\-=._:\/@]/

        AwsPolicy() { super('aws') }

        @Override
        String sanitizeKey(String key) {
            return scrub(key, 128)
        }

        @Override
        String sanitizeValue(String value) {
            return scrub(value, 256)
        }

        private static String scrub(String value, int max) {
            if( !value )
                return ''
            final result = INVALID.matcher(value).replaceAll('_')
            // a trailing blank left over by the truncation would be dropped by the API anyway
            return truncate(result, max).stripTrailing()
        }
    }

    /**
     * Google Cloud label syntax: lowercase letters, digits, underscores and dashes;
     * the key must start with a lowercase letter; key and value up to 63 characters.
     * Therefore {@code nextflow.io/runName} becomes {@code nextflow_io_runname}.
     */
    @CompileStatic
    private static class GooglePolicy extends ResourceLabelPolicy {

        private static final Pattern INVALID = ~/[^a-z0-9_-]/

        GooglePolicy() { super('google') }

        @Override
        String sanitizeKey(String key) {
            if( !key )
                return ''
            final scrubbed = INVALID.matcher(key.toLowerCase()).replaceAll('_')
            // the key must begin with a letter, hence drop whatever comes before the first one
            int start = 0
            while( start < scrubbed.length() && !Character.isLetter(scrubbed.charAt(start)) )
                start++
            return truncate(scrubbed.substring(start), 63)
        }

        @Override
        String sanitizeValue(String value) {
            if( !value )
                return ''
            return truncate(INVALID.matcher(value.toLowerCase()).replaceAll('_'), 63)
        }

        @Override
        protected boolean dropEmptyValues() { true }
    }

    /**
     * Kubernetes label syntax: the key is an optional DNS subdomain prefix followed by
     * {@code /} and a name matching {@code [A-Za-z0-9]([-_.A-Za-z0-9]*)} of at most 63
     * characters — so {@code nextflow.io/runName} is applied verbatim. Values follow the
     * name syntax, therefore a URL value such as {@code https://github.com/foo/bar} is
     * stripped of its scheme and of its slashes.
     */
    @CompileStatic
    private static class K8sPolicy extends ResourceLabelPolicy {

        private static final Pattern SCHEME = ~/^[A-Za-z][A-Za-z0-9+.\-]*:\/\//
        private static final Pattern NAME_INVALID = ~/[^A-Za-z0-9_.\-]/
        private static final Pattern PREFIX_INVALID = ~/[^a-z0-9.\-]/

        K8sPolicy() { super('k8s') }

        @Override
        String sanitizeKey(String key) {
            if( !key )
                return ''
            final p = key.indexOf('/')
            if( p == -1 )
                return name0(key)
            final prefix = prefix0(key.substring(0, p))
            final name = name0(key.substring(p+1))
            if( !name )
                return ''
            return prefix ? prefix + '/' + name : name
        }

        @Override
        String sanitizeValue(String value) {
            if( !value )
                return ''
            return name0(SCHEME.matcher(value).replaceFirst(''))
        }

        @Override
        protected boolean dropEmptyValues() { true }

        private static String name0(String value) {
            final scrubbed = trimNonAlphanumeric(NAME_INVALID.matcher(value).replaceAll('_'))
            // trim again, since the truncation can uncover a non alphanumeric trailing char
            return trimNonAlphanumeric(truncate(scrubbed, 63))
        }

        private static String prefix0(String value) {
            final scrubbed = trimNonAlphanumeric(PREFIX_INVALID.matcher(value.toLowerCase()).replaceAll('-'))
            return trimNonAlphanumeric(truncate(scrubbed, 253))
        }
    }

    /**
     * Azure Batch pool metadata syntax: the labels are applied as-is, apart from the
     * {@code microsoft} name prefix which is reserved by the service.
     */
    @CompileStatic
    private static class AzurePolicy extends ResourceLabelPolicy {

        private static final Pattern RESERVED = ~/(?i)^microsoft[\s._\-\/]*/

        AzurePolicy() { super('azure') }

        @Override
        String sanitizeKey(String key) {
            if( !key )
                return ''
            return RESERVED.matcher(key).replaceFirst('')
        }

        @Override
        String sanitizeValue(String value) {
            return value ?: ''
        }
    }

    /**
     * Identity policy for a backend applying the labels verbatim
     */
    @CompileStatic
    private static class IdentityPolicy extends ResourceLabelPolicy {

        IdentityPolicy() { super('identity') }

        @Override
        String sanitizeKey(String key) {
            return key ?: ''
        }

        @Override
        String sanitizeValue(String value) {
            return value ?: ''
        }

        @Override
        Map<String,String> sanitize(Map<String,String> labels) {
            return labels != null ? labels : Collections.<String,String>emptyMap()
        }
    }

}
