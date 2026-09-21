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

package nextflow.script.dsl

import java.util.regex.Pattern

import groovy.transform.CompileStatic

/**
 * Stateless selector matching shared by process directives and agent-only options.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
final class ConfigSelectorResolver {

    private static final String LABEL_PREFIX = 'withLabel:'
    private static final String NAME_PREFIX = 'withName:'

    static class SelectorMatch {
        final String rule
        final Object settings

        SelectorMatch(String rule, Object settings) {
            this.rule = rule
            this.settings = settings
        }
    }

    private ConfigSelectorResolver() {}

    /**
     * Return matching selector bodies in precedence order: labels first, followed by each
     * distinct process/agent name from least to most specific.
     */
    static List<Object> matchingSettings(Map<String,?> scope, List<String> labels,
            String baseName, String simpleName, String fullyQualifiedName) {
        final result = new ArrayList<Object>()
        result.addAll(matchingLabelSettings(scope, labels))
        for( final name : distinctNames(baseName, simpleName, fullyQualifiedName) )
            result.addAll(matchingNameSettings(scope, name))
        return result
    }

    static List<Object> matchingLabelSettings(Map<String,?> scope, List<String> labels) {
        return matchingLabelSelectors(scope, labels).collect { it.settings }
    }

    static List<SelectorMatch> matchingLabelSelectors(Map<String,?> scope, List<String> labels) {
        final result = new ArrayList<SelectorMatch>()
        for( final entry : scope.entrySet() ) {
            final rule = entry.key.toString()
            if( rule.startsWith(LABEL_PREFIX)
                    && matchesLabels(labels, rule.substring(LABEL_PREFIX.length()).trim()) )
                result.add(new SelectorMatch(rule, entry.value))
        }
        return result
    }

    static List<Object> matchingNameSettings(Map<String,?> scope, String name) {
        return matchingNameSelectors(scope, name).collect { it.settings }
    }

    static List<SelectorMatch> matchingNameSelectors(Map<String,?> scope, String name) {
        final result = new ArrayList<SelectorMatch>()
        for( final entry : scope.entrySet() ) {
            final rule = entry.key.toString()
            if( rule.startsWith(NAME_PREFIX)
                    && matchesName(name, rule.substring(NAME_PREFIX.length()).trim()) )
                result.add(new SelectorMatch(rule, entry.value))
        }
        return result
    }

    static List<String> distinctNames(String baseName, String simpleName, String fullyQualifiedName) {
        final result = new ArrayList<String>(3)
        for( final name : [baseName, simpleName, fullyQualifiedName] ) {
            if( name && !result.contains(name) )
                result.add(name)
        }
        return result
    }

    static boolean matchesLabels(List<String> labels, String pattern) {
        final isNegated = pattern.startsWith('!')
        if( isNegated )
            pattern = pattern.substring(1).trim()

        final regex = Pattern.compile(pattern)
        for( final label : labels ) {
            if( regex.matcher(label).matches() )
                return !isNegated
        }
        return isNegated
    }

    static boolean matchesName(String name, String pattern) {
        final isNegated = pattern.startsWith('!')
        if( isNegated )
            pattern = pattern.substring(1).trim()
        return Pattern.compile(pattern).matcher(name).matches() ^ isNegated
    }
}
