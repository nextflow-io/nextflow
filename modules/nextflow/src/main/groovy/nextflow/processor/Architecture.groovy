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

package nextflow.processor

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.Memoized

/**
 * Implements the {@code arch} process directive, to hold information on the
 * CPU (micro)architecture required by the process.
 *
 * <p>Supports multiple comma-separated architectures (e.g. {@code arch 'linux/amd64,linux/arm64'}).
 * Multi-arch is fully supported by selected executors (e.g. Seqera) via {@link #platforms()} and
 * {@link #containerPlatform()}. Other executors use {@link #getDockerArch()} and {@link #getSpackArch()},
 * which return only the first (primary) architecture.
 *
 * @author Marco De La Pierre <marco.delapierre@gmail.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@EqualsAndHashCode
@CompileStatic
class Architecture {

    @EqualsAndHashCode
    @CompileStatic
    static class ArchEntry {
        final String value

        static protected String normalize( String name ) {
            def chunks = name.tokenize('/')
            if( chunks.size() == 3 )
                return chunks[1] + '/' + chunks[2]
            else if( chunks.size() == 2 )
                return chunks[1]
            else
                return chunks[0]
        }

        static private void validate( String value, String name ) {
            if( value == 'x86_64' || value == 'amd64' )
                return
            if( value == 'aarch64' || value == 'arm64' || value == 'arm64/v8' )
                return
            if( value == 'arm64/v7' )
                return
            throw new IllegalArgumentException("Not a valid `arch` value: ${name}")
        }

        static ArchEntry parse(String name) {
            final value = normalize(name)
            validate(value, name)
            return new ArchEntry(value)
        }

        private ArchEntry(String value ) {
            this.value = value
        }

        @Override
        String toString() {
            return value
        }
    }

    private final List<ArchEntry> entries

    private final String target

    /**
     * @return all architectures as Docker platform strings (e.g. {@code ['linux/amd64', 'linux/arm64']})
     */
    List<String> platforms() {
        entries.collect(it -> toDockerArch(it))
    }

    /**
     * @return all architectures as a comma-separated Docker platform string
     *         (e.g. {@code 'linux/amd64,linux/arm64'})
     */
    String containerPlatform() {
        platforms().join(',')
    }

    static private String toDockerArch(ArchEntry arch) {
        final value = arch.value
        if( value == 'x86_64' || value == 'amd64' )
            return 'linux/amd64'
        if( value == 'aarch64' || value == 'arm64' || value == 'arm64/v8' )
            return 'linux/arm64'
        if( value == 'arm64/v7' )
            return 'linux/arm64/v7'
        return null
    }

    /**
     * @return the Docker platform string for the first (primary) architecture
     */
    @Memoized
    String getDockerArch() {
        return toDockerArch(entries[0])
    }

    /**
     * @return the Spack-compatible architecture name for the first (primary) architecture,
     *         or the explicit {@code target} microarchitecture if specified
     */
    @Memoized
    String getSpackArch() {
        if( target != null )
            return target
        final value = entries[0].value
        if( value == 'x86_64' || value == 'amd64' )
            return 'x86_64'
        if( value == 'aarch64' || value == 'arm64' || value == 'arm64/v8' )
            return 'aarch64'
        return null
    }

    Architecture( String value ) {
        this(Map.of('name', value))
    }

    Architecture( Map res ) {
        if( !res.name )
            throw new IllegalArgumentException("Missing architecture `name` attribute")

        this.target = res.target as String
        this.entries = (res.name as String)
            .tokenize(',')
            .collect(it -> ArchEntry.parse(it.trim()) )
    }

    @Override
    String toString() {
        return platforms().join(',')
    }

}
