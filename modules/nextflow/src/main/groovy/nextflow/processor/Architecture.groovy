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
 * @author Marco De La Pierre <marco.delapierre@gmail.com>
 */
@EqualsAndHashCode
@CompileStatic
class Architecture {

    @EqualsAndHashCode
    @CompileStatic
    static class Arch {
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

        static Arch parse( String name ) {
            final value = normalize(name)
            validate(value, name)
            return new Arch(value)
        }

        private Arch( String value ) {
            this.value = value
        }

        @Override
        String toString() {
            return value
        }
    }

    private final List<Arch> archs
    private final String target

    List<String> platforms() {
        archs.collect(it -> toDockerArch(it))
    }

    String containerPlatform() {
        platforms().join(',')
    }

    static private String toDockerArch( Arch arch ) {
        final value = arch.value
        if( value == 'x86_64' || value == 'amd64' )
            return 'linux/amd64'
        if( value == 'aarch64' || value == 'arm64' || value == 'arm64/v8' )
            return 'linux/arm64'
        if( value == 'arm64/v7' )
            return 'linux/arm64/v7'
        return null
    }

    @Memoized
    String getDockerArch() {
        return toDockerArch(archs[0])
    }

    @Memoized
    String getSpackArch() {
        if( target != null )
            return target
        final value = archs[0].value
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
        this.archs = (res.name as String)
            .tokenize(',')
            .collect(it -> Arch.parse(it.trim()) )
    }

    @Override
    String toString() {
        return platforms().join(',')
    }

}
