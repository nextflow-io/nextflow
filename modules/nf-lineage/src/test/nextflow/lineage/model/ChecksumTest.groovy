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

package nextflow.lineage.model

import java.nio.file.Files
import java.nio.file.Path

import nextflow.lineage.model.v1beta1.Checksum
import nextflow.util.CacheHelper
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ChecksumTest extends Specification {

    def 'should create a checksum'() {
        given:
        def checksum = new Checksum(algorithm: 'sha1', value: '1234567890abcdef', mode: 'hex')

        expect:
        checksum.algorithm == 'sha1'
        checksum.value == '1234567890abcdef'
        checksum.mode == 'hex'
    }

    def 'should create a checksum with of factory method'() {
        given:
        def checksum1 = Checksum.of('1234567890abcdef','sha1', CacheHelper.HashMode.DEFAULT())

        expect:
        checksum1.algorithm == 'sha1'
        checksum1.value == '1234567890abcdef'
        checksum1.mode == 'standard'
    }

    def 'should create checksum with ofNextflow factory method'() {
        given:
        def checksum1 = Checksum.ofNextflow('1234567890abcdef')

        expect:
        checksum1.algorithm == 'nextflow'
        checksum1.value == CacheHelper.hasher('1234567890abcdef').hash().toString()
        checksum1.mode == 'standard'
    }

    @Unroll
    def 'should compute ofNextflow checksum with the mode it reports - #mode'() {
        given: 'a file and NXF_CACHE_MODE set to the given mode'
        Path path = Files.createTempFile('checksum', '.txt')
        path.text = 'Hello world'
        setDefaultHashMode(mode)

        when:
        def checksum = Checksum.ofNextflow(path)

        then: 'the reported mode is the configured one'
        checksum.algorithm == 'nextflow'
        checksum.mode == mode.toString().toLowerCase()

        and: 'the value is the hash for that mode, not for the standard mode'
        checksum.value == CacheHelper.hasher(path, mode).hash().toString()

        cleanup:
        setDefaultHashMode(null)
        Files.deleteIfExists(path)

        where:
        mode << CacheHelper.HashMode.values()
    }

    /**
     * Override the hash mode default, which {@link CacheHelper.HashMode} reads from
     * NXF_CACHE_MODE in a static initialiser and therefore cannot be changed by
     * setting the environment variable from a test.
     */
    private static void setDefaultHashMode(CacheHelper.HashMode mode) {
        final field = CacheHelper.HashMode.class.getDeclaredField('defaultValue')
        field.setAccessible(true)
        field.set(null, mode)
    }
}
