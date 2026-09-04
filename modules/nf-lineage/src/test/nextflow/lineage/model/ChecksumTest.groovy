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

import nextflow.lineage.model.v1beta1.Checksum
import nextflow.util.CacheHelper
import spock.lang.Specification

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

    def 'should hash a path in the mode it reports'() {
        given:
        def folder = Files.createTempDirectory('test')
        def file = folder.resolve('foo.txt'); file.text = 'Hello world'

        when:
        def checksum = Checksum.ofNextflow(file, MODE)

        then:
        checksum.algorithm == 'nextflow'
        checksum.mode == EXPECTED
        checksum.value == CacheHelper.hasher(file, MODE).hash().toString()

        cleanup:
        folder?.deleteDir()

        where:
        MODE                          | EXPECTED
        CacheHelper.HashMode.STANDARD | 'standard'
        CacheHelper.HashMode.LENIENT  | 'lenient'
        CacheHelper.HashMode.DEEP     | 'deep'
        CacheHelper.HashMode.SHA256   | 'sha256'
    }

    def 'should give identical files the same deep checksum'() {
        given:
        def folder = Files.createTempDirectory('test')
        def one = folder.resolve('one.txt'); one.text = 'Hello world'
        def two = folder.resolve('two.txt'); two.text = 'Hello world'

        expect:
        Checksum.ofNextflow(one, CacheHelper.HashMode.DEEP).value ==
            Checksum.ofNextflow(two, CacheHelper.HashMode.DEEP).value

        cleanup:
        folder?.deleteDir()
    }
}
