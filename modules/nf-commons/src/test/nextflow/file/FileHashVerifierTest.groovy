/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.file

import java.nio.file.Files
import java.nio.file.Path
import java.security.MessageDigest
import java.io.IOException
import java.nio.file.NoSuchFileException

import spock.lang.Specification
import spock.lang.Unroll

/**
 * Tests for FileHashVerifier
 *
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
class FileHashVerifierTest extends Specification {

    def 'should return file when no hash is provided'() {
        given:
        def file = Files.createTempFile('test', '.txt')
        file.text = 'test content'

        when:
        def result = FileHashVerifier.verifyHash(file, null)

        then:
        result == file

        cleanup:
        file?.deleteDir()
    }

    def 'should verify file hash with MD5'() {
        given:
        def content = 'test content'
        def file = Files.createTempFile('test', '.txt')
        file.text = content
        def md5Hash = MessageDigest.getInstance("MD5").digest(content.bytes).encodeHex().toString()

        when:
        def result = FileHashVerifier.verifyHash(file, "md5:${md5Hash}")

        then:
        result == file

        cleanup:
        file?.deleteDir()
    }

    def 'should verify file hash with SHA-256'() {
        given:
        def content = 'test content'
        def file = Files.createTempFile('test', '.txt')
        file.text = content
        def sha256Hash = MessageDigest.getInstance("SHA-256").digest(content.bytes).encodeHex().toString()

        when:
        def result = FileHashVerifier.verifyHash(file, "sha256:${sha256Hash}")

        then:
        result == file

        cleanup:
        file?.deleteDir()
    }

    def 'should throw exception on hash mismatch'() {
        given:
        def file = Files.createTempFile('test', '.txt')
        file.text = 'test content'

        when:
        FileHashVerifier.verifyHash(file, "md5:invalidhash")

        then:
        def e = thrown(IllegalArgumentException)
        e.message.startsWith("Hash verification failed for file:")

        cleanup:
        file?.deleteDir()
    }

    def 'should throw exception on invalid hash format'() {
        given:
        def file = Files.createTempFile('test', '.txt')
        file.text = 'test content'

        when:
        FileHashVerifier.verifyHash(file, "invalid:hash")

        then:
        def e = thrown(IllegalArgumentException)
        e.message.startsWith("Unsupported hash algorithm: invalid")

        cleanup:
        file?.deleteDir()
    }

    @Unroll
    def 'should support different hash algorithms: #algorithm'() {
        given:
        def content = 'test content'
        def file = Files.createTempFile('test', '.txt')
        file.text = content
        def hash = MessageDigest.getInstance(javaAlgorithm).digest(content.bytes).encodeHex().toString()

        when:
        def result = FileHashVerifier.verifyHash(file, "${algorithm}:${hash}")

        then:
        result == file

        cleanup:
        file?.deleteDir()

        where:
        algorithm | javaAlgorithm
        'md5'     | 'MD5'
        'sha1'    | 'SHA-1'
        'sha256'  | 'SHA-256'
        'sha384'  | 'SHA-384'
        'sha512'  | 'SHA-512'
    }

    def 'should handle large files'() {
        given:
        def file = Files.createTempFile('test', '.txt')
        def content = new StringBuilder()
        10000.times { content.append('test content\n') }
        file.text = content.toString()
        def sha256Hash = MessageDigest.getInstance("SHA-256").digest(content.toString().bytes).encodeHex().toString()

        when:
        def result = FileHashVerifier.verifyHash(file, "sha256:${sha256Hash}")

        then:
        result == file

        cleanup:
        file?.deleteDir()
    }

    def 'should throw exception when verifying hash of a directory'() {
        given:
        def dir = Files.createTempDirectory('test-dir')
        def content = 'test content'
        def file1 = dir.resolve('file1.txt')
        def file2 = dir.resolve('file2.txt')
        file1.text = content
        file2.text = content
        def sha256Hash = MessageDigest.getInstance("SHA-256").digest(content.bytes).encodeHex().toString()

        when:
        FileHashVerifier.verifyHash(dir, "sha256:${sha256Hash}")

        then:
        def e = thrown(IllegalArgumentException)
        e.message == "Cannot verify hash of a directory: ${dir.toString()}"

        cleanup:
        dir?.deleteDir()
    }

    def 'should throw exception for non-existent file'() {
        given:
        def file = Files.createTempFile('test', '.txt')
        Files.delete(file)
        def sha256Hash = MessageDigest.getInstance("SHA-256").digest('any content'.bytes).encodeHex().toString()

        when:
        FileHashVerifier.verifyHash(file, "sha256:${sha256Hash}")

        then:
        def e = thrown(NoSuchFileException)
        e.file == file.toString()

        cleanup:
        file?.deleteDir()
    }

    def 'should verify hash of empty file'() {
        given:
        def file = Files.createTempFile('test', '.txt')
        def emptyHash = MessageDigest.getInstance("SHA-256").digest(''.bytes).encodeHex().toString()

        when:
        def result = FileHashVerifier.verifyHash(file, "sha256:${emptyHash}")

        then:
        result == file

        cleanup:
        file?.deleteDir()
    }
}
