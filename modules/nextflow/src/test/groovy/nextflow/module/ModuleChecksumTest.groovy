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

package nextflow.module

import java.nio.file.Files
import java.nio.file.Path

import spock.lang.Specification

/**
 * Test suite for ModuleChecksum
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleChecksumTest extends Specification {

    Path tempDir

    def setup() {
        tempDir = Files.createTempDirectory('nf-checksum-test-')
    }

    def cleanup() {
        tempDir?.deleteDir()
    }

    def 'should compute checksum for directory'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        // Create test files
        moduleDir.resolve('main.nf').text = 'process TEST { }'
        moduleDir.resolve('meta.yml').text = 'name: test\nversion: 1.0.0'
        moduleDir.resolve('README.md').text = '# Test Module'

        when:
        def checksum = ModuleChecksum.compute(moduleDir)

        then:
        checksum != null
        checksum.length() == 64  // SHA-256 produces 64 hex characters
        checksum ==~ /^[a-f0-9]{64}$/
    }

    def 'should produce consistent checksums for same content'() {
        given:
        def moduleDir1 = tempDir.resolve('module1')
        def moduleDir2 = tempDir.resolve('module2')
        Files.createDirectories(moduleDir1)
        Files.createDirectories(moduleDir2)

        // Create identical content in both directories
        ['main.nf', 'meta.yml', 'README.md'].each { filename ->
            moduleDir1.resolve(filename).text = "content of ${filename}"
            moduleDir2.resolve(filename).text = "content of ${filename}"
        }

        when:
        def checksum1 = ModuleChecksum.compute(moduleDir1)
        def checksum2 = ModuleChecksum.compute(moduleDir2)

        then:
        checksum1 == checksum2
    }

    def 'should produce different checksums for different content'() {
        given:
        def moduleDir1 = tempDir.resolve('module1')
        def moduleDir2 = tempDir.resolve('module2')
        Files.createDirectories(moduleDir1)
        Files.createDirectories(moduleDir2)

        // Create different content
        moduleDir1.resolve('main.nf').text = 'process TEST1 { }'
        moduleDir2.resolve('main.nf').text = 'process TEST2 { }'

        when:
        def checksum1 = ModuleChecksum.compute(moduleDir1)
        def checksum2 = ModuleChecksum.compute(moduleDir2)

        then:
        checksum1 != checksum2
    }

    def 'should exclude .checksum file from computation'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        moduleDir.resolve('main.nf').text = 'process TEST { }'

        // Compute initial checksum
        def checksum1 = ModuleChecksum.compute(moduleDir)

        // Add .checksum file
        moduleDir.resolve('.checksum').text = 'some-checksum-value'

        // Compute checksum again
        def checksum2 = ModuleChecksum.compute(moduleDir)

        expect:
        checksum1 == checksum2  // Should be the same, .checksum is ignored
    }

    def 'should include subdirectories in checksum'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        moduleDir.resolve('main.nf').text = 'process TEST { }'

        // Compute checksum without subdirectory
        def checksum1 = ModuleChecksum.compute(moduleDir)

        // Add subdirectory with file
        def subDir = moduleDir.resolve('templates')
        Files.createDirectories(subDir)
        subDir.resolve('script.sh').text = '#!/bin/bash\necho "test"'

        // Compute checksum with subdirectory
        def checksum2 = ModuleChecksum.compute(moduleDir)

        expect:
        checksum1 != checksum2  // Checksums should differ
    }

    def 'should save checksum to .checksum file'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)
        def checksumValue = 'abc123def456'

        when:
        ModuleChecksum.save(moduleDir, checksumValue)

        then:
        def checksumFile = moduleDir.resolve('.checksum')
        Files.exists(checksumFile)
        checksumFile.text.trim() == checksumValue
    }

    def 'should load checksum from .checksum file'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)
        def checksumFile = moduleDir.resolve('.checksum')
        checksumFile.text = 'abc123def456'

        when:
        def checksum = ModuleChecksum.load(moduleDir)

        then:
        checksum == 'abc123def456'
    }

    def 'should return null when loading non-existent checksum file'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        when:
        def checksum = ModuleChecksum.load(moduleDir)

        then:
        checksum == null
    }

    def 'should handle empty directory'() {
        given:
        def moduleDir = tempDir.resolve('empty-module')
        Files.createDirectories(moduleDir)

        when:
        def checksum = ModuleChecksum.compute(moduleDir)

        then:
        checksum != null
        checksum.length() == 64
    }

    def 'should handle files with special characters in names'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        // Create files with special characters (but valid on filesystem)
        moduleDir.resolve('main.nf').text = 'process TEST { }'
        moduleDir.resolve('file with spaces.txt').text = 'content'
        moduleDir.resolve('file-with-dashes.txt').text = 'content'

        when:
        def checksum = ModuleChecksum.compute(moduleDir)

        then:
        checksum != null
        checksum.length() == 64
    }

    def 'should handle binary files'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        // Create text and binary files
        moduleDir.resolve('main.nf').text = 'process TEST { }'
        def binaryFile = moduleDir.resolve('data.bin')
        binaryFile.bytes = [0x00, 0x01, 0x02, 0xFF] as byte[]

        when:
        def checksum = ModuleChecksum.compute(moduleDir)

        then:
        checksum != null
        checksum.length() == 64
    }

    def 'should sort files consistently for checksum computation'() {
        given:
        def moduleDir = tempDir.resolve('module')
        Files.createDirectories(moduleDir)

        // Create files in arbitrary order
        moduleDir.resolve('zzz.nf').text = 'content'
        moduleDir.resolve('aaa.nf').text = 'content'
        moduleDir.resolve('mmm.nf').text = 'content'

        def checksum1 = ModuleChecksum.compute(moduleDir)

        // Create another directory with files in different order
        def moduleDir2 = tempDir.resolve('module2')
        Files.createDirectories(moduleDir2)
        moduleDir2.resolve('aaa.nf').text = 'content'
        moduleDir2.resolve('mmm.nf').text = 'content'
        moduleDir2.resolve('zzz.nf').text = 'content'

        def checksum2 = ModuleChecksum.compute(moduleDir2)

        expect:
        checksum1 == checksum2  // Order shouldn't matter
    }

    // Tests for computeFile() method

    def 'should compute checksum for a single file with default algorithm'() {
        given:
        def testFile = tempDir.resolve('test.txt')
        testFile.text = 'Hello, World!'

        when:
        def checksum = ModuleChecksum.computeFile(testFile)

        then:
        checksum != null
        checksum.length() == 64  // SHA-256 produces 64 hex characters
        checksum ==~ /^[a-f0-9]{64}$/
    }

    def 'should produce consistent checksums for same file content'() {
        given:
        def file1 = tempDir.resolve('file1.txt')
        def file2 = tempDir.resolve('file2.txt')
        def content = 'Same content in both files'
        file1.text = content
        file2.text = content

        when:
        def checksum1 = ModuleChecksum.computeFile(file1)
        def checksum2 = ModuleChecksum.computeFile(file2)

        then:
        checksum1 == checksum2
    }

    def 'should produce different checksums for different file content'() {
        given:
        def file1 = tempDir.resolve('file1.txt')
        def file2 = tempDir.resolve('file2.txt')
        file1.text = 'Content A'
        file2.text = 'Content B'

        when:
        def checksum1 = ModuleChecksum.computeFile(file1)
        def checksum2 = ModuleChecksum.computeFile(file2)

        then:
        checksum1 != checksum2
    }

    def 'should compute checksum for empty file'() {
        given:
        def emptyFile = tempDir.resolve('empty.txt')
        emptyFile.text = ''

        when:
        def checksum = ModuleChecksum.computeFile(emptyFile)

        then:
        checksum != null
        checksum.length() == 64
        // SHA-256 of empty string
        checksum == 'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855'
    }

    def 'should compute checksum for binary file'() {
        given:
        def binaryFile = tempDir.resolve('binary.dat')
        binaryFile.bytes = [0x00, 0xFF, 0x42, 0xAB, 0xCD, 0xEF] as byte[]

        when:
        def checksum = ModuleChecksum.computeFile(binaryFile)

        then:
        checksum != null
        checksum.length() == 64
        checksum ==~ /^[a-f0-9]{64}$/
    }

    def 'should support different hash algorithms'() {
        given:
        def testFile = tempDir.resolve('test.txt')
        testFile.text = 'Test content for different algorithms'

        when:
        def sha256 = ModuleChecksum.computeFile(testFile, 'SHA-256')
        def sha512 = ModuleChecksum.computeFile(testFile, 'SHA-512')

        then:
        sha256 != null
        sha512 != null
        sha256.length() == 64   // SHA-256: 256 bits = 64 hex chars
        sha512.length() == 128  // SHA-512: 512 bits = 128 hex chars
        sha256 != sha512
    }

    def 'should handle case-insensitive algorithm names'() {
        given:
        def testFile = tempDir.resolve('test.txt')
        testFile.text = 'Test content'

        when:
        def checksum1 = ModuleChecksum.computeFile(testFile, 'sha-256')
        def checksum2 = ModuleChecksum.computeFile(testFile, 'SHA-256')
        def checksum3 = ModuleChecksum.computeFile(testFile, 'Sha-256')

        then:
        checksum1 == checksum2
        checksum2 == checksum3
    }

    def 'should throw exception for non-existent file'() {
        given:
        def nonExistentFile = tempDir.resolve('does-not-exist.txt')

        when:
        ModuleChecksum.computeFile(nonExistentFile)

        then:
        thrown(IllegalArgumentException)
    }

    def 'should throw exception for directory instead of file'() {
        given:
        def directory = tempDir.resolve('subdir')
        Files.createDirectories(directory)

        when:
        ModuleChecksum.computeFile(directory)

        then:
        thrown(IllegalArgumentException)
    }

    def 'should compute known SHA-256 checksum correctly'() {
        given:
        def testFile = tempDir.resolve('known.txt')
        testFile.text = 'abc'

        when:
        def checksum = ModuleChecksum.computeFile(testFile)

        then:
        // Known SHA-256 hash of "abc"
        checksum == 'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad'
    }

    def 'should handle large file checksum computation'() {
        given:
        def largeFile = tempDir.resolve('large.txt')
        // Create a file with ~1MB of data
        def content = 'x' * 1024 * 1024
        largeFile.text = content

        when:
        def checksum = ModuleChecksum.computeFile(largeFile)

        then:
        checksum != null
        checksum.length() == 64
    }
}
