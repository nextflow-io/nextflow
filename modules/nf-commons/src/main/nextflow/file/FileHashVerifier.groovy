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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Implements file hash verification functionality
 *
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
@Slf4j
@CompileStatic
class FileHashVerifier {

    private static final Map<String, String> SUPPORTED_ALGORITHMS = [
        'md5': 'MD5',
        'sha1': 'SHA-1',
        'sha256': 'SHA-256',
        'sha384': 'SHA-384',
        'sha512': 'SHA-512'
    ]

    /**
     * Verify a file's hash matches the expected value
     *
     * @param file The file to verify
     * @param knownHash The expected hash in format 'algorithm:hash'
     * @return The file if verification succeeds
     * @throws IllegalArgumentException if hash format is invalid or verification fails
     */
    static Path verifyHash(Path file, String knownHash) {
        if (!knownHash) {
            return file
        }

        if (Files.isDirectory(file)) {
            throw new IllegalArgumentException("Cannot verify hash of a directory: ${file.toString()}")
        }

        def parts = knownHash.split(':', 2)
        if (parts.length != 2) {
            throw new IllegalArgumentException("Invalid hash format. Expected 'algorithm:hash', got: ${knownHash}")
        }

        String algorithm = parts[0].toLowerCase()
        String expectedHash = parts[1].toLowerCase()

        if (!SUPPORTED_ALGORITHMS.containsKey(algorithm)) {
            throw new IllegalArgumentException("Unsupported hash algorithm: ${algorithm}. Supported algorithms: ${SUPPORTED_ALGORITHMS.keySet().join(', ')}")
        }

        String actualHash = calculateHash(file, SUPPORTED_ALGORITHMS[algorithm])
        if (actualHash != expectedHash) {
            throw new IllegalArgumentException(
                "Hash verification failed for file: ${file.toString()}\n" +
                "Expected: ${expectedHash}\n" +
                "Actual: ${actualHash}"
            )
        }

        return file
    }

    /**
     * Calculate hash for a file using the specified algorithm
     *
     * @param file The file to hash
     * @param algorithm The hash algorithm to use
     * @return The calculated hash as a hex string
     */
    private static String calculateHash(Path file, String algorithm) {
        def digest = MessageDigest.getInstance(algorithm)
        Files.newInputStream(file).withStream { input ->
            byte[] buffer = new byte[8192]
            int read
            while ((read = input.read(buffer)) != -1) {
                digest.update(buffer, 0, read)
            }
        }
        return digest.digest().encodeHex().toString()
    }
}
