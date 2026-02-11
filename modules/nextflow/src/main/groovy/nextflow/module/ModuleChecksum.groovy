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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

import java.nio.file.Files
import java.nio.file.Path
import java.security.MessageDigest

/**
 * Utility class for computing SHA-256 checksums of module directories
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class ModuleChecksum {

    public static final String CHECKSUM_ALGORITHM = "SHA-256"
    public static final String CHECKSUM_FILE = ".checksum"

    /**
     * Compute the SHA-256 checksum of a module directory
     *
     * @param moduleDir The module directory path
     * @return The hex-encoded SHA-256 checksum
     */
    static String compute(Path moduleDir) {
        if (!Files.exists(moduleDir) || !Files.isDirectory(moduleDir)) {
            throw new IllegalArgumentException("Module directory does not exist or is not a directory: ${moduleDir}")
        }

        try {
            def digest = MessageDigest.getInstance(CHECKSUM_ALGORITHM)

            // Collect all files in sorted order for consistent checksums
            List<Path> files = []
            Files.walk(moduleDir)
                .filter { Path path -> Files.isRegularFile(path) }
                .filter { Path path -> !path.fileName.toString().equals(CHECKSUM_FILE) }
                .sorted()
                .each { Path path -> files.add(path) }

            // Compute checksum over all file contents
            for (Path file : files) {
                // Include relative path in checksum for directory structure integrity
                def relativePath = moduleDir.relativize(file).toString()
                digest.update(relativePath.bytes)

                // Include file contents
                def bytes = Files.readAllBytes(file)
                digest.update(bytes)
            }

            def hashBytes = digest.digest()
            return bytesToHex(hashBytes)
        }
        catch (Exception e) {
            log.error("Failed to compute checksum for module directory: ${moduleDir}", e)
            throw new RuntimeException("Failed to compute module checksum", e)
        }
    }

    /**
     * Save a checksum to the .checksum file in the module directory
     *
     * @param moduleDir The module directory path
     * @param checksum The checksum to save
     */
    static void save(Path moduleDir, String checksum) {
        def checksumFile = moduleDir.resolve(CHECKSUM_FILE)
        Files.writeString(checksumFile, checksum)
    }

    /**
     * Load a checksum from the .checksum file in the module directory
     *
     * @param moduleDir The module directory path
     * @return The checksum, or null if file doesn't exist
     */
    static String load(Path moduleDir) {
        def checksumFile = moduleDir.resolve(CHECKSUM_FILE)
        if (!Files.exists(checksumFile)) {
            return null
        }
        return checksumFile.text
    }

    /**
     * Verify that a module directory matches the expected checksum
     *
     * @param moduleDir The module directory path
     * @param expectedChecksum The expected checksum
     * @return true if checksums match, false otherwise
     */
    static boolean verify(Path moduleDir, String expectedChecksum) {
        def actualChecksum = compute(moduleDir)
        return actualChecksum == expectedChecksum
    }

    /**
     * Compute the checksum of a single file
     *
     * @param file The file path
     * @param type checksum algorithm (sha-256 if not provided)
     * @return The hex-encoded checksum
     */
    static String computeFile(Path file, String type = CHECKSUM_ALGORITHM) {
        if (!Files.exists(file) || !Files.isRegularFile(file)) {
            throw new IllegalArgumentException("File does not exist or is not a regular file: ${file}")
        }

        try {
            final digest = MessageDigest.getInstance(type)
            final bytes = Files.readAllBytes(file)
            digest.update(bytes)
            final hashBytes = digest.digest()
            return bytesToHex(hashBytes)
        }
        catch (Exception e) {
            log.error("Failed to compute checksum for file: ${file}", e)
            throw new RuntimeException("Failed to compute file checksum", e)
        }
    }

    /**
     * Convert byte array to hex string
     *
     * @param bytes The byte array
     * @return Hex-encoded string
     */
    private static String bytesToHex(byte[] bytes) {
        def hexChars = new char[bytes.length * 2]
        for (int i = 0; i < bytes.length; i++) {
            int v = bytes[i] & 0xFF
            hexChars[i * 2] = Character.forDigit(v >>> 4, 16)
            hexChars[i * 2 + 1] = Character.forDigit(v & 0x0F, 16)
        }
        return new String(hexChars)
    }
}
