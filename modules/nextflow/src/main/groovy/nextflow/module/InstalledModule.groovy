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
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import io.seqera.npr.api.schema.v1.ModuleMetadata
import org.yaml.snakeyaml.Yaml

import java.nio.file.Files
import java.nio.file.Path

/**
 * Represents a module installed in the local modules/ directory
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
@ToString(includeNames = true)
class InstalledModule {

    ModuleReference reference
    Path directory
    Path mainFile
    Path manifestFile
    Path checksumFile
    String installedVersion
    String expectedChecksum

    /**
     * Get the integrity status of this installed module
     *
     * @return ModuleIntegrity status
     */
    ModuleIntegrity getIntegrity() {
        // Check if main.nf exists
        if (!Files.exists(mainFile) || !Files.exists(manifestFile)) {
            return ModuleIntegrity.CORRUPTED
        }

        // Check if checksum file exists
        if (!Files.exists(checksumFile)) {
            return ModuleIntegrity.MISSING_CHECKSUM
        }

        try {
            // Compute actual checksum
            def actualChecksum = ModuleChecksum.compute(directory)

            // Compare with expected
            if (actualChecksum == expectedChecksum) {
                return ModuleIntegrity.VALID
            } else {
                log.debug("Actual: $actualChecksum, expected: $expectedChecksum")
                return ModuleIntegrity.MODIFIED
            }
        } catch (Exception e) {
            log.warn "Failed to compute checksum for module ${reference.nameWithoutPrefix}: ${e.message}"
            return ModuleIntegrity.CORRUPTED
        }
    }
}

/**
 * Module integrity status
 */
@CompileStatic
enum ModuleIntegrity {
    VALID,              // Checksum matches
    MODIFIED,           // Checksum mismatch (local changes)
    MISSING_CHECKSUM,   // No .checksum file
    CORRUPTED           // Missing required files
}
