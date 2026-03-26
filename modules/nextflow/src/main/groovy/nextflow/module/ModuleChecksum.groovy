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
import io.seqera.npr.client.ModuleChecksum as SharedChecksum

import java.nio.file.Path

/**
 * Nextflow-specific checksum utility that delegates to the shared
 * {@link io.seqera.npr.client.ModuleChecksum} implementation and adds
 * convenience methods for .module-info persistence.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
class ModuleChecksum {

    public static final String CHECKSUM_ALGORITHM = SharedChecksum.ALGORITHM

    /**
     * Compute the SHA-256 checksum of a module directory
     */
    static String compute(Path moduleDir) {
        return SharedChecksum.computeDirectory(moduleDir)
    }

    /**
     * Save a checksum to the .module-info file
     */
    static void save(Path moduleDir, String checksum) {
        ModuleInfo.save(moduleDir, 'checksum', checksum)
    }

    /**
     * Load a checksum from the .module-info file
     */
    static String load(Path moduleDir) {
        ModuleInfo.load(moduleDir, 'checksum')
    }

    /**
     * Verify that a module directory matches the expected checksum
     */
    static boolean verify(Path moduleDir, String expectedChecksum) {
        def actualChecksum = compute(moduleDir)
        return actualChecksum == expectedChecksum
    }

    /**
     * Compute the checksum of a single file
     */
    static String computeFile(Path file, String type = CHECKSUM_ALGORITHM) {
        return SharedChecksum.computeFile(file, type)
    }
}
