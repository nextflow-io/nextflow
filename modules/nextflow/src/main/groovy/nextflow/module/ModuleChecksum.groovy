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

import java.nio.charset.StandardCharsets
import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes
import java.security.MessageDigest
import java.util.concurrent.ConcurrentHashMap

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

    // a whole pipeline tree is costly to hash and an include resolves it several times per run
    private static final Map<Path,String> PIPELINE_CHECKSUMS = new ConcurrentHashMap<>()

    /**
     * Compute the SHA-256 checksum of a module directory
     */
    static String compute(Path moduleDir) {
        return SharedChecksum.computeDirectory(moduleDir)
    }

    /**
     * Compute the SHA-256 checksum of an installed pipeline, skipping {@code .git},
     * the {@code .module-info} file and nested folders with their own {@code .module-info}.
     * The result is memoized per real path until {@link #clearPipeline} is called.
     */
    static String computePipeline(Path dir) {
        return PIPELINE_CHECKSUMS.computeIfAbsent(realPath(dir)) { Path real -> computePipeline0(real) }
    }

    /**
     * Forget the memoized checksum of a pipeline, e.g. when it is installed or removed
     */
    static void clearPipeline(Path dir) {
        PIPELINE_CHECKSUMS.remove(realPath(dir))
    }

    private static Path realPath(Path dir) {
        return Files.exists(dir) ? dir.toRealPath() : dir.toAbsolutePath().normalize()
    }

    private static String computePipeline0(Path dir) {
        final files = new ArrayList<Path>()
        Files.walkFileTree(dir, new SimpleFileVisitor<Path>() {
            @Override
            FileVisitResult preVisitDirectory(Path d, BasicFileAttributes attrs) {
                if( d.fileName.toString() == '.git' )
                    return FileVisitResult.SKIP_SUBTREE
                if( d != dir && Files.exists(d.resolve(ModuleInfo.MODULE_INFO_FILE)) )
                    return FileVisitResult.SKIP_SUBTREE
                return FileVisitResult.CONTINUE
            }

            @Override
            FileVisitResult visitFile(Path file, BasicFileAttributes attrs) {
                if( attrs.isRegularFile() && file.fileName.toString() != ModuleInfo.MODULE_INFO_FILE )
                    files.add(file)
                return FileVisitResult.CONTINUE
            }
        })
        Collections.sort(files)

        final digest = MessageDigest.getInstance('SHA-256')
        for( final file : files ) {
            final content = Files.readAllBytes(file)
            // the path and length delimit each file
            digest.update("${dir.relativize(file)}\0${content.length}\0".toString().getBytes(StandardCharsets.UTF_8))
            digest.update(content)
        }
        return digest.digest().encodeHex().toString()
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
