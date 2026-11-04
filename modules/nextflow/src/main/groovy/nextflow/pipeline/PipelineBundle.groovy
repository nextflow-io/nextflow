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

package nextflow.pipeline

import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes

import groovy.transform.CompileStatic
import nextflow.module.ModuleStorage
import nextflow.module.spi.RemoteModuleResolver
import org.apache.commons.compress.archivers.tar.TarArchiveEntry
import org.apache.commons.compress.archivers.tar.TarArchiveOutputStream
import org.apache.commons.compress.compressors.gzip.GzipCompressorOutputStream
import org.eclipse.jgit.ignore.IgnoreNode
import org.eclipse.jgit.lib.Repository
import org.eclipse.jgit.lib.RepositoryBuilder

/**
 * Creates and extracts the {@code .tar.gz} bundle of a pipeline
 */
@CompileStatic
class PipelineBundle {

    /** Relative paths to bundle, sorted */
    static List<Path> files(Path dir) {
        final result = isGit(dir) ? tracked(dir) : walk(dir)
        result.removeIf { Path p -> p.fileName.toString() == RemoteModuleResolver.PIPELINE_INFO_FILE }
        Collections.sort(result)
        return result
    }

    /** Write the bundle and return its checksum as {@code sha256:<hex>} */
    static String create(Path dir, Path target) {
        if( target.parent )
            Files.createDirectories(target.parent)
        Files.newOutputStream(target).withCloseable { fos ->
            new GzipCompressorOutputStream(fos).withCloseable { gzos ->
                new TarArchiveOutputStream(gzos).withCloseable { tos ->
                    for( final rel : files(dir) ) {
                        final file = dir.resolve(rel)
                        final entry = new TarArchiveEntry(file.toFile(), rel.toString())
                        entry.setSize(Files.size(file))
                        tos.putArchiveEntry(entry)
                        Files.copy(file, tos)
                        tos.closeArchiveEntry()
                    }
                }
            }
        }
        return 'sha256:' + ModuleStorage.computeBundleChecksum(target)
    }

    static void extract(Path bundle, Path targetDir) {
        // pipelines legitimately ship modules/ and tests/ dirs
        ModuleStorage.extractTarGz(bundle, targetDir, false)
    }

    private static boolean isGit(Path dir) {
        // only a repo rooted at dir: index paths are relative to the repo root
        return Files.exists(dir.resolve('.git'))
    }

    private static List<Path> tracked(Path dir) {
        final builder = new RepositoryBuilder().setWorkTree(dir.toFile()).findGitDir(dir.toFile())
        final result = new ArrayList<Path>()
        try (Repository repo = builder.build()) {
            final cache = repo.readDirCache()
            for( int i = 0; i < cache.entryCount; i++ ) {
                final rel = dir.fileSystem.getPath(cache.getEntry(i).pathString)
                if( Files.isRegularFile(dir.resolve(rel)) )
                    result.add(rel)
            }
        }
        return result
    }

    // ponytail: root .gitignore only outside git; nested .gitignore files are ignored
    private static List<Path> walk(Path dir) {
        final ignore = new IgnoreNode()
        final gitignore = dir.resolve('.gitignore')
        if( Files.exists(gitignore) )
            gitignore.withInputStream { is -> ignore.parse(is) }
        final result = new ArrayList<Path>()
        Files.walkFileTree(dir, new SimpleFileVisitor<Path>() {
            @Override
            FileVisitResult preVisitDirectory(Path d, BasicFileAttributes attrs) {
                if( d == dir )
                    return FileVisitResult.CONTINUE
                if( d.fileName.toString() == '.git' || PipelineInfo.isNested(d) )
                    return FileVisitResult.SKIP_SUBTREE
                return ignored(dir, d, true) ? FileVisitResult.SKIP_SUBTREE : FileVisitResult.CONTINUE
            }

            @Override
            FileVisitResult visitFile(Path file, BasicFileAttributes attrs) {
                if( attrs.isRegularFile() && !ignored(dir, file, false) )
                    result.add(dir.relativize(file))
                return FileVisitResult.CONTINUE
            }

            private boolean ignored(Path root, Path p, boolean isDir) {
                return ignore.isIgnored(root.relativize(p).toString(), isDir) == IgnoreNode.MatchResult.IGNORED
            }
        })
        return result
    }
}
