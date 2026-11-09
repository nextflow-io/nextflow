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
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleStorage
import nextflow.module.spi.RemoteModuleResolver
import org.apache.commons.compress.archivers.tar.TarArchiveEntry
import org.apache.commons.compress.archivers.tar.TarArchiveOutputStream
import org.apache.commons.compress.compressors.gzip.GzipCompressorOutputStream
import org.eclipse.jgit.ignore.IgnoreNode
import org.eclipse.jgit.lib.FileMode
import org.eclipse.jgit.lib.Repository
import org.eclipse.jgit.lib.RepositoryBuilder

/**
 * Creates and extracts the {@code .tar.gz} bundle of a pipeline
 */
@Slf4j
@CompileStatic
class PipelineBundle {

    /** Relative paths to bundle, sorted. Symlinks are never followed */
    static List<Path> files(Path dir) {
        return new ArrayList<Path>(entries(dir).keySet())
    }

    /** Relative path to executable flag, sorted by path */
    private static Map<Path,Boolean> entries(Path dir) {
        final repo = openRepository(dir)
        final result = new TreeMap<Path,Boolean>()
        try {
            result.putAll(repo != null ? tracked(repo, dir) : walk(dir))
        }
        finally {
            repo?.close()
        }
        result.keySet().removeIf { Path p -> p.fileName.toString() == RemoteModuleResolver.PIPELINE_INFO_FILE }
        return result
    }

    /** Write the bundle and return its checksum as {@code sha256:<hex>} */
    static String create(Path dir, Path target) {
        try {
            if( target.parent )
                Files.createDirectories(target.parent)
            Files.newOutputStream(target).withCloseable { fos ->
                new GzipCompressorOutputStream(fos).withCloseable { gzos ->
                    new TarArchiveOutputStream(gzos).withCloseable { tos ->
                        tos.setLongFileMode(TarArchiveOutputStream.LONGFILE_POSIX)
                        tos.setBigNumberMode(TarArchiveOutputStream.BIGNUMBER_POSIX)
                        for( final it : entries(dir).entrySet() ) {
                            final file = dir.resolve(it.key)
                            // fixed mtime and no owner, so the same tree always gives the same checksum
                            final entry = new TarArchiveEntry(it.key.toString())
                            entry.setSize(Files.size(file))
                            entry.setMode(it.value ? 0755 : 0644)
                            entry.setModTime(0L)
                            tos.putArchiveEntry(entry)
                            Files.copy(file, tos)
                            tos.closeArchiveEntry()
                        }
                    }
                }
            }
        }
        catch( Exception e ) {
            Files.deleteIfExists(target)
            if( e instanceof AbortOperationException )
                throw e
            throw new AbortOperationException("Failed to create pipeline bundle", e)
        }
        return 'sha256:' + ModuleStorage.computeBundleChecksum(target)
    }

    static void extract(Path bundle, Path targetDir) {
        // pipelines legitimately ship modules/ and tests/ dirs
        ModuleStorage.extractTarGz(bundle, targetDir, false)
    }

    /** Open the Git repo whose work tree contains dir, or null when dir is not in a checkout. The caller closes it */
    static Repository openRepository(Path dir) {
        final gitDir = new RepositoryBuilder().findGitDir(dir.toFile()).gitDir
        if( gitDir == null )
            return null
        final repo = new RepositoryBuilder().setGitDir(gitDir).build()
        if( repo.bare ) {
            repo.close()
            return null
        }
        return repo
    }

    private static Map<Path,Boolean> tracked(Repository repo, Path dir) {
        final root = repo.workTree.toPath().toRealPath()
        final base = root.relativize(dir.toRealPath()).toString()
        final prefix = base.isEmpty() ? '' : base + '/'
        final result = new LinkedHashMap<Path,Boolean>()
        final skipped = new ArrayList<String>()
        final cache = repo.readDirCache()
        for( int i = 0; i < cache.entryCount; i++ ) {
            final entry = cache.getEntry(i)
            final name = entry.pathString
            if( !name.startsWith(prefix) )
                continue
            final mode = entry.fileMode
            if( mode == FileMode.SYMLINK || mode == FileMode.GITLINK ) {
                skipped.add(name)
                continue
            }
            final rel = dir.fileSystem.getPath(name.substring(prefix.length()))
            final file = dir.resolve(rel)
            if( Files.isRegularFile(file, LinkOption.NOFOLLOW_LINKS) )
                result.put(rel, mode == FileMode.EXECUTABLE_FILE)
        }
        if( skipped )
            log.warn "Not bundling tracked symlinks or submodules: ${skipped.join(', ')}"
        return result
    }

    // ponytail: root .gitignore only outside git; nested .gitignore files are ignored
    private static Map<Path,Boolean> walk(Path dir) {
        final ignore = new IgnoreNode()
        final gitignore = dir.resolve('.gitignore')
        if( Files.exists(gitignore) )
            gitignore.withInputStream { is -> ignore.parse(is) }
        final result = new LinkedHashMap<Path,Boolean>()
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
                if( attrs.isRegularFile() && !attrs.isSymbolicLink() && !ignored(dir, file, false) )
                    result.put(dir.relativize(file), Files.isExecutable(file))
                return FileVisitResult.CONTINUE
            }

            private boolean ignored(Path root, Path p, boolean isDir) {
                return ignore.isIgnored(root.relativize(p).toString(), isDir) == IgnoreNode.MatchResult.IGNORED
            }
        })
        return result
    }
}
