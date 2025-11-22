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
 *
 */

package io.seqera.wave.plugin.packer

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileTime

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.wave.plugin.ContainerLayer
import io.seqera.wave.plugin.util.DigestFunctions
import nextflow.file.FileHelper
import org.apache.commons.compress.archivers.tar.TarArchiveEntry
import org.apache.commons.compress.archivers.tar.TarArchiveOutputStream
import org.apache.commons.compress.compressors.gzip.GzipCompressorOutputStream

/**
 * Utility class to create container layer packages
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class Packer {

    /**
     * See {@link TarArchiveEntry#DEFAULT_DIR_MODE}
     */
    private static final int DIR_MODE = 040000;

    /**
     * See {@link TarArchiveEntry#DEFAULT_FILE_MODE}
     */
    private static final int FILE_MODE = 0100000;

    /**
     * Timestamp when {@link #preserveFileTimestamp} is false
     */
    private static final FileTime ZERO = FileTime.fromMillis(0)

    /**
     * Whenever the timestamp of compressed files should be preserved.
     * By default it is used the timestamp ZERO to guarantee reproducible builds
     */
    private boolean preserveFileTimestamp

    Packer withPreserveTimestamp(boolean value) {
        this.preserveFileTimestamp = value
        return this
    }

    def <T extends OutputStream> T makeTar(Path root, List<Path> files, T target) {
        final entries = new HashMap<String,Path>()
        for( Path it : files ) {
            final name = root.relativize(it).toString()
            entries.put(name, it)
        }
        return makeTar(entries, target)
    }

    def <T extends OutputStream> T makeTar(Map<String,Path> entries, T target) {
        try ( final archive = new TarArchiveOutputStream(target) ) {
            final sorted = new TreeSet<>(entries.keySet())
            for (String name : sorted ) {
                final targetPath = entries.get(name)
                final ftime = preserveFileTimestamp
                    ? Files.readAttributes(targetPath, BasicFileAttributes.class).lastModifiedTime()
                    : ZERO
                final entry = new TarArchiveEntry(targetPath, name)
                entry.setIds(0,0)
                entry.setGroupName("root")
                entry.setUserName("root")
                entry.setModTime(ftime)
                entry.setMode(getMode(targetPath))
                // file permissions
                archive.putArchiveEntry(entry)
                if( !targetPath.isDirectory()) {
                    Files.copy(targetPath, archive)
                }
                archive.closeArchiveEntry()
            }
            archive.finish()
        }

        return target
    }

    private int getMode(Path path) {
        final mode = path.isDirectory() ? DIR_MODE : FILE_MODE
        return mode + path.getPermissionsMode()
    }

    protected <T extends OutputStream> T makeGzip(InputStream source, T target) {
        try (final compressed = new GzipCompressorOutputStream(target)) {
            source.transferTo(compressed)
            compressed.flush()
        }
        return target
    }

    ContainerLayer layer(Path root, List<Path> files) {
        final Map<String,Path> entries = new HashMap<>()
        for( Path it : files )
            entries.put(root.relativize(it).toString(), it)
        return layer(entries)
    }

    ContainerLayer layer(Map<String,Path> entries) {
        final tar = makeTar(entries, new ByteArrayOutputStream()).toByteArray()
        final tarDigest = DigestFunctions.digest(tar)
        final gzipStream = new ByteArrayOutputStream()
        makeGzip(new ByteArrayInputStream(tar), gzipStream); gzipStream.close()
        final gzipBytes = gzipStream.toByteArray()
        final gzipSize = gzipBytes.length
        final gzipDigest = DigestFunctions.digest(gzipBytes)
        final data = 'data:' + gzipBytes.encodeBase64()

        return new ContainerLayer(
                location: data,
                tarDigest: tarDigest,
                gzipSize: gzipSize,
                gzipDigest: gzipDigest )
    }

    ContainerLayer createContainerPack(Path root, String targetName) {
        final opts = [type: 'any', hidden: true, relative: false]
        final files = new ArrayList(100)
        FileHelper.visitFiles(opts, root, '**') { files.add(it) }
        return createContainerPack0(root, files, targetName)
    }

    ContainerLayer createContainerPack0(Path root, List<Path> files, String targetName=null) {
        final name = targetName ?: root.getName()
        // create the tar file
        final tarFilePath = root.resolveSibling(name + '.tar')
        makeTar(root, files, new FileOutputStream(tarFilePath.toFile()))
            .close()

        // compute the tar digest
        final tarDigest = DigestFunctions.digest(tarFilePath)

        // target gzipped file
        final gzipFilePath = root.resolveSibling(name + '.tar.gz')
        makeGzip(tarFilePath.newInputStream(), new FileOutputStream(gzipFilePath.toFile()))
                .close()
        final gzipSize = Files.size(gzipFilePath)

        // compute gzip digest
        final gzipDigest = DigestFunctions.digest(gzipFilePath)

        // remove tar file
        tarFilePath.delete()

        // finally return the expected  
        return new ContainerLayer(
                location: gzipFilePath.toUri().toString(),
                tarDigest: tarDigest,
                gzipSize: gzipSize,
                gzipDigest: gzipDigest )
    }

}
