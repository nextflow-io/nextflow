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
import java.security.MessageDigest

import groovy.transform.CompileStatic
import nextflow.module.ModuleInfo
import nextflow.module.spi.RemoteModuleResolver

/**
 * Metadata of an installed pipeline, stored in its {@code .pipeline-info} file
 */
@CompileStatic
class PipelineInfo {

    private static final List<String> FIELDS = ['name', 'version', 'workflowPath', 'repositoryUrl', 'revision', 'commitId', 'checksum', 'registryUrl']

    String name
    String version
    String workflowPath
    String repositoryUrl
    String revision
    String commitId
    String checksum
    String registryUrl

    static PipelineInfo load(Path dir) {
        final file = dir.resolve(RemoteModuleResolver.PIPELINE_INFO_FILE)
        if( !Files.exists(file) )
            return null
        final props = new Properties()
        file.withInputStream { is -> props.load(is) }
        final result = new PipelineInfo()
        for( final f : FIELDS )
            result.setProperty(f, props.getProperty(f))
        return result
    }

    void save(Path dir) {
        final props = new Properties()
        for( final f : FIELDS ) {
            final value = getProperty(f)
            if( value != null )
                props.setProperty(f, value.toString())
        }
        dir.resolve(RemoteModuleResolver.PIPELINE_INFO_FILE).withOutputStream { os -> props.store(os, null) }
    }

    boolean isModified(Path dir) {
        return checksum(dir) != this.checksum
    }

    /**
     * Same digest as {@code ModuleChecksum.computeDirectory}, skipping {@code .git}
     * and nested modules/pipelines (they carry their own info file).
     */
    static String checksum(Path dir) {
        final files = new ArrayList<Path>()
        Files.walkFileTree(dir, new SimpleFileVisitor<Path>() {
            @Override
            FileVisitResult preVisitDirectory(Path d, BasicFileAttributes attrs) {
                if( d.fileName.toString() == '.git' )
                    return FileVisitResult.SKIP_SUBTREE
                if( d != dir && (Files.exists(d.resolve(ModuleInfo.MODULE_INFO_FILE)) || Files.exists(d.resolve(RemoteModuleResolver.PIPELINE_INFO_FILE))) )
                    return FileVisitResult.SKIP_SUBTREE
                return FileVisitResult.CONTINUE
            }

            @Override
            FileVisitResult visitFile(Path file, BasicFileAttributes attrs) {
                if( attrs.isRegularFile() && file.fileName.toString() != RemoteModuleResolver.PIPELINE_INFO_FILE )
                    files.add(file)
                return FileVisitResult.CONTINUE
            }
        })
        Collections.sort(files)

        final digest = MessageDigest.getInstance('SHA-256')
        for( final file : files ) {
            digest.update(dir.relativize(file).toString().getBytes())
            digest.update(Files.readAllBytes(file))
        }
        return digest.digest().encodeHex().toString()
    }
}
