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

package io.seqera.tower.plugin.fs

import nextflow.file.FileHelper

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.file.FileSystemPathFactory

/**
 * PF4J extension that registers the {@code seqera://} URI scheme with Nextflow's file helper,
 * allowing pipeline scripts to use {@code file('seqera://org/ws/datasets/name')} transparently.
 *
 * Registered as a PF4J extension point via {@code extensionPoints} in {@code build.gradle}.
 *
 * @author Seqera Labs
 */
@Slf4j
@CompileStatic
class SeqeraPathFactory extends FileSystemPathFactory {


    @Override
    Path parseUri(String str) {
        return str?.startsWith(SeqeraPath.PROTOCOL) ? create(str) : null
    }

    @Override
    String toUriString(Path path) {
        if (path instanceof SeqeraPath)
            return path.toUri().toString()
        return null
    }

    @Override
    String getBashLib(Path target) {
        // No bash-level staging for seqera:// — handled via NIO newInputStream/copy
        return null
    }

    @Override
    String getUploadCmd(String source, Path target) {
        // No bash upload command — upload handled via NIO newOutputStream
        return null
    }

    static SeqeraPath create(String path) {
        final uri = SeqeraPath.asUri(path)
        return (SeqeraPath) FileHelper.getOrCreateFileSystemFor(uri).provider().getPath(uri)
    }

}
