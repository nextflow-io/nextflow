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

import spock.lang.Specification

/**
 * Guards that the generic NIO layer does not reach into resource-type-specific packages.
 *
 * {@link SeqeraPath}, {@link SeqeraFileSystem}, {@link SeqeraFileAttributes} must not
 * depend on {@code dataset/}, {@code datalink/}, or {@code fs/handler/}. Dispatch goes
 * through the {@link ResourceTypeHandler} interface.
 *
 * {@link SeqeraFileSystemProvider} is the dispatch point: it wires handlers and routes
 * calls to them, so it is *expected* to import the handler packages. The guard only
 * applies to the generic classes above.
 */
class ResourceTypeAbstractionTest extends Specification {

    static final Class[] GENERIC_CLASSES = [SeqeraPath, SeqeraFileSystem, SeqeraFileAttributes]

    private static File srcRoot() {
        // Gradle test cwd may be the plugin module dir or the repo root.
        final candidates = [
                'src/main/io/seqera/tower/plugin/fs',
                'plugins/nf-tower/src/main/io/seqera/tower/plugin/fs'
        ]
        for (String c : candidates) {
            final f = new File(c)
            if (f.isDirectory()) return f
        }
        throw new IllegalStateException("Cannot locate plugin source directory from ${new File('.').absolutePath}")
    }

    static final File SRC_ROOT = srcRoot()

    def "generic fs classes do not import resource-type-specific packages"() {
        expect:
        GENERIC_CLASSES.each { Class c ->
            final src = new File(SRC_ROOT, "${c.simpleName}.groovy").text
            assert !src.contains('io.seqera.tower.plugin.datalink.'), "${c.simpleName} must not import datalink package"
            assert !src.contains('io.seqera.tower.plugin.fs.handler.'), "${c.simpleName} must not import handler package"
            assert !src.contains('DataLink'), "${c.simpleName} must not reference data-link types"
            assert !src.contains('DatasetDto'), "${c.simpleName} must not reference DatasetDto"
            assert !src.contains('DatasetVersionDto'), "${c.simpleName} must not reference DatasetVersionDto"
        }
    }

    def "generic fs classes do not carry resource-type string literals"() {
        expect:
        GENERIC_CLASSES.each { Class c ->
            final src = new File(SRC_ROOT, "${c.simpleName}.groovy").text
            assert !src.contains("'datasets'"), "${c.simpleName} must not hard-code the 'datasets' resource type"
            assert !src.contains('"datasets"'), "${c.simpleName} must not hard-code the 'datasets' resource type"
            assert !src.contains("'data-links'"), "${c.simpleName} must not hard-code the 'data-links' resource type"
            assert !src.contains('"data-links"'), "${c.simpleName} must not hard-code the 'data-links' resource type"
        }
    }

    def "both handlers implement the ResourceTypeHandler interface"() {
        expect:
        ResourceTypeHandler.isAssignableFrom(io.seqera.tower.plugin.fs.handler.DatasetsResourceHandler)
        ResourceTypeHandler.isAssignableFrom(io.seqera.tower.plugin.fs.handler.DataLinksResourceHandler)
    }

    def "handlers do not reference each other's resource type"() {
        expect:
        final datasetSrc = new File(SRC_ROOT, 'handler/DatasetsResourceHandler.groovy').text
        final dataLinkSrc = new File(SRC_ROOT, 'handler/DataLinksResourceHandler.groovy').text
        !datasetSrc.contains('DataLink')
        !datasetSrc.contains('datalink')
        !dataLinkSrc.contains('DatasetDto')
        !dataLinkSrc.contains('DatasetVersionDto')
    }
}
