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

package nextflow.plugin

import java.nio.file.Files
import java.nio.file.Path
import java.util.jar.Manifest

import org.pf4j.ManifestPluginDescriptorFinder
/**
 * Plugin finder specialised for test environment. It looks for
 * the plugin 'MANIFEST.MF' file in the 'testFixtures' resources directory
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TestPluginDescriptorFinder extends ManifestPluginDescriptorFinder {

    @Override
    protected Manifest readManifestFromDirectory(Path pluginPath) {
        if( !Files.isDirectory(pluginPath) )
            return null

        final manifestPath = pluginPath.resolve('build/resources/testFixtures/META-INF/MANIFEST.MF')
        if( !Files.exists(manifestPath) )
            return null

        final input = Files.newInputStream(manifestPath)
        return new Manifest(input)
    }
}
