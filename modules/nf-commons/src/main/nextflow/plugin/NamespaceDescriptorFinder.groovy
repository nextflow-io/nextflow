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

package nextflow.plugin

import java.nio.file.Path

import groovy.transform.CompileStatic
import org.pf4j.DefaultPluginDescriptor
import org.pf4j.ManifestPluginDescriptorFinder
import org.pf4j.PluginDescriptor
import org.pf4j.PluginRuntimeException

/**
 * Reads the plugin ID as `<namespace>/<Plugin-Id>`, taking the namespace from
 * `Plugin-Namespace` or, for a release built without it, from the store folder name
 */
@CompileStatic
class NamespaceDescriptorFinder extends ManifestPluginDescriptorFinder {

    static final String PLUGIN_NAMESPACE = 'Plugin-Namespace'

    @Override
    PluginDescriptor find(Path pluginPath) {
        return find(pluginPath, pluginPath.fileName.toString())
    }

    /**
     * @param folderName The store folder name, which can differ from the path name
     */
    PluginDescriptor find(Path pluginPath, String folderName) {
        final manifest = readManifest(pluginPath)
        final result = (Descriptor) createPluginDescriptor(manifest)
        final id = result.pluginId
        if( !id || id.contains('/') )
            return result
        final declared = manifest.mainAttributes.getValue(PLUGIN_NAMESPACE)
        if( declared != null && !isValid(declared, id) )
            throw new PluginRuntimeException("Invalid {} '{}' in plugin {}", PLUGIN_NAMESPACE, declared, pluginPath)
        final ns = declared ?: folderNamespace(folderName)
        if( ns )
            result.setPluginId("$ns/$id")
        return result
    }

    private static boolean isValid(String ns, String id) {
        try {
            PluginRef.parse("$ns/$id")
            return true
        }
        catch( IllegalArgumentException e ) {
            return false
        }
    }

    @Override
    protected DefaultPluginDescriptor createPluginDescriptorInstance() {
        return new Descriptor()
    }

    private static String folderNamespace(String folderName) {
        final id = PluginRef.fromDirName(folderName)
        return id.contains('/') ? id.substring(0, id.indexOf('/')) : null
    }

    // pf4j keeps `setPluginId` protected
    static class Descriptor extends DefaultPluginDescriptor {
        @Override
        protected DefaultPluginDescriptor setPluginId(String id) { super.setPluginId(id) }
    }
}
