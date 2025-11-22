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
 */
package nextflow.ui.console

import java.util.jar.JarFile
import java.util.zip.ZipException

import groovy.transform.CompileStatic
import nextflow.plugin.BasePlugin
import org.codehaus.groovy.reflection.CachedClass
import org.codehaus.groovy.reflection.ClassInfo
import org.codehaus.groovy.runtime.m12n.ExtensionModuleScanner
import org.codehaus.groovy.runtime.metaclass.MetaClassRegistryImpl
import org.pf4j.PluginClassLoader
import org.pf4j.PluginWrapper
import org.slf4j.Logger
import org.slf4j.LoggerFactory
/**
 * Nextflow plugin for Console extension
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ConsolePlugin extends BasePlugin {

    private static Logger log = LoggerFactory.getLogger(ConsolePlugin)

    ConsolePlugin(PluginWrapper wrapper) {
        super(wrapper)
    }

    @Override
    void start() {
        super.start()
        // The console UI requires some Groovy extensions method that are not loaded automatically
        // using the plugin system classloader.
        // see
        // - org.apache.groovy.swing.extensions.SwingExtensions
        // - META-INF/groovy/org.apache.groovy.runtime.ExtensionModule in the 'groovy-swing' JAR
        loadExtensions()
    }

    protected void loadExtensions() {
        if( wrapper.getPluginClassLoader() !instanceof PluginClassLoader )
            return

        final pcl = (PluginClassLoader) wrapper.getPluginClassLoader()
        for( URL it  : pcl.getURLs() ) {
            log.trace "Checking console lib for groovy module extensions: $it"
            processCategoryMethods( wrapper.pluginClassLoader, new File(it.getFile()) )
        }
    }

    /**
     * Install extension method dynamically. Taken from {@link groovy.grape.GrapeIvy#processCategoryMethods(java.lang.ClassLoader, java.io.File)}
     *
     * @param loader Plugin class loader
     * @param file Extension library jar file
     */
    private void processCategoryMethods(ClassLoader loader, File file) {
        // register extension methods if jar
        if (file.name.toLowerCase().endsWith('.jar')) {
            def mcRegistry = GroovySystem.metaClassRegistry
            if (mcRegistry instanceof MetaClassRegistryImpl) {
                try (JarFile jar = new JarFile(file)) {
                    def entry = jar.getEntry(ExtensionModuleScanner.MODULE_META_INF_FILE)
                    if (!entry) {
                        entry = jar.getEntry(ExtensionModuleScanner.LEGACY_MODULE_META_INF_FILE)
                    }
                    if (entry) {
                        Properties props = new Properties()

                        try (InputStream is = jar.getInputStream(entry)) {
                            props.load(is)
                        }

                        Map<CachedClass, List<MetaMethod>> metaMethods = new HashMap<CachedClass, List<MetaMethod>>()
                        mcRegistry.registerExtensionModuleFromProperties(props, loader, metaMethods)
                        // add old methods to the map
                        metaMethods.each { CachedClass c, List<MetaMethod> methods ->
                            // GROOVY-5543: if a module was loaded using grab, there are chances that subclasses
                            // have their own ClassInfo, and we must change them as well!
                            Set<CachedClass> classesToBeUpdated = [c].toSet()
                            ClassInfo.onAllClassInfo { ClassInfo info ->
                                if (c.theClass.isAssignableFrom(info.cachedClass.theClass)) {
                                    classesToBeUpdated << info.cachedClass
                                }
                            }
                            classesToBeUpdated*.addNewMopMethods(methods)
                        }
                    }
                } catch (ZipException zipException) {
                    throw new RuntimeException("Grape could not load jar '$file'", zipException)
                }
            }
        }
    }

}
