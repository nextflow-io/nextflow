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

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.BuildInfo
import nextflow.exception.AbortOperationException
import org.pf4j.Plugin
import org.pf4j.PluginWrapper
import org.slf4j.Logger
import org.slf4j.LoggerFactory
/**
 * Base class for NF plugins
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
abstract class BasePlugin extends Plugin {

    private static Logger log = LoggerFactory.getLogger(BasePlugin)

    BasePlugin(PluginWrapper wrapper) {
        super(wrapper)
    }

    @PackageScope boolean verMatches(String requires, String current=BuildInfo.version) {
        return getWrapper()
                .getPluginManager()
                .getVersionManager()
                .checkVersionConstraint(current, requires)
    }

    @Override
    void start() {
        final desc = getWrapper().getDescriptor()
        final name = "${desc.pluginId}@${desc.version}"
        if( desc.requires && !verMatches(desc.requires)) {
            throw new AbortOperationException("Failed requirement - Plugin $name requires Nextflow version $desc.requires (current $BuildInfo.version)")
        }
        log.debug "Plugin started $name"
    }

    @Override
    void stop() {
        log.debug "Plugin stopped ${wrapper.descriptor.pluginId}"
    }
}
