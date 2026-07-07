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

package nextflow.install2r

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.plugin.BasePlugin
import org.pf4j.PluginWrapper

/**
 * Nextflow install2.r Package Manager Plugin
 */
@Slf4j
@CompileStatic
class Install2rPlugin extends BasePlugin {

    Install2rPlugin(PluginWrapper wrapper) {
        super(wrapper)
    }

    @Override
    void start() {
        log.debug "Starting R install2.r package manager plugin"
        super.start()
    }

    @Override
    void stop() {
        log.debug "Stopping R install2.r package manager plugin"
        super.stop()
    }
}
