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

package nextflow.pixi

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.plugin.BasePlugin
import org.pf4j.PluginWrapper

/**
 * Nextflow Pixi Package Manager Plugin
 * 
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
@Slf4j
@CompileStatic
class PixiPlugin extends BasePlugin {

    PixiPlugin(PluginWrapper wrapper) {
        super(wrapper)
    }

    @Override
    void start() {
        log.info "Starting Pixi package manager plugin"
        super.start()
    }

    @Override
    void stop() {
        log.info "Stopping Pixi package manager plugin"
        super.stop()
    }
}