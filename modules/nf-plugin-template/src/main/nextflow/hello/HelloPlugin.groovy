/*
 * Copyright 2021, Seqera Labs
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

package nextflow.hello

import groovy.transform.CompileStatic
import nextflow.plugin.BasePlugin
import nextflow.plugin.Scoped
import org.pf4j.PluginWrapper

/**
 * Implements the Hello plugins entry point
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class HelloPlugin extends BasePlugin {

    HelloPlugin(PluginWrapper wrapper) {
        super(wrapper)
    }
}
