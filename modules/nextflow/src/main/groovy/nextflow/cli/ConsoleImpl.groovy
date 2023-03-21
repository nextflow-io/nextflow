/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.cli

import groovy.transform.CompileStatic
import nextflow.plugin.Plugins
import nextflow.ui.console.ConsoleExtension

/**
 * CLI `console` sub-command
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ConsoleImpl {

    interface Options {
        String getScript()
    }

    @Delegate
    private Options options

    ConsoleImpl(Options options) {
        this.options = options
    }

    void run() {
        Plugins.setup()
        Plugins.start('nf-console')
        final console = Plugins.getExtension(ConsoleExtension)
        if( !console )
            throw new IllegalStateException("Failed to find Nextflow Console extension")
        console.run(script)
    }
}
