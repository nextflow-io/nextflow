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

package nextflow.cli

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.plugin.Plugins
import nextflow.ui.console.ConsoleExtension

/**
 * Launch the Nextflow Console plugin
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CmdConsole {

    interface Options {
        String getScript()
    }

    @Parameters(commandDescription = "Launch Nextflow interactive console")
    static class V1 extends CmdBase implements Options {

        @Parameter(description = 'Nextflow console arguments')
        List<String> args = []

        @Override
        String getScript() {
            args.size() > 0 ? args[0] : null
        }

        @Override
        String getName() { 'console' }

        @Override
        void run() {
            new CmdConsole(this).run()
        }
    }

    @Delegate
    private Options options

    CmdConsole(Options options) {
        this.options = options
    }

    void run() {
        Plugins.init()
        Plugins.start('nf-console')
        final console = Plugins.getExtension(ConsoleExtension)
        if( !console )
            throw new IllegalStateException("Failed to find Nextflow Console extension")
        console.run(script)
    }
}
