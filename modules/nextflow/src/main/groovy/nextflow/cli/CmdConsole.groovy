/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

/**
 * Launch the Nextflow Console plugin
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Command(name = 'console', description = "Launch Nextflow interactive console")
class CmdConsole extends CmdBase {

    @Parameters(arity = '0..1', description = 'script filename')
    String script

    @Override
    void run() {
        Plugins.setup()
        Plugins.start('nf-console')
        final console = Plugins.getExtension(ConsoleExtension)
        if( !console )
            throw new IllegalStateException("Failed to find Nextflow Console extension")
        console.run(script)
    }
}
