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
@Parameters(commandDescription = "Launch Nextflow interactive console")
class CmdConsole extends CmdBase {

    @Parameter(description = 'Nextflow console arguments')
    List<String> args

    String getName() { 'console' }

    void run() {
        Plugins.init()
        Plugins.start('nf-console')
        final console = Plugins.getExtension(ConsoleExtension)
        if( !console )
            throw new IllegalStateException("Failed to find Nextflow Console extension")
        // normalise the console args prepending the `console` command itself
        if( args == null )
            args = []
        args.add(0, 'console')
        // go !
        console.run(args as String[])
    }
}
