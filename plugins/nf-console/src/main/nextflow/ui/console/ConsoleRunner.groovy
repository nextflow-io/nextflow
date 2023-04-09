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
package nextflow.ui.console

import javax.swing.UIManager

import groovy.util.logging.Slf4j
import nextflow.cli.CliOptions
import nextflow.util.LoggerHelper
import org.codehaus.groovy.runtime.StackTraceUtils

/**
 * Implement the {@link ConsoleExtension} to launch the NF console app.
 *
 * See {@link nextflow.cli.CmdConsole#run()}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ConsoleRunner implements ConsoleExtension {

    /**
     * Nextflow REPL entry point
     *
     * @param args
     */
    @Override
    void run(String... args) {
        CliOptions opts = new CliOptions()
        opts.logFile = '.nextflow-console.log'
        new LoggerHelper(opts).setup()

        if (args.length == 2 && args[1] == '--help') {
            println 'usage: nextflow console [filename]'
            return
        }

        // full stack trace should not be logged to the output window - GROOVY-4663
        java.util.logging.Logger.getLogger(StackTraceUtils.STACK_LOG_NAME).useParentHandlers = false

        //when starting via main set the look and feel to system
        UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName())

        def console = new Nextflow(ConsoleRunner.getClassLoader())
        console.useScriptClassLoaderForScriptExecution = true
        console.run()
        if (args.length == 2)
            try {
                console.loadScriptFile(args[1] as File)
            }
            catch( IOException e ) {
                log.warn("Can't open script file: ${args[1]}" )
            }

    }
}
