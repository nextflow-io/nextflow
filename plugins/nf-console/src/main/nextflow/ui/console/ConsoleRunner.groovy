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
