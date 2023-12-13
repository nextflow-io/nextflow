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

package nextflow.cli.v2

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.cli.CliOptions
import nextflow.exception.AbortOperationException
import nextflow.exception.AbortRunException
import nextflow.exception.ConfigParseException
import nextflow.exception.ScriptCompilationException
import nextflow.exception.ScriptRuntimeException
import nextflow.util.Escape
import nextflow.util.LoggerHelper
import nextflow.util.ProxyHelper
import org.eclipse.jgit.api.errors.GitAPIException
import picocli.CommandLine
import picocli.CommandLine.ArgGroup
import picocli.CommandLine.Command
import picocli.CommandLine.HelpCommand
import picocli.CommandLine.IExecutionExceptionHandler
import picocli.CommandLine.Option
import picocli.CommandLine.ParseResult

/**
 * Main application entry point. It parses the command line and
 * launches the pipeline execution.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
@Command(
    name = 'nf',
    description = 'Nextflow CLI v2',
    subcommands = [
        CleanCmd.class,
        CloneCmd.class,
        ConfigCmd.class,
        ConsoleCmd.class,
        DropCmd.class,
        FsCmd.class,
        HelpCommand.class,
        InfoCmd.class,
        InspectCmd.class,
        ListCmd.class,
        LogCmd.class,
        NodeCmd.class,
        PluginCmd.class,
        PullCmd.class,
        RunCmd.class,
        SecretsCmd.class,
        SelfUpdateCmd.class,
        ViewCmd.class
    ]
)
class Launcher extends AbstractCmd {

    @ArgGroup(validate = false)
    private CliOptionsV2 options

    private String cliString

    private boolean daemonMode

    Launcher() {
        this.options = new CliOptionsV2()
    }

    protected int executionStrategy(ParseResult parseResult) {
        if( parseResult.subcommand() )
            parseResult = parseResult.subcommand()

        def command = parseResult.commandSpec().commandLine().getCommand()
        def args = parseResult.originalArgs() as String[]

        // make command line string
        this.cliString = makeCli(System.getenv('NXF_CLI'), args)

        // whether is running a daemon
        this.daemonMode = command instanceof NodeCmd

        // set the log file name
        if( !options.logFile ) {
            if( isDaemon() )
                options.logFile = System.getenv('NXF_LOG_FILE') ?: '.node-nextflow.log'
            else if( command instanceof RunCmd || options.debug || options.trace )
                options.logFile = System.getenv('NXF_LOG_FILE') ?: '.nextflow.log'
        }

        LoggerHelper.configureLogger(options, isDaemon())

        // setup proxy environment
        ProxyHelper.setupEnvironment()

        // launch the command
        log.debug '$> ' + cliString

        int exitCode = new CommandLine.RunLast().execute(parseResult)

        if( log.isTraceEnabled() )
            log.trace "Exit\n" + dumpThreads()

        return exitCode
    }

    protected String makeCli(String cli, String... args) {
        if( !cli )
            cli = 'nf'
        if( !args )
            return cli
        def cmd = ' ' + args[0]
        int p = cli.indexOf(cmd)
        if( p!=-1 )
            cli = cli.substring(0,p)
        if( cli.endsWith('nf') )
            cli = 'nf'
        cli += ' ' + Escape.cli(args)
        return cli
    }

    CliOptionsV2 getOptions() { options }

    String getCliString() { cliString }

    boolean isDaemon() { daemonMode }

    @Override
    void run() {
        // -- print out the version number, then exit
        if( options.version ) {
            println getVersion(false)
            return
        }

        if( options.fullVersion ) {
            println getVersion(true)
            return
        }

        // -- print out the program help, then exit
        spec.commandLine().usage(System.err)
    }

    /**
     * Dump the stack trace of current running threads
     */
    private String dumpThreads() {

        def buffer = new StringBuffer()
        Map<Thread, StackTraceElement[]> m = Thread.getAllStackTraces()
        for(Map.Entry<Thread,  StackTraceElement[]> e : m.entrySet()) {
            buffer.append('\n').append(e.getKey().toString()).append('\n')
            for (StackTraceElement s : e.getValue()) {
                buffer.append("  " + s).append('\n')
            }
        }

        return buffer.toString()
    }

    static class ExecutionExceptionHandler implements IExecutionExceptionHandler {
        int handleExecutionException(Exception ex, CommandLine cmd, ParseResult parseResult) {
            // bold red error message
            cmd.getErr().println(cmd.getColorScheme().errorText(ex.getMessage() ?: ''))

            return cmd.getExitCodeExceptionMapper() != null
                ? cmd.getExitCodeExceptionMapper().getExitCode(ex)
                : cmd.getCommandSpec().exitCodeOnExecutionException()
        }
    }

    /**
     * Hey .. Nextflow starts here!
     *
     * @param args The program options as specified by the user on the CLI
     */
    static void main(String[] args) {
        try {
            // create launcher
            def launcher = new Launcher()
            def cmd = new CommandLine(launcher)
                .setExecutionStrategy(launcher::executionStrategy)
                .setExecutionExceptionHandler(new ExecutionExceptionHandler())
                .setAllowSubcommandsAsOptionParameters(true)
                .setUnmatchedOptionsArePositionalParams(true)

            // when the first argument is a file, it's supposed to be a script to be executed
            if( args.length > 0 && !cmd.getCommandSpec().subcommands().containsKey(args[0]) && new File(args[0]).isFile() ) {
                def argsList = args as List<String>
                argsList.add(0, 'run')
                args = argsList as String[]
            }

            // execute command
            System.exit(cmd.execute(args))
        }

        catch( AbortRunException e ) {
            System.exit(1)
        }

        catch( AbortOperationException e ) {
            def message = e.getMessage()
            if( message ) System.err.println(message)
            log.debug ("Operation aborted", e.cause ?: e)
            System.exit(1)
        }

        catch( GitAPIException e ) {
            System.err.println (e.getMessage() ?: e.toString())
            log.debug ("Operation aborted", e.cause ?: e)
            System.exit(1)
        }

        catch( ConfigParseException e )  {
            def message = e.message
            if( e.cause?.message ) {
                message += "\n\n${e.cause.message.toString().indent('  ')}"
            }
            log.error(message, e.cause ?: e)
            System.exit(1)
        }

        catch( ScriptCompilationException e ) {
            log.error(e.message, e)
            System.exit(1)
        }

        catch( ScriptRuntimeException | IllegalArgumentException e) {
            log.error(e.message, e)
            System.exit(1)
        }

        catch( IOException e ) {
            log.error(e.message, e)
            System.exit(1)
        }

        catch( Throwable fail ) {
            log.error("@unknown", fail)
            System.exit(1)
        }
    }


    /**
     * Print the version number.
     *
     * @param full When {@code true} prints full version number including build timestamp
     */
    static String getVersion(boolean full = false) {

        if ( full ) {
            Const.SPLASH
        }
        else {
            "${Const.APP_NAME} version ${Const.APP_VER}.${Const.APP_BUILDNUM}"
        }

    }

}
