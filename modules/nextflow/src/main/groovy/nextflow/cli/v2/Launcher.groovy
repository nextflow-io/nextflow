/*
 * Copyright 2023, Seqera Labs
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
import nextflow.cli.ILauncherOptions
import nextflow.exception.AbortOperationException
import nextflow.exception.AbortRunException
import nextflow.exception.ConfigParseException
import nextflow.exception.ScriptCompilationException
import nextflow.exception.ScriptRuntimeException
import nextflow.secret.SecretsLoader
import nextflow.util.LoggerHelper
import nextflow.util.ProxyHelper
import org.eclipse.jgit.api.errors.GitAPIException
import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.HelpCommand
import picocli.CommandLine.Option
import picocli.CommandLine.ParseResult
import picocli.CommandLine.ScopeType

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
        ListCmd.class,
        LogCmd.class,
        NodeCmd.class,
        PluginCmd.class,
        PullCmd.class,
        RunCmd.class,
        SelfUpdateCmd.class,
        ViewCmd.class
    ]
)
class Launcher extends AbstractCmd implements ILauncherOptions {

    boolean ansiLogCli

    void setAnsiLog(boolean value) { ansiLogCli = value }

    @Option(names = ['--bg'], arity = '0', description = 'Execute nextflow in background')
    boolean background

    @Option(names = ['-C'], description = 'Use the specified configuration file(s), ignoring any defaults')
    List<String> config

    @Option(names = ['-c','--config'], description = 'Add the specified file to configuration set')
    List<String> userConfig

    @Option(names = ['--config-ignore-includes'], description = 'Disable the parsing of config includes')
    boolean ignoreConfigIncludes

    @Option(names = ['-D'], description = 'Set JVM properties')
    Map<String,String> jvmOpts = [:]

    @Option(names = ['--debug'], description = 'Enable DEBUG level logging for the specified package name -- multiple packages can be provided as a comma-separated list (e.g. \'-debug nextflow,io.seqera\')', hidden = true)
    List<String> debug

    @Option(names = ['-d','--dockerize'], arity = '0', description = 'Launch Nextflow via Docker (experimental)')
    boolean dockerize

    @Option(names = ['--log'], description = 'Set the log file path')
    String logFile

    @Option(names = ['-q','--quiet'], description = 'Do not print information messages')
    boolean quiet

    @Option(names = ['--self-update'], arity = '0', description = 'Update Nextflow to the latest version', hidden = true)
    boolean selfUpdate

    @Option(names = ['--syslog'], description = 'Send logs to syslog server (e.g. localhost:514)')
    String syslog

    @Option(names = ['--trace'], description = 'Enable TRACE level logging for the specified package name -- multiple packages can be provided as a comma-separated list (e.g. \'-trace nextflow,io.seqera\')')
    List<String> trace

    @Option(names = ['-v'], description = 'Print the version number and exit')
    boolean version

    @Option(names = ['-V','--version'], description = 'Print the full version info and exit')
    boolean fullVersion

    boolean isDaemon() {
        spec.commandLine().getCommand() == NodeCmd
    }

    private int executionStrategy(ParseResult parseResult) {
        init()
        new CommandLine.RunLast().execute(parseResult)
    }

    private void init() {
        // setup logging
        if( !logFile ) {
            if( isDaemon() )
                logFile = System.getenv('NXF_LOG_FILE') ?: '.node-nextflow.log'
            else if( spec.commandLine().getCommand() == RunCmd || debug || trace )
                logFile = System.getenv('NXF_LOG_FILE') ?: ".nextflow.log"
        }

        LoggerHelper.configureLogger(this, isDaemon())

        // setup proxy environment
        ProxyHelper.setupEnvironment()

        // log command string
        log.debug '$>' + cliString
    }

    @Override
    Integer call() {
        // print version if specified
        if ( version ) {
            println "${Const.APP_NAME} version ${Const.APP_VER}.${Const.APP_BUILDNUM}"
            return 0
        }

        // print full version if specified
        if ( fullVersion ) {
            println Const.SPLASH
            return 0
        }

        // print usage and exit
        spec.commandLine().usage(System.err)
        return -1
    }

    /**
     * Hey .. Nextflow starts here!
     *
     * @param args The program options as specified by the user on the CLI
     */
    static void main(String[] args) {
        try {
            // create launcher
            Launcher launcher = new Launcher()
            CommandLine commandLine = new CommandLine(launcher)
                .setExecutionStrategy(launcher::executionStrategy)

            // add secrets command if enabled
            if( SecretsLoader.isEnabled() )
                commandLine.addSubcommand('secrets', new SecretsCmd())

            // execute command
            System.exit(commandLine.execute(args))
        }

        catch( AbortOperationException e ) {
            System.err.println (e.message ?: "Unknown abort reason")
            System.exit(1)
        }

        catch( Throwable e ) {
            e.printStackTrace(System.err)
            System.exit(1)
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

}
