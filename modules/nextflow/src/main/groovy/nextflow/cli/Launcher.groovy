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

import java.lang.reflect.Field

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.exception.AbortOperationException
import nextflow.exception.AbortRunException
import nextflow.exception.ConfigParseException
import nextflow.exception.ScriptCompilationException
import nextflow.exception.ScriptRuntimeException
import nextflow.secret.SecretsLoader
import nextflow.util.Escape
import nextflow.util.LoggerHelper
import nextflow.util.ProxyConfig
import org.eclipse.jgit.api.errors.GitAPIException
import picocli.CommandLine
import picocli.CommandLine.ArgGroup
import picocli.CommandLine.Command
import picocli.CommandLine.HelpCommand
import picocli.CommandLine.ParseResult

/**
 * Main application entry point. It parses the command line and
 * launches the pipeline execution.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Command(
    name = 'nextflow',
    description = 'Nextflow CLI',
    subcommands = [
        CmdClean.class,
        CmdClone.class,
        CmdConfig.class,
        CmdConsole.class,
        CmdDrop.class,
        CmdFs.class,
        HelpCommand.class,
        CmdInfo.class,
        CmdKubeRun.class,
        CmdList.class,
        CmdLog.class,
        CmdNode.class,
        CmdPlugin.class,
        CmdPull.class,
        CmdRun.class,
        CmdSelfUpdate.class,
        CmdView.class
    ]
)
class Launcher extends CmdBase {

    @ArgGroup
    private CliOptions options

    private String cliString

    private boolean daemonMode

    Launcher() {
        this.options = new CliOptions()
    }

    private int executionStrategy(ParseResult parseResult) {

        // make command line string
        final args = parseResult.originalArgs() as String[]

        this.cliString = makeCli(System.getenv('NXF_CLI'), args)

        // set whether Nextflow is running as a daemon
        this.daemonMode = spec.commandLine().getCommand() == CmdNode

        // setup proxy environment
        setupEnvironment()

        // setup logging
        checkLogFileName()

        LoggerHelper.configureLogger(this)

        // launch the command
        log.debug '$> ' + cliString

        int exitCode = new CommandLine.RunLast().execute(parseResult)

        if( log.isTraceEnabled() )
            log.trace "Exit\n" + dumpThreads()

        return exitCode
    }

    private String makeCli(String cli, String... args) {
        if( !cli )
            cli = 'nextflow'
        if( !args )
            return cli
        def cmd = ' ' + args[0]
        int p = cli.indexOf(cmd)
        if( p!=-1 )
            cli = cli.substring(0,p)
        if( cli.endsWith('nextflow') )
            cli = 'nextflow'
        cli += ' ' + Escape.cli(args)
        return cli
    }

    private void checkLogFileName() {
        if( !options.logFile ) {
            if( isDaemon() )
                options.logFile = System.getenv('NXF_LOG_FILE') ?: '.node-nextflow.log'
            else if( spec.commandLine().getCommand() == CmdRun || options.debug || options.trace )
                options.logFile = System.getenv('NXF_LOG_FILE') ?: ".nextflow.log"
        }
    }

    CliOptions getOptions() { options }

    boolean isDaemon() { daemonMode }

    @Override
    void run() {
        // -- print out the version number, then exit
        if ( options.version ) {
            println getVersion(false)
            return
        }

        if ( options.fullVersion ) {
            println getVersion(true)
            return
        }

        // -- print out the program help, then exit
        spec.commandLine().usage(System.err)
    }

    /**
     * Dump th stack trace of current running threads
     * @return
     */
    private String dumpThreads() {

        def buffer = new StringBuffer()
        Map<Thread, StackTraceElement[]> m = Thread.getAllStackTraces();
        for(Map.Entry<Thread,  StackTraceElement[]> e : m.entrySet()) {
            buffer.append('\n').append(e.getKey().toString()).append('\n')
            for (StackTraceElement s : e.getValue()) {
                buffer.append("  " + s).append('\n')
            }
        }

        return buffer.toString()
    }

    /**
     * set up environment and system properties. It checks the following
     * environment variables:
     * <li>http_proxy</li>
     * <li>https_proxy</li>
     * <li>ftp_proxy</li>
     * <li>HTTP_PROXY</li>
     * <li>HTTPS_PROXY</li>
     * <li>FTP_PROXY</li>
     * <li>NO_PROXY</li>
     */
    private void setupEnvironment() {

        setProxy('HTTP',System.getenv())
        setProxy('HTTPS',System.getenv())
        setProxy('FTP',System.getenv())

        setProxy('http',System.getenv())
        setProxy('https',System.getenv())
        setProxy('ftp',System.getenv())

        setNoProxy(System.getenv())
    }

    /**
     * Set no proxy property if defined in the launching env
     *
     * See for details
     * https://docs.oracle.com/javase/8/docs/technotes/guides/net/proxies.html
     *
     * @param env
     */
    @PackageScope
    static void setNoProxy(Map<String,String> env) {
        final noProxy = env.get('NO_PROXY') ?: env.get('no_proxy')
        if(noProxy) {
            System.setProperty('http.nonProxyHosts', noProxy.tokenize(',').join('|'))
        }
    }


    /**
     * Setup proxy system properties and optionally configure the network authenticator
     *
     * See:
     * http://docs.oracle.com/javase/6/docs/technotes/guides/net/proxies.html
     * https://github.com/nextflow-io/nextflow/issues/24
     *
     * @param qualifier Either {@code http/HTTP} or {@code https/HTTPS}.
     * @param env The environment variables system map
     */
    @PackageScope
    static void setProxy(String qualifier, Map<String,String> env ) {
        assert qualifier in ['http','https','ftp','HTTP','HTTPS','FTP']
        def str = null
        def var = "${qualifier}_" + (qualifier.isLowerCase() ? 'proxy' : 'PROXY')

        // -- setup HTTP proxy
        try {
            final proxy = ProxyConfig.parse(str = env.get(var.toString()))
            if( proxy ) {
                // set the expected protocol
                proxy.protocol = qualifier.toLowerCase()
                log.debug "Setting $qualifier proxy: $proxy"
                System.setProperty("${qualifier.toLowerCase()}.proxyHost", proxy.host)
                if( proxy.port )
                    System.setProperty("${qualifier.toLowerCase()}.proxyPort", proxy.port)
                if( proxy.authenticator() ) {
                    log.debug "Setting $qualifier proxy authenticator"
                    Authenticator.setDefault(proxy.authenticator())
                }
            }
        }
        catch ( MalformedURLException e ) {
            log.warn "Not a valid $qualifier proxy: '$str' -- Check the value of variable `$var` in your environment"
        }

    }

    /**
     * Hey .. Nextflow starts here!
     *
     * @param args The program options as specified by the user on the CLI
     */
    static void main(String... args)  {

        try {
            // create launcher
            def launcher = new Launcher()
            def cmd = new CommandLine(launcher)
                .setExecutionStrategy(launcher::executionStrategy)

            // add secrets command if enabled
            if( SecretsLoader.isEnabled() )
                cmd.addSubcommand(new CmdSecret())

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
            System.err.println e.getMessage() ?: e.toString()
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
     * Print the application version number
     * @param full When {@code true} prints full version number including build timestamp
     * @return The version number string
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
