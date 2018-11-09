/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
import static nextflow.Const.APP_BUILDNUM
import static nextflow.Const.APP_NAME
import static nextflow.Const.APP_VER
import static nextflow.Const.SPLASH

import java.lang.reflect.Field

import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.JCommander
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.beust.jcommander.Parameters
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.exception.AbortRunException
import nextflow.exception.ConfigParseException
import nextflow.trace.GraphObserver
import nextflow.trace.ReportObserver
import nextflow.trace.TimelineObserver
import nextflow.trace.TraceFileObserver
import nextflow.util.LoggerHelper
import org.codehaus.groovy.control.CompilationFailedException
import org.eclipse.jgit.api.errors.GitAPIException
/**
 * Main application entry point. It parses the command line and
 * launch the pipeline execution.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class Launcher {

    /**
     * Create the application command line parser
     *
     * @return An instance of {@code CliBuilder}
     */

    private JCommander jcommander

    private CliOptions options

    private boolean fullVersion

    private CmdBase command

    private String cliString

    private List<CmdBase> allCommands

    private List<String> normalizedArgs

    private boolean daemonMode

    private String colsString

    /**
     * Create a launcher object and parse the command line parameters
     *
     * @param args The command line arguments provided by the user
     */
    Launcher() {
        init()
    }

    protected void init() {
        allCommands = (List<CmdBase>)[
                new CmdClean(),
                new CmdClone(),
                new CmdCloud(),
                new CmdFs(),
                new CmdHistory(),
                new CmdInfo(),
                new CmdList(),
                new CmdLog(),
                new CmdLs(),
                new CmdPull(),
                new CmdRun(),
                new CmdKubeRun(),
                new CmdDrop(),
                new CmdConfig(),
                new CmdNode(),
                new CmdView(),
                new CmdHelp(),
                new CmdSelfUpdate()
        ]

        options = new CliOptions()
        jcommander = new JCommander(options)
        allCommands.each { cmd ->
            cmd.launcher = this;
            jcommander.addCommand(cmd.name, cmd)
        }
        jcommander.setProgramName( APP_NAME )
    }

    /**
     * Create the Jcommander 'interpreter' and parse the command line arguments
     */
    @PackageScope
    Launcher parseMainArgs(String... args) {
        this.cliString = System.getenv('NXF_CLI')
        this.colsString = System.getenv('COLUMNS')

        def cols = getColumns()
        if( cols )
            jcommander.setColumnSize(cols)

        normalizedArgs = normalizeArgs(args)
        jcommander.parse( normalizedArgs as String[] )
        fullVersion = '-version' in normalizedArgs
        command = allCommands.find { it.name == jcommander.getParsedCommand()  }
        // whether is running a daemon
        daemonMode = command instanceof CmdNode
        // set the log file name
        checkLogFileName()

        return this
    }

    private void checkLogFileName() {
        if( !options.logFile ) {
            if( isDaemon() )
                options.logFile = '.node-nextflow.log'
            else if( command instanceof CmdRun || options.debug || options.trace )
                options.logFile = ".nextflow.log"
        }
    }

    private short getColumns() {
        if( !colsString ) {
            return 0
        }

        try {
            colsString.toShort()
        }
        catch( Exception e ) {
            log.debug "Oops .. not a valid \$COLUMNS value: $colsString"
            return 0
        }
    }

    CliOptions getOptions() { options }

    List<String> getNormalizedArgs() { normalizedArgs }

    String getCliString() { cliString }

    boolean isDaemon() { daemonMode }

    /**
     * normalize the command line arguments to handle some corner cases
     */
    @PackageScope
    List<String> normalizeArgs( String ... args ) {

        def normalized = []
        int i=0
        while( true ) {
            if( i==args.size() ) { break }

            def current = args[i++]
            normalized << current

            // when the first argument is a file, it's supposed to be a script to be executed
            if( i==1 && !allCommands.find { it.name == current } && new File(current).isFile()  ) {
                normalized.add(0,CmdRun.NAME)
            }

            else if( current == '-resume' ) {
                if( i<args.size() && !args[i].startsWith('-') && (args[i]=='last' || args[i] =~~ /[0-9a-f]{8}\-[0-9a-f]{4}\-[0-9a-f]{4}\-[0-9a-f]{4}\-[0-9a-f]{8}/) ) {
                    normalized << args[i++]
                }
                else {
                    normalized << 'last'
                }
            }
            else if( current == '-test' && (i==args.size() || args[i].startsWith('-'))) {
                normalized << '%all'
            }

            else if( current == '-with-trace' && (i==args.size() || args[i].startsWith('-'))) {
                normalized << TraceFileObserver.DEF_FILE_NAME
            }

            else if( current == '-with-report' && (i==args.size() || args[i].startsWith('-'))) {
                normalized << ReportObserver.DEF_FILE_NAME
            }

            else if( current == '-with-timeline' && (i==args.size() || args[i].startsWith('-'))) {
                normalized << TimelineObserver.DEF_FILE_NAME
            }

            else if( current == '-with-dag' && (i==args.size() || args[i].startsWith('-'))) {
                normalized << GraphObserver.DEF_FILE_NAME
            }

            else if( current == '-with-docker' && (i==args.size() || args[i].startsWith('-'))) {
                normalized << '-'
            }

            else if( current == '-with-singularity' && (i==args.size() || args[i].startsWith('-'))) {
                normalized << '-'
            }

            else if( current == '-with-weblog' && (i==args.size() || args[i].startsWith('-'))) {
                normalized << '-'
            }

            else if( (current == '-N' || current == '-with-notification') && (i==args.size() || args[i].startsWith('-'))) {
                normalized << 'true'
            }

            else if( (current == '-K' || current == '-with-k8s') && (i==args.size() || args[i].startsWith('-'))) {
                normalized << 'true'
            }

            else if( current == '-syslog' && (i==args.size() || args[i].startsWith('-') || allCommands.find { it.name == args[i] } )) {
                normalized << 'localhost'
            }

            else if( current == '-dump-channels' && (i==args.size() || args[i].startsWith('-'))) {
                normalized << '*'
            }

            else if( current ==~ /^\-\-[a-zA-Z\d].*/ && !current.contains('=') ) {
                current += '='
                current += ( i<args.size() && isValue(args[i]) ? args[i++] : 'true' )
                normalized[-1] = current
            }

            else if( current ==~ /^\-process\..+/ && !current.contains('=')) {
                current += '='
                current += ( i<args.size() && isValue(args[i]) ? args[i++] : 'true' )
                normalized[-1] = current
            }

            else if( current ==~ /^\-cluster\..+/ && !current.contains('=')) {
                current += '='
                current += ( i<args.size() && isValue(args[i]) ? args[i++] : 'true' )
                normalized[-1] = current
            }

            else if( current ==~ /^\-executor\..+/ && !current.contains('=')) {
                current += '='
                current += ( i<args.size() && isValue(args[i]) ? args[i++] : 'true' )
                normalized[-1] = current
            }

            else if( current == 'run' && i<args.size() && args[i] == '-' ) {
                i++
                normalized << '-stdin'
            }
        }

        return normalized
    }

    static private boolean isValue( String x ) {
        if( !x ) return false                   // an empty string -> not a value
        if( x.size() == 1 ) return true         // a single char is not an option -> value true
        !x.startsWith('-') || x.isNumber() || x.contains(' ')
    }

    CmdBase findCommand( String cmdName ) {
        allCommands.find { it.name == cmdName }
    }

    /**
     * Print the usage string for the given command - or -
     * the main program usage string if not command is specified
     *
     * @param command The command for which get help or {@code null}
     * @return The usage string
     */
    void usage(String command = null ) {

        if( command ) {
            def exists = allCommands.find { it.name == command } != null
            if( !exists ) {
                println "Asking help for unknown command: $command"
                return
            }

            jcommander.usage(command)
            return
        }

        println "Usage: nextflow [options] COMMAND [arg...]\n"
        printOptions(CliOptions)
        printCommands(allCommands)
    }

    @CompileDynamic
    protected void printOptions(Class clazz) {
        List params = []
        for( Field f : clazz.getDeclaredFields() ) {
            def p = f.getAnnotation(Parameter)
            if(!p)
                p = f.getAnnotation(DynamicParameter)

            if( p && !p.hidden() && p.description() && p.names() )
                params.add(p)

        }

        params.sort(true) { it -> it.names()[0] }

        println "Options:"
        for( def p : params ) {
            println "  ${p.names().join(', ')}"
            println "     ${p.description()}"
        }
    }

    protected void printCommands(List<CmdBase> commands) {
        println "\nCommands:"

        int len = 0
        def all = new TreeMap<String,String>()
        new ArrayList<CmdBase>(commands).each {
            def description = it.getClass().getAnnotation(Parameters)?.commandDescription()
            if( description ) {
                all[it.name] = description
                if( it.name.size()>len ) len = it.name.size()
            }
        }

        all.each { String name, String desc ->
            print '  '
            print name.padRight(len)
            print '   '
            println desc
        }
        println ''
    }

    Launcher command( String[] args ) {
        /*
         * CLI argument parsing
         */
        try {
            parseMainArgs(args)
            LoggerHelper.configureLogger(this)
        }
        catch( ParameterException e ) {
            // print command line parsing errors
            // note: use system.err.println since if an exception is raised
            //       parsing the cli params the logging is not configured
            System.err.println "${e.getMessage()} -- Check the available commands and options and syntax with 'help'"
            System.exit(1)

        }
        catch( Throwable e ) {
            e.printStackTrace(System.err)
            System.exit(1)
        }
        return this
    }

    protected void checkForHelp() {
        if( options.help || !command || command.help ) {
            if( command instanceof UsageAware ) {
                (command as UsageAware).usage()
                // reset command to null to skip default execution
                command = null
                return
            }

            // replace the current command with the `help` command
            def target = command?.name
            command = allCommands.find { it instanceof CmdHelp }
            if( target ) {
                (command as CmdHelp).args = [target]
            }
        }

    }

    /**
     * Launch the pipeline execution
     */
    int run() {

        /*
         * setup environment
         */
        setupEnvironment()

        /*
         * Real execution starts here
         */
        try {
            log.debug '$> ' + cliString

            // -- print out the version number, then exit
            if ( options.version ) {
                println getVersion(fullVersion)
                return 0
            }

            // -- print out the program help, then exit
            checkForHelp()

            // launch the command
            command?.run()

            log.trace "Exit\n" + dumpThreads()
            return 0
        }

        catch( AbortRunException e ) {
            return(1)
        }

        catch ( AbortOperationException e ) {
            def message = e.getMessage()
            if( message ) System.err.println(message)
            log.debug ("Operation aborted", e.cause ?: e)
            return(1)
        }

        catch ( GitAPIException e ) {
            System.err.println e.getMessage() ?: e.toString()
            log.debug ("Operation aborted", e.cause ?: e)
            return(1)
        }

        catch( ConfigParseException e )  {
            def message = e.message
            if( e.cause?.message ) {
                message += "\n\n${e.cause.message.toString().indent('  ')}"
            }
            log.error(message, e.cause ?: e)
            return(1)
        }

        catch( CompilationFailedException e ) {
            log.error e.message
            return(1)
        }

        catch( IOException e ) {
            log.error(e.message, e)
            return(1)
        }

        catch( Throwable fail ) {
            log.error("@unknown", fail)
            return(1)
        }

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
     * <li>HTTP_PROXY</li>
     * <li>HTTPS_PROXY</li>
     */
    private void setupEnvironment() {

        setProxy('HTTP',System.getenv())
        setProxy('HTTPS',System.getenv())

        setProxy('http',System.getenv())
        setProxy('https',System.getenv())

    }

    /**
     * Setup proxy system properties
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
        assert qualifier in ['http','https','HTTP','HTTPS']
        def str = null
        def var = "${qualifier}_" + (qualifier.isLowerCase() ? 'proxy' : 'PROXY')

        // -- setup HTTP proxy
        try {
            List<String> proxy = parseProxy(str = env.get(var.toString()))
            if( proxy ) {
                log.debug "Setting $qualifier proxy: $proxy"
                System.setProperty("${qualifier.toLowerCase()}.proxyHost", proxy[0])
                if( proxy[1] ) System.setProperty("${qualifier.toLowerCase()}.proxyPort", proxy[1])
            }
        }
        catch ( MalformedURLException e ) {
            log.warn "Not a valid $qualifier proxy: '$str' -- Check the value of variable `$var` in your environment"
        }

    }

    /**
     * Parse a proxy URL string retrieving the host and port components
     *
     * @param value A proxy string e.g. {@code hostname}, {@code hostname:port}, {@code http://hostname:port}
     * @return A list object containing at least the host name and optionally a second entry for the port.
     *      An empty list if the specified value is empty
     *
     * @throws MalformedURLException when the specified value is not a valid proxy url
     */
    @PackageScope
    static List parseProxy( String value ) {
        List<String> result = []
        int p

        if( !value ) return result

        if( value.contains('://') ) {
            def url = new URL(value)
            result.add(url.host)
            if( url.port > 0 )
                result.add(url.port as String)

        }
        else if( (p=value.indexOf(':')) != -1 ) {
            result.add( value.substring(0,p) )
            result.add( value.substring(p+1) )
        }
        else {
            result.add( value )
        }

        return result
    }

    /**
     * Hey .. Nextflow starts here!
     *
     * @param args The program options as specified by the user on the CLI
     */
    public static void main(String... args)  {

        final launcher = DripMain.LAUNCHER ?: new Launcher()
        final status = launcher .command(args) .run()
        if( status )
            System.exit(status)
    }


    /**
     * Print the application version number
     * @param full When {@code true} prints full version number including build timestamp
     * @return The version number string
     */
    static String getVersion(boolean full = false) {

        if ( full ) {
            SPLASH
        }
        else {
            "${APP_NAME} version ${APP_VER}.${APP_BUILDNUM}"
        }

    }


}
