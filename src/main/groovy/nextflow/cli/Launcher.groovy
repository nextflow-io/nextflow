/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.cli
import static Const.APP_BUILDNUM
import static Const.APP_NAME
import static Const.APP_VER
import static Const.SEE_LOG_FOR_DETAILS
import static Const.SPLASH

import com.beust.jcommander.JCommander
import com.beust.jcommander.ParameterException
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.ExitCode
import nextflow.daemon.DaemonLauncher
import nextflow.exception.AbortOperationException
import nextflow.exception.ConfigParseException
import nextflow.executor.ServiceName
import nextflow.script.ConfigBuilder
import nextflow.util.LoggerHelper
import nextflow.util.ServiceDiscover
import org.eclipse.jgit.api.errors.GitAPIException
/**
 * Main application entry point. It parses the command line and
 * launch the pipeline execution.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class Launcher implements ExitCode {

    /**
     * Create the application command line parser
     *
     * @return An instance of {@code CliBuilder}
     */

    private JCommander jcommander

    private CliOptions options

    private boolean fullVersion

    private CmdBase command

    private String[] args

    private String cliString

    private List<CmdBase> allCommands

    private List<String> normalizedArgs

    /**
     * Create a launcher object and parse the command line parameters
     *
     * @param args The command line arguments provided by the user
     */
    protected Launcher(String... args) {
        this.args = args
        this.cliString = System.getenv('NXF_CLI')

        try {
            // -- parse the program arguments - and - configure the logger
            parseMainArgs()
            LoggerHelper.configureLogger(options)
        }
        catch( ParameterException e ) {
            // print command line parsing errors
            // note: use system.err.println since if an exception is raised
            //       parsing the cli params the logging is not configured
            System.err.println "${e.getMessage()} -- Check the available commands and options and syntax with 'help'"
            System.exit( INVALID_COMMAND_LINE_PARAMETER )

        }
        catch( Throwable e ) {
            e.printStackTrace(System.err)
            System.exit( UNKNOWN_ERROR )
        }
    }


    /** ONLY FOR TEST */
    protected Launcher () { }


    /**
     * Create the Jcommander 'interpreter' and parse the command line arguments
     */
    @PackageScope
    Launcher parseMainArgs() {

        allCommands = (List<CmdBase>)[
                new CmdClone(),
                new CmdHistory(),
                new CmdInfo(),
                new CmdList(),
                new CmdPull(),
                new CmdRun(),
                new CmdDrop(),
                new CmdConfig(),
                new CmdNode(),
                new CmdHelp()
        ]

        def cols = getColumns()
        normalizedArgs = normalizeArgs(args)
        options = new CliOptions()
        jcommander = new JCommander(options)
        allCommands.each { cmd ->
            cmd.launcher = this;
            jcommander.addCommand(cmd.name, cmd)
        }
        jcommander.setProgramName( APP_NAME )
        if( cols )
            jcommander.setColumnSize(cols)
        jcommander.parse( normalizedArgs as String[] )
        fullVersion = '-version' in normalizedArgs
        command = allCommands.find { it.name == jcommander.getParsedCommand()  }

        return this
    }

    private int getColumns() {
        try {
            System.getenv('COLUMNS').toInteger()
        }
        catch( Exception e ) {
            return 0
        }
    }

    JCommander getJcommander() { jcommander }

    CliOptions getOptions() { options }

    List<String> getNormalizedArgs() { normalizedArgs }

    String getCliString() { cliString }

    /**
     * normalize the command line arguments to handle some corner cases
     */
    @PackageScope
    static List<String> normalizeArgs( String ... args ) {

        def normalized = []
        int i=0
        while( true ) {
            if( i==args.size() ) { break }

            def current = args[i++]
            normalized << current

            if( current == '-resume' ) {
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

            else if( current == '-with-drmaa' && (i==args.size() || args[i].startsWith('-'))) {
                normalized << '-'
            }

            else if( current == '-with-trace' && (i==args.size() || args[i].startsWith('-'))) {
                normalized << 'trace.csv'
            }

            else if( current ==~ /^\-\-[a-zA-Z\d].*/ && !current.contains('=')) {
                current += '='
                current += ( i<args.size() ? args[i++] : 'true' )
                normalized[-1] = current
            }

            else if( current ==~ /^\-process\..+/ && !current.contains('=')) {
                current += '='
                current += ( i<args.size() ? args[i++] : 'true' )
                normalized[-1] = current
            }

            else if( current ==~ /^\-daemon\..+/ && !current.contains('=')) {
                current += '='
                current += ( i<args.size() ? args[i++] : 'true' )
                normalized[-1] = current
            }

            else if( current ==~ /^\-executor\..+/ && !current.contains('=')) {
                current += '='
                current += ( i<args.size() ? args[i++] : 'true' )
                normalized[-1] = current
            }

            else if( current == 'run' && i<args.size() && args[i] == '-' ) {
                i++
                normalized << '-stdin'
            }
        }

        return normalized
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

        def list = new ArrayList<CmdBase>(allCommands).findAll { it.name != CmdNode.NAME }
        println "Usage: nextflow [options] COMMAND [arg...]\n"
        println "Commands: "

        int len = 0
        def all = new TreeMap<String,String>()
        list.each {
            def description = it.getClass().getAnnotation(Parameters)?.commandDescription()
            all[it.name] = description ?: '-'
            if( it.name.size()>len ) len = it.name.size()
        }

        all.each { String name, String desc ->
            print '  '
            print name.padRight(len)
            print '   '
            println desc
        }
        println ''
    }

    /**
     * Launch the pipeline execution
     */
    protected void run() {

        try {
            log.debug '$> ' + cliString

            // -- print out the version number, then exit
            if ( options.version ) {
                println getVersion(fullVersion)
                System.exit(OK)
            }

            // -- print out the program help, then exit
            if( options.help || !command || command.help ) {
                def target = command?.name
                command = allCommands.find { it instanceof CmdHelp }
                if( target )
                    (command as CmdHelp).args = [target]
            }

            // launch the command
            command.run()

        }

        catch ( GitAPIException | AbortOperationException e ) {
            System.err.println e.getMessage() ?: e.toString()
            log.debug ("Operation aborted", e.cause ?: e)
            System.exit(COMMAND_RUNTIME_ERROR)
        }

        catch( ConfigParseException e )  {
            log.error("${e.message}\n\n${indent(e.cause?.message?.toString(), '  ')}\n  ${SEE_LOG_FOR_DETAILS}\n", e.cause ?: e)
            System.exit(INVALID_CONFIG)
        }

        catch( Throwable fail ) {
            log.error("${fail.toString()} ${SEE_LOG_FOR_DETAILS}", fail)
            System.exit(UNKNOWN_ERROR)
        }

    }

    /**
     * Hey .. Nextflow starts here!
     *
     * @param args The program options as specified by the user on the CLI
     */
    public static void main(String... args)  {

        new Launcher(args).run()

    }


    public static String indent( String str, String separator = ' ') {
        def result = new StringBuilder()
        str?.eachLine { result << separator << it << '\n' }
        return result.toString()
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
