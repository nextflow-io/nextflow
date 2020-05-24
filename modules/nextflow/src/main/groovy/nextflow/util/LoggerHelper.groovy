/*
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

package nextflow.util

import static nextflow.Const.*

import java.lang.reflect.Field
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.util.concurrent.atomic.AtomicBoolean
import java.util.regex.Pattern

import ch.qos.logback.classic.Level
import ch.qos.logback.classic.LoggerContext
import ch.qos.logback.classic.encoder.PatternLayoutEncoder
import ch.qos.logback.classic.net.SyslogAppender
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.classic.spi.ThrowableProxy
import ch.qos.logback.core.Appender
import ch.qos.logback.core.AppenderBase
import ch.qos.logback.core.ConsoleAppender
import ch.qos.logback.core.CoreConstants
import ch.qos.logback.core.FileAppender
import ch.qos.logback.core.LayoutBase
import ch.qos.logback.core.encoder.LayoutWrappingEncoder
import ch.qos.logback.core.filter.Filter
import ch.qos.logback.core.joran.spi.NoAutoStart
import ch.qos.logback.core.rolling.FixedWindowRollingPolicy
import ch.qos.logback.core.rolling.RollingFileAppender
import ch.qos.logback.core.rolling.TriggeringPolicyBase
import ch.qos.logback.core.spi.FilterReply
import ch.qos.logback.core.util.FileSize
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Global
import nextflow.Session
import nextflow.cli.CliOptions
import nextflow.cli.Launcher
import nextflow.exception.AbortOperationException
import nextflow.exception.ProcessException
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.OpCall
import nextflow.file.FileHelper
import nextflow.script.BaseScript
import nextflow.script.ChainableDef
import nextflow.script.ComponentDef
import nextflow.script.CompositeDef
import nextflow.script.FunctionDef
import nextflow.script.ScriptMeta
import nextflow.script.WorkflowBinding
import nextflow.script.WorkflowDef
import org.apache.commons.lang.exception.ExceptionUtils
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import org.slf4j.Marker
import org.slf4j.MarkerFactory
/**
 * Helper class to setup the logging subsystem
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class LoggerHelper {

    static private Logger log = LoggerFactory.getLogger(LoggerHelper)

    static public Marker STICKY = MarkerFactory.getMarker('sticky')

    static private String STARTUP_ERROR = 'startup failed:\n'

    static private String logFileName

    private CliOptions opts

    private boolean rolling = false

    private String syslog

    private boolean daemon = false

    private int minIndex = 1

    private int maxIndex = 9

    private Map<String,Level> packages = [:]

    private LoggerContext loggerContext

    private Appender consoleAppender

    private Appender syslogAppender

    private Appender fileAppender

    LoggerHelper setRolling( boolean value ) {
        this.rolling = value
        return this
    }

    LoggerHelper setSyslog( String value ) {
        if(value) {
            this.syslog = value
        }
        return this
    }

    LoggerHelper setDaemon( boolean value ) {
        this.daemon = value
        return this
    }

    LoggerHelper setMinIndex( int value ) {
        this.minIndex = value
        return this
    }

    LoggerHelper setMaxIndex( int value ) {
        this.maxIndex = value
        return this
    }

    LoggerHelper(CliOptions opts) {
        this.opts = opts
        this.loggerContext = (LoggerContext) LoggerFactory.getILoggerFactory()
    }

    private boolean containsClassName(List<String> debug, List<String> trace, String prefix) {
        for( String it : debug ) {
            if( it.startsWith(prefix) ) return true
        }
        for( String it : trace ) {
            if( it.startsWith(prefix) ) return true
        }
        return false
    }

    void setup() {
        logFileName = opts.logFile ?: System.getenv('NXF_LOG_FILE')

        final boolean quiet = opts.quiet
        final List<String> debugConf = opts.debug ?: new ArrayList<String>()
        final List<String> traceConf = opts.trace ?: ( System.getenv('NXF_TRACE')?.tokenize(', ') ?: new ArrayList<String>())

        // Reset all the logger
        final root = loggerContext.getLogger('ROOT')
        root.detachAndStopAllAppenders()

        // -- define the console appender
        packages[MAIN_PACKAGE] = quiet ? Level.WARN : Level.INFO

        // -- add the S3 uploader by default
        if( !containsClassName(debugConf,traceConf, 'com.upplication.s3fs') )
            debugConf << S3_UPLOADER_CLASS
        if( !containsClassName(debugConf,traceConf, 'io.seqera') )
            debugConf << 'io.seqera'

        for( def it : debugConf ) {
            packages[it] = Level.DEBUG
        }
        for( def it : traceConf ) {
            packages[it] = Level.TRACE
        }

        // -- the console appender
        this.consoleAppender = createConsoleAppender()

        // -- the syslog appender
        this.syslogAppender = createSyslogAppender()

        // -- the file appender
        this.fileAppender = rolling ? createRollingAppender() : createFileAppender()

        // -- configure the ROOT logger
        root.setLevel(Level.INFO)
        if( fileAppender )
            root.addAppender(fileAppender)
        if( consoleAppender )
            root.addAppender(consoleAppender)
        if( syslogAppender )
            root.addAppender(syslogAppender)

        // -- main package logger
        def mainLevel = packages[MAIN_PACKAGE] == Level.TRACE ? Level.TRACE : Level.DEBUG
        def logger = createLogger(MAIN_PACKAGE, mainLevel)

        // -- set AWS lib level to WARN to reduce noise in the log file
        final AWS = 'com.amazonaws'
        if( !debugConf.contains(AWS) && !traceConf.contains(AWS)) {
            createLogger(AWS, Level.WARN)
        }

        // -- debug packages specified by the user
        for( String clazz : debugConf ) {
            createLogger(clazz, Level.DEBUG)
        }
        // -- trace packages specified by the user
        for( String clazz : traceConf ) {
            createLogger(clazz, Level.TRACE)
        }

        if(!consoleAppender)
            logger.debug "Console appender: disabled"
    }

    protected Logger createLogger(String clazz, Level level ) {
        def logger = loggerContext.getLogger( clazz )
        logger.setLevel(level)
        logger.setAdditive(false)
        if( fileAppender )
            logger.addAppender(fileAppender)
        if( consoleAppender )
            logger.addAppender(consoleAppender)
        if( syslogAppender )
            logger.addAppender(syslogAppender)

        return logger
    }

    protected Appender createConsoleAppender() {

        final Appender result = daemon && opts.isBackground() ? null : ( opts.ansiLog ? new CaptureAppender() : new ConsoleAppender())
        if( result )  {
            final filter = new ConsoleLoggerFilter( packages )
            filter.setContext(loggerContext)
            filter.start()

            result.setContext(loggerContext)
            if( result instanceof ConsoleAppender )
                result.setEncoder( new LayoutWrappingEncoder( layout: new PrettyConsoleLayout() ) )
            result.addFilter(filter)
            result.start()
        }

        return result
    }

    protected SyslogAppender createSyslogAppender() {
        SyslogAppender result = syslog ? new SyslogAppender() : null
        if( result )  {
            // Multiple settings were sent or its just the hostname
            final conf = syslog.contains(':') ? syslog.tokenize(':') : [syslog]
            def sysLogHost = conf[0]
            def sysLogPort = conf[1] ? conf[1].toInteger() : 514
            def sysLogFacility = conf[2] ? conf[2] : "LOCAL0"
            result.setSyslogHost(sysLogHost)
            result.setPort(sysLogPort)
            result.setFacility(sysLogFacility)
            result.setSuffixPattern("nextflow: %-5level [%thread] %logger{0} - %msg")
            result.setContext(loggerContext)
            result.setThrowableExcluded(false)
            result.setStackTracePattern('\t')
            result.start()
        }

        return result
    }

    protected RollingFileAppender createRollingAppender() {

        RollingFileAppender result = logFileName ? new RollingFileAppender() : null
        if( result ) {
            result.file = logFileName

            def rollingPolicy = new  FixedWindowRollingPolicy( )
            rollingPolicy.fileNamePattern = "${logFileName}.%i"
            rollingPolicy.setContext(loggerContext)
            rollingPolicy.setParent(result)
            rollingPolicy.setMinIndex(1)
            rollingPolicy.setMaxIndex(9)
            rollingPolicy.start()

            result.rollingPolicy = rollingPolicy
            result.encoder = createEncoder()
            result.setContext(loggerContext)
            result.setTriggeringPolicy(new RollOnStartupPolicy())
            result.triggeringPolicy.start()
            result.start()
        }

        return result
    }

    protected FileAppender createFileAppender() {

        FileAppender result = logFileName ? new FileAppender() : null
        if( result ) {
            result.file = logFileName
            result.encoder = createEncoder()
            result.setContext(loggerContext)
            result.bufferSize = FileSize.valueOf('64KB')
            result.start()
        }

        return result
    }

    protected PatternLayoutEncoder createEncoder() {
        def result = new PatternLayoutEncoder()
        result.setPattern('%d{MMM-dd HH:mm:ss.SSS} [%t] %-5level %logger{36} - %msg%n')
        result.setContext(loggerContext)
        result.start()
        return result
    }

    /**
     * Configure the underlying logging subsystem. The logging information are saved
     * in a file named '.nextflow.log' in the current path, plus to the terminal console.
     * <p>
     *     By default only the INFORMATION level logging messages are visualized to the console,
     *     instead in the file are saved the DEBUG level messages.
     *
     * @param logFileName The file where save the application log
     * @param quiet When {@code true} only Warning and Error messages are visualized to teh console
     * @param debugConf The list of packages for which use a Debug logging level
     * @param traceConf The list of packages for which use a Trace logging level
     */

    static void configureLogger( Launcher launcher ) {
        new LoggerHelper(launcher.options)
                .setDaemon(launcher.isDaemon())
                .setRolling(true)
                .setSyslog(launcher.options.syslog)
                .setup()
    }

    static void configureLogger( final CliOptions opts, boolean daemon = false ) {
        new LoggerHelper(opts)
                .setDaemon(daemon)
                .setRolling(true)
                .setSyslog(opts.syslog)
                .setup()
    }

    /*
     * Filters the logging event based on the level assigned to a specific 'package'
     */
    static class ConsoleLoggerFilter extends Filter<ILoggingEvent> {

        private final Map<String,Level> allLevels
        private final Set<String> packages
        private final int len

        ConsoleLoggerFilter( Map<String,Level> levels )  {
            this.allLevels = levels
            this.packages = levels.keySet()
            this.len = packages.size()
        }

        @Override
        FilterReply decide(ILoggingEvent event) {

            if (!isStarted()) {
                return FilterReply.NEUTRAL;
            }

            def logger = event.getLoggerName()
            def level = event.getLevel()
            for( int i=0; i<len; i++ ) {
                final key=packages[i]
                if ( logger.startsWith(key) && level.isGreaterOrEqual(Level.INFO) && level.isGreaterOrEqual(allLevels[key]) ) {
                    return FilterReply.NEUTRAL
                }
            }

            return FilterReply.DENY
        }
    }

    /*
     * Do not print INFO level prefix, used to print logging information
     * to the application stdout
     */
    static private class PrettyConsoleLayout extends LayoutBase<ILoggingEvent> {

        @Override
        String doLayout(ILoggingEvent event) {
            final session = (Session)Global.session
            fmtEvent(event, session)
        }
    }

    /**
     * Format a logging event
     *
     * @param event The logging event
     * @param session Nextflow session object
     * @return the formatted logging message
     */
    static protected String fmtEvent(ILoggingEvent event, Session session) {
        final buffer = new StringBuilder(512);
        if( event.level == Level.INFO ) {
            buffer .append(event.getFormattedMessage()) .append(CoreConstants.LINE_SEPARATOR)
        }
        else if( event.level == Level.ERROR ) {
            def error = ( event.getThrowableProxy() instanceof ThrowableProxy
                    ? (event.getThrowableProxy() as ThrowableProxy).throwable
                    : null) as Throwable

            appendFormattedMessage(buffer, event, error, session)
        }
        else {
            buffer
                    .append( event.getLevel().toString() ) .append(": ")
                    .append(event.getFormattedMessage())
                    .append(CoreConstants.LINE_SEPARATOR)
        }

        return buffer.toString()
    }

    /**
     * Find out the script line where the error has thrown
     */
    static protected void appendFormattedMessage( StringBuilder buffer, ILoggingEvent event, Throwable fail, Session session) {
        final className = session?.script?.getClass()?.getName()
        final message = event.getFormattedMessage()
        final quiet = fail instanceof AbortOperationException || fail instanceof ProcessException || ScriptRuntimeException
        final normalize = { String str -> str ?. replace("${className}.", '')}
        List error = fail ? findErrorLine(fail, ScriptMeta.allScriptNames()) : null

        // error string is not shown for abort operation
        if( !quiet ) {
            buffer.append("ERROR ~ ")
        }

        // add the log message
        if( fail instanceof MissingPropertyException ) {
            def name = fail.property ?: getDetailMessage(fail)
            buffer.append("No such variable: ${name}")
        }
        else if( fail instanceof MissingMethodException ) {
            buffer.append(getMissingMethodMessage(fail))
        }
        else if( fail instanceof NoSuchFileException ) {
            buffer.append("No such file: ${normalize(fail.message)}")
        }
        else if( message && message.startsWith(STARTUP_ERROR))  {
            buffer.append(formatStartupErrorMessage(message))
        }
        else if( message && !message.startsWith('@') ) {
            buffer.append(normalize(message))
        }
        else if( fail ) {
            buffer.append( normalize(fail.message) ?: "Unexpected error [${fail.class.simpleName}]" )
        }
        else {
            buffer.append("Unexpected error")
        }
        buffer.append(CoreConstants.LINE_SEPARATOR)
        buffer.append(CoreConstants.LINE_SEPARATOR)

        // extra formatting
        if( error ) {
            buffer.append(" -- Check script '${error[0]}' at line: ${error[1]} or see '${logFileName}' file for more details")
            buffer.append(CoreConstants.LINE_SEPARATOR)
        }
        else if( logFileName && !quiet ) {
            buffer.append(" -- Check '${logFileName}' file for details")
            buffer.append(CoreConstants.LINE_SEPARATOR)
        }

    }

    @PackageScope
    static String formatStartupErrorMessage( String message ) {
        message.replace(STARTUP_ERROR,'').replaceFirst(/^_nf_script_[a-z0-9]+: *[0-9]+: */,'')
    }

    static String getDetailMessage(Throwable error) {
        try {
            def clazz = error.class
            while(  clazz != Throwable && clazz )
                clazz = clazz.getSuperclass()
            Field field = clazz.getDeclaredField('detailMessage')
            field.setAccessible(true)
            return field.get(error)
        }
        catch( Throwable e ) {
            return null
        }
    }

    static String getMissingMethodMessage(MissingMethodException error) {
        try {
            return getMissingMethodMessage0(error)
        }
        catch( Throwable e ) {
            return error?.message
        }
    }

    private static String getMissingMethodMessage0(MissingMethodException error) {
        def name = error.getMethod()
        def type = error.getType()
        def args = error.getArguments()

        String left = fmtType(type)

        String right = args?.collect { fmtValue(it) }?.join(', ')

        if( name == 'or' )
            return "Invalid invocation of operator | (pipe) with left operand: $left and right operand: $right"

        String msg
        def hasLeftTarget = !WorkflowBinding.isAssignableFrom(type) && !BaseScript.isAssignableFrom(type)
        def found = type.getMethods().find { it.name == name }
        if( found ) {
            msg = "Invalid method invocation `$name` with arguments: $right"
            if( hasLeftTarget ) msg += " on $left"
        }
        else {
            msg = "Unknown method invocation `$name`"
            if( hasLeftTarget ) msg += " on $left"
            def tips = type.getMethods().collect { it.name }.closest(name)
            if( tips )
                msg += " -- Did you mean?\n" + tips.collect { "  $it"}.join('\n')
        }

        return msg
    }

    static @PackageScope List<String> findErrorLine( Throwable e, Map<String, Path> allNames ) {
        def lines = getErrorLines(e)
        List error = null
        for( String str : lines ) {
            if( (error=getErrorLine(str,allNames))) {
                break
            }
        }
        return error
    }

    static @PackageScope String[] getErrorLines(Throwable e) {
        try {
            return ExceptionUtils.getStackTrace(e).split('\n')
        }
        catch( Throwable t ) {
            log.warn "Oops .. something wrong formatting the error stack trace | ${t.message ?: t}", e
            return Collections.emptyList() as String[]
        }
    }

    static private Pattern ERR_LINE_REGEX = ~/\((Script_[0-9a-f]{8}):(\d*)\)$/

    @PackageScope
    static List<String> getErrorLine( String str, Map<String,Path> allNames ) {
        def m = ERR_LINE_REGEX.matcher(str)
        if( m.find() ) {
            def name = m.group(1)
            def line = m.group(2)
            if( allNames[name] ) {
                final pwd = System.getProperty("user.dir")
                def script = allNames[name].toUriString()
                if( script.startsWith(pwd) ) {
                    script = script.substring(pwd.length())
                    while( script.startsWith('/') )
                        script = script.substring(1)
                }

                return [ script, line ]
            }
        }
        return null
    }

    /**
     * Trigger a new log file on application start-up
     *
     * See here http://stackoverflow.com/a/2647471/395921
     */
    @NoAutoStart
    static class RollOnStartupPolicy<E> extends TriggeringPolicyBase<E> {

        private final AtomicBoolean firstTime = new AtomicBoolean(true);

        @Override
        boolean isTriggeringEvent(File activeFile, E event) {

            if ( !firstTime.get() ) { // fast path
                return false;
            }

            if (firstTime.getAndSet(false) && !FileHelper.empty(activeFile) ) {
                return true;
            }
            return false;
        }

    }

    /**
     * Capture logging events and forward them to
     */
    static private class CaptureAppender extends AppenderBase<ILoggingEvent> {

        @Override
        protected void append(ILoggingEvent event) {
            final session = (Session)Global.session

            try {
                final message = fmtEvent(event, session).trim()
                final renderer = session?.ansiLogObserver
                if( !renderer )
                    System.out.println(message)

                else if( event.marker == STICKY )
                    renderer.appendSticky(message)

                else if( event.level==Level.ERROR )
                    renderer.appendError(message)

                else if( event.level==Level.WARN )
                    renderer.appendWarning(message)

                else
                    renderer.appendInfo(message)
            }
            catch (Throwable e) {
                e.printStackTrace()
            }

        }
    }

    static private char OPEN_CH = '[' as char
    static private char CLOSE_CH = ']' as char
    static private char SLASH_CH = '/' as char
    static private int ZERO_CH = '0' as char
    static private int NINE_CH = '9' as char
    static private int ALPHA_CH = 'a' as char
    static private int EFFE_CH = 'f' as char

    static boolean isHashLogPrefix(String str) {
        if( str?.length()<10 )
            return false
        if( str.charAt(0) != OPEN_CH )
            return false
        if( str.charAt(3) != SLASH_CH )
            return false
        if( str.charAt(10) != CLOSE_CH )
            return false
        if( !isHex(str.charAt(1)) )
            return false
        if( !isHex(str.charAt(2)) )
            return false

        for( int i=4; i<10; i++ ) {
            if( !isHex(str.charAt(i)))
                return false
        }

        return true
    }

    static protected boolean isHex(char ch) {
        (ch >= ZERO_CH && ch <= NINE_CH) || (ch >= ALPHA_CH && ch <= EFFE_CH)
    }


    static String fmtValue( obj ) {
        if( obj instanceof ComponentDef )
            return obj.getSignature()
        if( obj != null )
            return "$obj (${obj.getClass().getName()})"
        return 'null'
    }

    static String fmtType( ChainableDef obj ) {
        if( obj instanceof ComponentDef )
            return obj.getSignature()
        if( obj != null )
            return "$obj.type '$obj.name'"
        else
            return "component null"
    }

    static String fmtType( OpCall obj ) {
        return obj ? "operator '${obj?.getMethodName()}'" : "operator null"
    }

    static String fmtType( Object type ) {
        if( type instanceof Class ) {
            if( DataflowWriteChannel.isAssignableFrom(type) || DataflowReadChannel .isAssignableFrom(type) )
                return 'channel type'
            if( ComponentDef.isAssignableFrom(type) )
                return 'process type'
            if( FunctionDef.isAssignableFrom(type) )
                return 'function type'
            if( WorkflowDef.isAssignableFrom(type) )
                return 'workflow type'
            if( CompositeDef.isAssignableFrom(type) )
                return 'expression type'
            if( OpCall.isAssignableFrom(type) )
                return 'operator type'

            return type.getSimpleName() + ' type'
        }

        if( type instanceof ChainableDef )
            return fmtType((ChainableDef)type)

        if( type instanceof OpCall )
            return fmtType((OpCall)type)

        if( type instanceof DataflowWriteChannel || type instanceof DataflowReadChannel )
            return "channel object"

        type == null ? "null object" : "${type.getClass().getSimpleName()} object"
    }
}
