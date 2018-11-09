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

package nextflow.extension

import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.text.DateFormat
import java.text.SimpleDateFormat
import java.util.concurrent.locks.Lock
import java.util.regex.Pattern

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import nextflow.file.FileHelper
import nextflow.util.CheckHelper
import nextflow.util.Duration
import nextflow.file.FileMutex
import nextflow.util.MemoryUnit
import nextflow.util.RateUnit
import org.apache.commons.lang.StringUtils
import org.codehaus.groovy.runtime.DefaultGroovyMethods
import org.codehaus.groovy.runtime.GStringImpl
import org.codehaus.groovy.runtime.ResourceGroovyMethods
import org.codehaus.groovy.runtime.StringGroovyMethods
import org.slf4j.Logger

/**
 * Generic extensions
 *
 * See more about extension methods
 * http://docs.codehaus.org/display/GROOVY/Creating+an+extension+module
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Bolts {

    static public final String DATETIME_FORMAT = 'dd-MM-yyyy HH:mm'

    static private Pattern PATTERN_RIGHT_TRIM = ~/\s+$/

    static private Pattern PATTERN_LEFT_TRIM = ~/^\s+/

    @Memoized
    static private ThreadLocal<DateFormat> getLocalDateFormat(String fmt, TimeZone tz) {

        return new ThreadLocal<DateFormat>() {
            @Override
            protected DateFormat initialValue() {
                def result = new SimpleDateFormat(fmt)
                if(tz) result.setTimeZone(tz)
                return result
            }
        }
    }

    /**
     * Format a {@link Date} object
     *
     * @param self The {@link Date} object to format
     * @param format The date format to use eg. {@code dd-MM-yyyy HH:mm}.
     * @param tz The timezone to be used eg. {@code UTC}. If {@code null} the current timezone is used.
     * @return The date-time formatted as a string
     */
    static format(Date self, String format=null, String tz=null) {
        TimeZone zone = tz ? TimeZone.getTimeZone(tz) : null
        getLocalDateFormat(format ?: DATETIME_FORMAT, zone).get().format(self)
    }

    /**
     * Format a {@link Date} object
     *
     * @param self The {@link Date} object to format
     * @param format The date format to use eg. {@code dd-MM-yyyy HH:mm}
     * @param tz The timezone to be used. If {@code null} the current timezone is used.
     * @return The date-time formatted as a string
     */
    static format(Date self, String format, TimeZone tz) {
        getLocalDateFormat(format ?: DATETIME_FORMAT, tz).get().format(self)
    }

    static List pairs(Map self, Map opts=null) {
        def flat = opts?.flat == true
        def result = []
        for( Map.Entry entry : self.entrySet() ) {
            if( flat && entry.value instanceof Collection )
                entry.value.iterator().each { result << [entry.key, it] }
            else
                result << [entry.key, entry.value]
        }

        return result
    }

    /**
     * Remove the left side after a dot (including it) e.g.
     * <pre>
     *     0.10     => 0
     *     10000.00 => 10000
     * </pre>
     *
     * @param self
     * @return
     */
    static String trimDotZero(String self) {
        int p = self?.indexOf('.')
        p!=-1 ? self.substring(0,p) : self
    }

    /**
     * Remove blank chars at the end of the string
     *
     * @param self The string itself
     * @return The string with blanks removed
     */

    static String rightTrim(String self) {
        self.replaceAll(PATTERN_RIGHT_TRIM,'')
    }

    /**
     * Remove blank chars at string beginning
     *
     * @param self The string itself
     * @return The string with blanks removed
     */
    static String leftTrim( String self ) {
        self.replaceAll(PATTERN_LEFT_TRIM,'')
    }

    /**
     * <p>Strips any of a set of characters from the start and end of a String.
     * This is similar to {@link String#trim()} but allows the characters
     * to be stripped to be controlled.</p>
     *
     * <p>A <code>null</code> input String returns <code>null</code>.
     * An empty string ("") input returns the empty string.</p>
     *
     * <p>If the stripChars String is <code>null</code>, whitespace is
     * stripped as defined by {@link Character#isWhitespace(char)}.
     * Alternatively use {@link #strip(String)}.</p>
     *
     * <pre>
     * StringUtils.strip(null, *)          = null
     * StringUtils.strip("", *)            = ""
     * StringUtils.strip("abc", null)      = "abc"
     * StringUtils.strip("  abc", null)    = "abc"
     * StringUtils.strip("abc  ", null)    = "abc"
     * StringUtils.strip(" abc ", null)    = "abc"
     * StringUtils.strip("  abcyx", "xyz") = "  abc"
     * </pre>
     *
     * @param str  the String to remove characters from, may be null
     * @param stripChars  the characters to remove, null treated as whitespace
     * @return the stripped String, <code>null</code> if null String input
     */
    static strip( String self, String stripChars = null ) {
        StringUtils.strip(self, stripChars)
    }

    /**
     * <p>Strips any of a set of characters from the start of a String.</p>
     *
     * <p>A <code>null</code> input String returns <code>null</code>.
     * An empty string ("") input returns the empty string.</p>
     *
     * <p>If the stripChars String is <code>null</code>, whitespace is
     * stripped as defined by {@link Character#isWhitespace(char)}.</p>
     *
     * <pre>
     * StringUtils.stripStart(null, *)          = null
     * StringUtils.stripStart("", *)            = ""
     * StringUtils.stripStart("abc", "")        = "abc"
     * StringUtils.stripStart("abc", null)      = "abc"
     * StringUtils.stripStart("  abc", null)    = "abc"
     * StringUtils.stripStart("abc  ", null)    = "abc  "
     * StringUtils.stripStart(" abc ", null)    = "abc "
     * StringUtils.stripStart("yxabc  ", "xyz") = "abc  "
     * </pre>
     *
     * @param str  the String to remove characters from, may be null
     * @param stripChars  the characters to remove, null treated as whitespace
     * @return the stripped String, <code>null</code> if null String input
     */
    static stripStart( String self, String stripChars = null ) {
        StringUtils.stripStart(self, stripChars)
    }

    /**
     * <p>Strips any of a set of characters from the end of a String.</p>
     *
     * <p>A <code>null</code> input String returns <code>null</code>.
     * An empty string ("") input returns the empty string.</p>
     *
     * <p>If the stripChars String is <code>null</code>, whitespace is
     * stripped as defined by {@link Character#isWhitespace(char)}.</p>
     *
     * <pre>
     * StringUtils.stripEnd(null, *)          = null
     * StringUtils.stripEnd("", *)            = ""
     * StringUtils.stripEnd("abc", "")        = "abc"
     * StringUtils.stripEnd("abc", null)      = "abc"
     * StringUtils.stripEnd("  abc", null)    = "  abc"
     * StringUtils.stripEnd("abc  ", null)    = "abc"
     * StringUtils.stripEnd(" abc ", null)    = " abc"
     * StringUtils.stripEnd("  abcyx", "xyz") = "  abc"
     * StringUtils.stripEnd("120.00", ".0")   = "12"
     * </pre>
     *
     * @param str  the String to remove characters from, may be null
     * @param stripChars  the set of characters to remove, null treated as whitespace
     * @return the stripped String, <code>null</code> if null String input
     */
    static stripEnd( String self, String stripChars = null ) {
        StringUtils.stripEnd(self, stripChars)
    }

    /**
     * <p>Capitalizes a String changing the first letter to title case as
     * per {@link Character#toTitleCase(char)}. No other letters are changed.</p>
     *
     * <p>For a word based algorithm, see {@link WordUtils#capitalize(String)}.
     * A <code>null</code> input String returns <code>null</code>.</p>
     *
     * <pre>
     * StringUtils.capitalize(null)  = null
     * StringUtils.capitalize("")    = ""
     * StringUtils.capitalize("cat") = "Cat"
     * StringUtils.capitalize("cAt") = "CAt"
     * </pre>
     *
     */
    static String capitalize(String self) {
        StringUtils.capitalize(self)
    }

    /**
     * <p>Uncapitalizes a String changing the first letter to title case as
     * per {@link Character#toLowerCase(char)}. No other letters are changed.</p>
     *
     * <p>For a word based algorithm, see {@link WordUtils#uncapitalize(String)}.
     * A <code>null</code> input String returns <code>null</code>.</p>
     *
     * <pre>
     * StringUtils.uncapitalize(null)  = null
     * StringUtils.uncapitalize("")    = ""
     * StringUtils.uncapitalize("Cat") = "cat"
     * StringUtils.uncapitalize("CAT") = "cAT"
     * </pre>
     *
     */
    static String uncapitalize(String self) {
        StringUtils.uncapitalize(self)
    }

    /**
     * Check if a alphabetic characters in a string are lowercase. Non alphabetic characters are ingored
     * @param self The string to check
     * @return {@true} if the string contains no uppercase characters, {@code false} otherwise
     */
    static boolean isLowerCase(String self) {
        if( self ) for( int i=0; i<self.size(); i++ ) {
            if( Character.isUpperCase(self.charAt(i)))
                return false
        }
        return true
    }

    /**
     * Check if a alphabetic characters in a string are uppercase. Non alphabetic characters are ignored
     * @param self The string to check
     * @return {@true} if the string contains no lowercase characters, {@code false} otherwise
     */
    static boolean isUpperCase(String self) {
        if( self ) for( int i=0; i<self.size(); i++ ) {
            if( Character.isLowerCase(self.charAt(i)))
                return false
        }
        return true
    }

    /**
     * Check if ALL characters in a string are lowercase.
     * @param self The string to check
     * @return {@true} when all characters are uppercase, {@code false} otherwise
     */
    static boolean isAllLowerCase(String self) {
        StringUtils.isAllLowerCase(self)
    }

    /**
     * Check if ALL characters in a string are uppercase.
     * @param self The string to check
     * @return {@true} when all characters are uppercase, {@code false} otherwise
     */
    static boolean isAllUpperCase(String self) {
        StringUtils.isAllUpperCase(self)
    }

    static private Pattern getPattern( obj ) {

        if( obj instanceof Map ) {
            if( obj.containsKey('pattern') )
                obj = obj.pattern
            else
                return null
        }

        if( obj instanceof Pattern ) {
            return (Pattern)obj
        }

        if( obj instanceof String ) {
            Pattern.compile(Pattern.quote(obj))
        }

        if( obj != null )
            throw new IllegalArgumentException()

        return null
    }



    /**
     * Invokes the specify closure including it with a lock/unlock calls pair
     *
     * @param self
     * @param interruptible
     * @param closure
     * @return the closure result
     */
    static <T> T withLock( Lock self, boolean interruptible = false, Closure<T> closure ) {
        // acquire the lock
        if( interruptible )
            self.lockInterruptibly()
        else
            self.lock()

        try {
            return closure.call()
        }
        finally {
            self.unlock();
        }
    }

    /**
     * Invokes the specify closure only if it is able to acquire a lock
     *
     * @param self
     * @param interruptible
     * @param closure
     * @return the closure result
     */
    static boolean tryLock( Lock self, Closure closure ) {
        if( !self.tryLock() )
            return false

        try {
            closure.call()
        }
        finally {
            self.unlock()
            return true
        }
    }

    /**
     * Creates a file system wide lock that prevent two or more JVM instances/process
     * to work on the same file
     *
     * Note: this does not protected against multiple-thread accessing the file in a
     * concurrent manner.
     *
     * @param
     *      self The file over which define the lock
     * @param
     *      timeout An option timeout elapsed which the a {@link InterruptedException} is thrown
     * @param
     *      closure The action to apply during the lock file spawn
     * @return
     *      The user provided {@code closure} result
     *
     * @throws
     *      InterruptedException if the lock cannot be acquired within the specified {@code timeout}
     */
    static withLock(File self, Duration timeout = null, Closure closure) {
        def locker = new FileMutex(self)
        if( timeout )
            locker.setTimeout(timeout)
        locker.lock(closure)
    }


    /**
     * Creates a file system wide lock that prevent two or more JVM instances/process
     * to work on the same file
     *
     * Note: this does not protected against multiple-thread accessing the file in a
     * concurrent manner.
     *
     * @param
     *      self The file over which define the lock
     * @param
     *      timeout An option timeout elapsed which the a {@link InterruptedException} is thrown
     * @param
     *      closure The action to apply during the lock file spawn
     * @return
     *      The user provided {@code closure} result
     *
     * @throws
     *      InterruptedException if the lock cannot be acquired within the specified {@code timeout}
     */
    static withLock( Path self, Duration timeout, Closure closure ) {
        def locker = new FileMutex(self.toFile())
        if( timeout )
            locker.setTimeout(timeout)
        locker.lock(closure)
    }


    /**
     * Converts a {@code String} to a {@code Duration} object
     *
     * @param self
     * @param type
     * @return
     */
    static def asType( String self, Class type ) {
        if( type == Duration ) {
            return new Duration(self)
        }
        else if( Path.isAssignableFrom(type) ) {
            return FileHelper.asPath(self)
        }
        else if( type == MemoryUnit ) {
            return new MemoryUnit(self)
        }
        else if( type == RateUnit ) {
            return new RateUnit(self)
        }

        StringGroovyMethods.asType(self, type);
    }

    /**
     * Converts a {@code GString} to a {@code Duration} object
     *
     * @param self
     * @param type
     * @return
     */
    static def asType( GString self, Class type ) {
        if( type == Duration ) {
            return new Duration(self.toString())
        }
        else if( Path.isAssignableFrom(type) ) {
            return FileHelper.asPath(self)
        }
        else if( type == MemoryUnit ) {
            return new MemoryUnit(self.toString())
        }

        StringGroovyMethods.asType(self, type);
    }

    /**
     * Converts a {@code Number} to a {@code Duration} object
     *
     * @param self
     * @param type
     * @return
     */
    static def asType( Number self, Class type ) {
        if( type == Duration ) {
            return new Duration(self.longValue())
        }
        if( type == MemoryUnit ) {
            return new MemoryUnit(self.longValue())
        }
        DefaultGroovyMethods.asType(self, type);
    }

    /**
     * Converts a {@code File} to a {@code Path} object
     *
     * @param self
     * @param type
     * @return
     */
    static def asType( File self, Class type ) {
        if( Path.isAssignableFrom(type) ) {
            return self.toPath()
        }

        ResourceGroovyMethods.asType(self, type);
    }


    static def <T> T getOrCreate( Map self, key, factory ) {

        if( self.containsKey(key) )
            return (T)self.get(key)

        synchronized (self) {
            if( self.containsKey(key) )
                return (T)self.get(key)

            def result = factory instanceof Closure ? factory.call(key) : factory
            self.put(key,result)
            return (T)result
        }

    }

    /**
     * Navigate a map of maps traversing multiple attribute using dot notation. For example:
     * {@code x.y.z }
     *
     * @param self The root map object
     * @param key A dot separated list of keys
     * @param closure An optional closure to be applied. Only if all keys exists
     * @return The value associated to the specified key(s) or null on first missing entry
     */
    static def navigate( Map self, String key, Closure closure = null ) {
        assert key
        def items = key.split(/\./)
        def current = self.get(items[0])

        for( int i=1; i<items.length; i++ ) {
            if( current instanceof Map ) {
                if( current.containsKey(items[i]))
                    current = current.get(items[i])
                else
                    return null
            }
            else if( !current ) {
                return null
            }
            else {
                throw new IllegalArgumentException("Cannot navigate map attribute: '$key' -- Content: $self")
            }
        }

        return closure ? closure(current) : current
    }


    static def navigate(Map self, String key, defValue) {
        def result = navigate(self,key)
        return result!=null ? result : defValue
    }

    /**
     * Converts {@code ConfigObject}s to a plain {@code Map}
     *
     * @param config
     * @return A normalized config object
     */
    static Map toMap( ConfigObject config ) {
        assert config != null
        (Map)normalize0((Map)config)
    }

    static ConfigObject toConfigObject(Map self) {

        def result = new ConfigObject()
        self.each { key, value ->
            if( value instanceof Map ) {
                result.put( key, toConfigObject((Map)value) )
            }
            else {
                result.put( key, value )
            }
        }

        return result
    }

    static private normalize0( config ) {

        if( config instanceof Map ) {
            Map result = new LinkedHashMap(config.size())
            config.keySet().each { name ->
                def value = (config as Map).get(name)
                result.put(name, normalize0(value))
            }
            return result
        }
        else if( config instanceof Collection ) {
            List result = new ArrayList(config.size())
            for( entry in config ) {
                result << normalize0(entry)
            }
            return result
        }
        else if( config instanceof GString ) {
            return config.toString()
        }
        else {
            return config
        }
    }

    /**
     * Indent each line in the given test by a specified prefix
     *
     * @param text
     * @param prefix
     * @return The string indented
     */
    public static String indent( String text, String prefix = ' ' ) {
        def result = new StringBuilder()
        def lines = text ? text.readLines() : Collections.emptyList()
        for( int i=0; i<lines.size(); i++ ) {
            result << prefix
            result << lines.get(i)
            result << '\n'
        }
        return result.toString()
    }

    /**
     * Find all the best matches for the given example string in a list of values
     *
     * @param options A list of string
     * @param sample The example string -- cannot be empty
     * @return The list of options that best matches to the specified example -- return an empty list if none match
     */
    static List<String> bestMatches( Collection<String> options, String sample ) {
        assert sample
        assert options

        // Otherwise look for the most similar
        Map<String,Integer> diffs = [:]
        options.each {
            diffs[it] = StringUtils.getLevenshteinDistance(sample, it)
        }

        // sort the Levenshtein Distance and get the fist entry
        def sorted = diffs.sort { Map.Entry<String,Integer> it -> it.value }
        def nearest = (Map.Entry<String,Integer>)sorted.find()
        int min = nearest.value
        int len = sample.length()

        int threshold = len<=3 ? 1 : ( len > 10 ? 5 : Math.floorDiv(len,2))

        List<String> result
        if( min <= threshold ) {
            result = (List<String>)sorted.findAll { it.value==min } .collect { it.key }
        }
        else {
            result = []
        }

        return result

    }

    static boolean isCamelCase(String str) {
        if( !str ) return false
        for( int i=0; i<str.size()-1; i++ )
            if( Character.getType(str.charAt(i)) == Character.LOWERCASE_LETTER && Character.getType(str.charAt(i+1)) == Character.UPPERCASE_LETTER)
                return true

        return false
    }

    /**
     * Clone a {@link Closure} object and set the specified map as delegate object
     * in the resulting closure
     *
     * @param self The {@link Closure} object to be cloned
     * @param binding The delegate object to set in the new closure
     * @return The cloned {@link Closure} object
     */
    static <T extends Closure> T cloneWith( T self, binding ) {

        def copy = (T)self.clone()
        if( binding != null ) {
            copy.setDelegate(binding)
            copy.setResolveStrategy( Closure.DELEGATE_FIRST )
        }

        return copy
    }

    /**
     * Create a copy of a {@link GString} object cloning all values that are instance of {@link Closure}
     *
     * @param self The {@link GString} itself
     * @param binding A {@link Map} object that is set as delegate object in the cloned closure.
     * @return The cloned {@link GString} instance
     */
    static GString cloneWith( GString self, binding ) {

        def values = new Object[ self.valueCount ]

        // clone the gstring setting the delegate for each closure argument
        for( int i=0; i<self.valueCount; i++ ) {
            values[i] = ( self.values[i] instanceof Closure
                    ? cloneWith(self.values[i] as Closure, binding)
                    : self.values[i]
            )
        }

        new GStringImpl(values, self.strings)
    }



    /**
     * Find all the best matches for the given example string in a list of values
     *
     * @param sample The example string -- cannot be empty
     * @param options A list of string
     * @return The list of options that best matches to the specified example -- return an empty list if none match
     */
    @CompileDynamic
    static List<String> closest(Collection<String> options, String sample ) {
        assert sample

        if( !options )
            return Collections.emptyList()

        // Otherwise look for the most similar
        def diffs = [:]
        options.each {
            diffs[it] = StringUtils.getLevenshteinDistance(sample, it)
        }

        // sort the Levenshtein Distance and get the fist entry
        def sorted = diffs.sort { it.value }
        def nearest = sorted.find()
        def min = nearest.value
        def len = sample.length()

        def threshold = len<=3 ? 1 : ( len > 10 ? 5 : Math.floor(len/2))

        def result
        if( min <= threshold ) {
            result = sorted.findAll { it.value==min } .collect { it.key }
        }
        else {
            result = []
        }

        return result
    }

    private static HashMap<Object,Long> LOGGER_CACHE = new LinkedHashMap<Object,Long>() {
        protected boolean removeEldestEntry(Map.Entry<Object, Long> eldest) {
            return size() > 10_000
        }
    }

    private static final Duration LOG_DFLT_THROTTLE = Duration.of('1min')

    static synchronized private checkLogCache( Object msg, Map params, Closure action ) {

        // -- check if this message has already been printed
        final String str = msg.toString()
        final Throwable error = params?.causedBy as Throwable
        final Duration throttle = params?.throttle as Duration ?: LOG_DFLT_THROTTLE
        final firstOnly = params?.firstOnly == true
        final key = params?.cacheKey ?: str

        long now = System.currentTimeMillis()
        Long ts = LOGGER_CACHE.get(key)
        if( ts && (now - ts <= throttle.toMillis() || firstOnly) ) {
            return
        }
        LOGGER_CACHE.put(key, now)

        action.call(str, error)
    }

    private static Map<String,?> LOGGER_PARAMS = [ cacheKey: Object, causedBy: Throwable, throttle: [String, Number, Duration], firstOnly: Boolean ]


    /**
     * Append a `trace` level entry in the application log.
     *
     * @param log
     *          The {@link Logger} object
     * @param params
     *          Optional named parameters
     *          - `causedBy`: A {@link Throwable} object that raised the error
     *          - `throttle`: When specified suppress identical logs within the specified time {@link Duration}
     * @param msg
     *          The message to print
     *
     */
    static void trace1(Logger log, Map params=null, Object msg ) {
        CheckHelper.checkParams('trace1', params, LOGGER_PARAMS)
        if( !log.isTraceEnabled() || msg==null ) return
        checkLogCache(msg,params) { String str, Throwable t -> t ? log.trace(str,t) : log.trace(str) }
    }

    /**
     * Append a `debug` level entry in the application log.
     *
     * @param log
     *          The {@link Logger} object
     * @param params
     *          Optional named parameters
     *          - `causedBy`: A {@link Throwable} object that raised the error
     *          - `throttle`: When specified suppress identical logs within the specified time {@link Duration}
     * @param msg
     *          The message to print
     *
     */
    static void debug1(Logger log, Map params=null, Object msg ) {
        CheckHelper.checkParams('debug1', params, LOGGER_PARAMS)
        if( !log.isDebugEnabled() || msg==null ) return
        checkLogCache(msg,params) { String str, Throwable t -> t ? log.debug(str,t) : log.debug(str) }
    }

    /**
     * Append a `info` level entry in the application log.
     *
     * @param log
     *          The {@link Logger} object
     * @param params
     *          Optional named parameters
     *          - `causedBy`: A {@link Throwable} object that raised the error
     *          - `throttle`: When specified suppress identical logs within the specified time {@link Duration}
     * @param msg
     *          The message to print
     *
     */
    static void info1(Logger log, Map params=null, Object msg ) {
        CheckHelper.checkParams('info1', params, LOGGER_PARAMS)
        if( !log.isInfoEnabled() || msg==null ) return
        checkLogCache(msg,params) { String str, Throwable t -> t ? log.info(str,t) : log.info(str) }
    }

    /**
     * Append a `warn` level entry in the application log.
     *
     * @param log
     *          The {@link Logger} object
     * @param params
     *          Optional named parameters
     *          - `causedBy`: A {@link Throwable} object that raised the error
     *          - `throttle`: When specified suppress identical logs within the specified time {@link Duration}
     * @param msg
     *          The message to print
     *
     */
    static void warn1(Logger log, Map params=null, Object msg ) {
        CheckHelper.checkParams('warn1', params, LOGGER_PARAMS)
        if( !log.isWarnEnabled() || msg==null ) return
        checkLogCache(msg,params) { String str, Throwable t -> t ? log.warn(str,t) : log.warn(str) }
    }

    /**
     * Append a `error` level entry in the application log.
     *
     * @param log
     *          The {@link Logger} object
     * @param params
     *          Optional named parameters
     *          - `causedBy`: A {@link Throwable} object that raised the error
     *          - `throttle`: When specified suppress identical logs within the specified time {@link Duration}
     * @param msg
     *          The message to print
     *
     */
    static void error1(Logger log, Map params=null, Object msg ) {
        CheckHelper.checkParams('error1', params, LOGGER_PARAMS)
        if( !log.isErrorEnabled() || msg==null ) return
        checkLogCache(msg,params) { String str, Throwable t -> t ? log.error(str,t) : log.error(str) }
    }

    static void trace(Logger log, Object msg) {
        if( log.isTraceEnabled() ) {
            log.trace(msg.toString())
        }
    }

    static void trace(Logger log, Object msg, Throwable e) {
        if( log.isTraceEnabled() ) {
            log.trace(msg.toString(),e)
        }
    }

    static String getErrMessage(Throwable e) {
        if( e instanceof NoSuchFileException ) {
            return "No such file: $e.message"
        }

        return e.message ?: e.toString()
    }
}
