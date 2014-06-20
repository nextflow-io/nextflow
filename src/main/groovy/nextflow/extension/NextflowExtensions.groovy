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

package nextflow.extension

import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.atomic.AtomicLong
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock
import java.util.regex.Pattern

import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import groovyx.gpars.agent.Agent
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import groovyx.gpars.dataflow.stream.DataflowStreamReadAdapter
import groovyx.gpars.group.DefaultPGroup
import groovyx.gpars.group.PGroup
import groovyx.gpars.scheduler.Pool
import nextflow.splitter.BytesSplitter
import nextflow.splitter.FastaSplitter
import nextflow.splitter.StringSplitter
import nextflow.splitter.TextSplitter
import nextflow.util.Duration
import nextflow.util.FileHelper
import org.apache.commons.lang.StringUtils
import org.codehaus.groovy.runtime.DefaultGroovyMethods
import org.codehaus.groovy.runtime.ResourceGroovyMethods
import org.codehaus.groovy.runtime.StringGroovyMethods
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Provides extension methods to chunk text and file
 *
 * See more about extension methods
 * http://docs.codehaus.org/display/GROOVY/Creating+an+extension+module
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NextflowExtensions {

    static private final Logger log = LoggerFactory.getLogger(NextflowExtensions)

    static private Pattern PATTERN_RIGHT_TRIM = ~/\s+$/

    static private Pattern PATTERN_LEFT_TRIM = ~/^\s+/


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
     * <p>For a word based algorithm, see {@link org.apache.commons.lang.WordUtils#capitalize(String)}.
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
     * <p>For a word based algorithm, see {@link org.apache.commons.lang.WordUtils#uncapitalize(String)}.
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


    static private Pattern getPattern( obj ) {

        if( obj instanceof Map ) {
            if( obj.containsKey('pattern') )
                obj = obj.pattern
            else
                return null
        }

        if( obj instanceof Pattern ) {
            return obj
        }

        if( obj instanceof String ) {
            Pattern.compile(Pattern.quote(obj))
        }

        if( obj != null )
            throw new IllegalArgumentException()

        return null
    }

    /// chopLines

    /**
     * Chop a text object by splitting into lines.
     * <p>
     *    The number of lines in each chunk is defined by an optional {@code count} parameter (default: 1).
     * <p>
     *    The lines can be 'forward' to a {@code Collection} or a Dataflow channel using the optional {@code into} named parameter.
     *    For example:
     *    <pre>
     *    def text = ..
     *    def list = text.chopLines( into: [] )
     *
     *    def channel = text.chopLines( count: 2, into: new DataflowQueue() )
     *    </pre>
     *
     *
     * @param obj The target text object, it can be an instance of the following classes:
     *      <li>{@code CharSequence}
     *      <li>{@code Path}
     *      <li>{@code File}
     *      <li>{@code Reader}
     *      <li>{@code InputStream}
     *
     *
     * @param options
     *      Optional named parameters:
     *      <li>{@code count}: The number of lines in each chopped text block. default: 1
     *      <li>{@code into}: The target container to which are sent the chopped lines.
     *          It can be an instance of {@code Collection} or {@code DataflowQueue}
     *      <li>{@code charset}: The {@code Charset} to use when reading a binary file
     *
     * @return
     *      The container object specified by using the {@code into} parameter, or the last
     *      chopped line(s) if container object has been specified.
     *
     */
    @Deprecated
    static chopLines( Object obj, Map options = [:] ) {
        chopLines(obj,options,null)
    }

    /**
     * Chop a text object by splitting into lines.
     * <p>
     *    The number of lines in each chunk is defined by an optional {@code count} parameter (default: 1).
     * <p>
     *    The lines can be 'forward' to a {@code Collection} or a Dataflow channel using the optional {@code into} named parameter.
     *    For example:
     *    <pre>
     *    def text = ..
     *    def list = text.chopLines( into: [] )
     *
     *    def channel = text.chopLines( count: 2, into: new DataflowQueue() )
     *    </pre>
     *
     *
     * @param obj The target text object, it can be an instance of the following classes:
     *      <li>{@code CharSequence}
     *      <li>{@code Path}
     *      <li>{@code File}
     *      <li>{@code Reader}
     *      <li>{@code InputStream}
     *
     *
     * @param options
     *      Optional named parameters:
     *      <li>{@code count}: The number of lines in each chopped text block. default: 1
     *      <li>{@code into}: The target container to which are sent the chopped lines.
     *          It can be an instance of {@code Collection} or {@code DataflowQueue}
     *      <li>{@code charset}: The {@code Charset} to use when reading a binary file
     *
     * @param closure
     *      For each chopped line(s) block this closure is invoked receiving the lines as
     *      argument.
     *      If the closure declares a second argument, it will contain the current line (or lines block) index
     *      The closure return value is added to the container specified by the {@code into} parameter,
     *      otherwise is ignored.
     *
     *
     * @return
     *      The container object specified by using the {@code into} parameter, or the last
     *      chopped line(s) if container object has been specified.
     *
     */
    @Deprecated
    static chopLines( Object obj, Map options = [:], Closure closure ) {
        log.warn "Method 'chopLines' has been deprecated and will be removed in a future release -- consider using 'splitText' instead"
        def opt = options != null ? options : [:]
        if( closure ) opt.each = closure
        new TextSplitter().options(options).target(obj).split()

    }

    /// chopString

    /**
     * Chop a string value by splitting into sub-strings having the same length.
     * <p>
     *    The number of characters in each piece is defined by an optional {@code count} parameter (default: 1).
     * <p>
     *    The pieces can be sent to a {@code Collection} or a Dataflow channel using the optional {@code into} named parameter.
     *    For example:
     *    <pre>
     *    def text = ..
     *    def list = text.chopString( into: [] )
     *
     *    def channel = text.chopString( count: 2, into: new DataflowQueue() )
     *    </pre>
     *
     *
     * @param obj The target string object, it can be an instance of the following classes:
     *      <li>{@code CharSequence}
     *      <li>{@code Path}
     *      <li>{@code File}
     *      <li>{@code Reader}
     *      <li>{@code InputStream}
     *
     *
     * @param options
     *      Optional named parameters:
     *      <li>{@code count}: The number of characters in each substring. default: 1
     *      <li>{@code into}: The target container to which are sent the chopped sub-strings.
     *          It can be an instance of {@code Collection} or {@code DataflowQueue}
     *      <li>{@code charset}: The {@code Charset} to use when reading a binary file
     *
     * @return
     *      The container object specified by using the {@code into} parameter, or the last
     *      chopped sub-string if container object has been specified.
     *
     */
    @Deprecated
    static chopString( Object object, Map options = [:] ) {
        chopString(object,options,null)
    }

    /**
     * Chop a string value by splitting into substrings having the same length.
     * <p>
     *    The number of characters in each piece is defined by an optional {@code count} parameter (default: 1).
     * <p>
     *    The pieces can be sent to a {@code Collection} or a Dataflow channel using the optional {@code into} named parameter.
     *    For example:
     *    <pre>
     *    def text = ..
     *    def list = text.chopString( into: [] )
     *
     *    def channel = text.chopString( count: 2, into: new DataflowQueue() )
     *    </pre>
     *
     *
     * @param obj The target string object, it can be an instance of the following classes:
     *      <li>{@code CharSequence}
     *      <li>{@code Path}
     *      <li>{@code File}
     *      <li>{@code Reader}
     *      <li>{@code InputStream}
     *
     *
     * @param options
     *      Optional named parameters:
     *      <li>{@code count}: The number of characters in each substring. default: 1
     *      <li>{@code into}: The target container to which are sent the chopped sub-strings.
     *          It can be an instance of {@code Collection} or {@code DataflowQueue}
     *      <li>{@code charset}: The {@code Charset} to use when reading a binary file
     *
     * @param closure
     *      For each chopped sub-string this closure is invoked receiving it as argument.
     *      If the closure declares a second argument, it will contain the current substring index
     *      The closure return value is added to the container specified by the {@code into} parameter,
     *      otherwise is ignored.
     *
     * @return
     *      The container object specified by using the {@code into} parameter, or the last
     *      chopped sub-string if container object has been specified.
     *
     */
    @Deprecated
    static chopString( Object obj, Map options, Closure closure ) {
        log.warn "Method 'chopString' has been deprecated and will be removed in a future release -- consider using 'splitString' instead"

        def opt = options != null ? options : [:]
        if( closure ) opt.each = closure
        new StringSplitter().options(options).target(obj).split()

    }


    ///  chopBytes

    /**
     * Chop a generic value by splitting into ${code byte[]} objects having the same length.
     * <p>
     *    The number of bytes in each array is defined by an optional {@code count} parameter (default: 1).
     * <p>
     *    The arrays can be sent to a {@code Collection} or a Dataflow channel by specifying the optional {@code into}
     *    named parameter.
     *    For example:
     *    <pre>
     *    def data = ..
     *    def list = data.chopBytes( into: [] )
     *
     *    def channel = data.chopBytes( count: 2, into: new DataflowQueue() )
     *    </pre>
     *
     *
     * @param obj The target data object, it can be an instance of the following classes:
     *      <li>{@code CharSequence}
     *      <li>{@code Path}
     *      <li>{@code File}
     *      <li>{@code InputStream}
     *
     * @param options
     *      Optional named parameters:
     *      <li>{@code count}: The number of byte in each byte array piece. default: 1
     *      <li>{@code into}: The target container to which are sent the chopped byte array .
     *          It can be an instance of {@code Collection} or {@code DataflowQueue}
     *      <li>{@code charset}: The {@code Charset} to use when reading a binary file
     *
     * @return
     *      The container object specified by using the {@code into} parameter, or the last
     *      byte array produced if container object has been specified.
     *
     */
    @Deprecated
    static chopBytes( Object obj, Map options = [:]) {
        chopBytes(obj,options,null)
    }

    /**
     * Chop a generic value by splitting into ${code byte[]} objects having the same length.
     * <p>
     *    The number of bytes in each array is defined by an optional {@code count} parameter (default: 1).
     * <p>
     *    The arrays can be sent to a {@code Collection} or a Dataflow channel by specifying the optional {@code into}
     *    named parameter.
     *    For example:
     *    <pre>
     *    def data = ..
     *    def list = data.chopBytes( into: [] )
     *
     *    def channel = data.chopBytes( count: 2, into: new DataflowQueue() )
     *    </pre>
     *
     *
     * @param obj The target data object, it can be an instance of the following classes:
     *      <li>{@code CharSequence}
     *      <li>{@code Path}
     *      <li>{@code File}
     *      <li>{@code InputStream}
     *
     * @param options
     *      Optional named parameters:
     *      <li>{@code count}: The number of byte in each byte array piece. default: 1
     *      <li>{@code into}: The target container to which are sent the chopped byte array .
     *          It can be an instance of {@code Collection} or {@code DataflowQueue}
     *      <li>{@code charset}: The {@code Charset} to use when reading a binary file
     *
     * @param closure
     *      For each byte array produced this closure is invoked receiving it as argument.
     *      If the closure declares a second argument, it will contain the current object index.
     *      The closure return value is added to the container specified by the {@code into} parameter,
     *      otherwise is ignored.
     *
     * @return
     *      The container object specified by using the {@code into} parameter, or the last
     *      byte array produced if container object has been specified.
     *
     */
    @Deprecated
    static chopBytes( Object obj, Map options, Closure closure ) {
        log.warn "Method 'chopBytes' has been deprecated and will be removed in a future release -- consider using 'splitBytes' instead"

        def opt = options != null ? options : [:]
        if( closure ) opt.each = closure
        new BytesSplitter().options(options).target(obj).split()

    }



    // FASTA

    /**
     * Chop a FASTA formatted text by splitting into chunks containing one or more *sequences*.
     * <p>
     *    The number of *sequences* in each chunk is defined by an optional {@code count} parameter (default: 1).
     * <p>
     *    Chunks can be sent to a {@code Collection} or a Dataflow channel by specifying the optional {@code into}
     *    named parameter.
     *    For example:
     *    <pre>
     *    def text = ..
     *    def list = text.chopFasta( into: [] )
     *
     *    def channel = text.chopFasta( count: 2, into: new DataflowQueue() )
     *    </pre>
     *
     *    <p>
     *    By default the sequence is returned as string object. When the *record* option is specified, the sequence
     *    is returned as map object. For example:
     *
     *    <pre>
     *    def text = ..
     *    def list = text.chopFasta( into: [], record: [id: true, seq: true] )
     *    </pre>
     *    The resulting {@code list} contains {@code Map} objects with two entries: the fasta id and the sequence string
     *
     *    <p>
     *    Option *count* cannot be greater than 1 when *record* is used.
     *
     *
     *
     * @param obj The target data object, it can be an instance of the following classes:
     *      <li>{@code CharSequence}
     *      <li>{@code Path}
     *      <li>{@code File}
     *      <li>{@code InputStream}
     *
     * @param options
     *      Optional named parameters:
     *      <li>{@code count}: The number of byte in each byte array piece. default: 1
     *      <li>{@code into}: The target container to which are sent the chopped byte array .
     *          It can be an instance of {@code Collection} or {@code DataflowQueue}
     *      <li>{@code charset}: The {@code Charset} to use when reading a binary file
     *      <li>{@code record}: Specifying the record option the sequence is returned
     *              as a map object instead of a string. This parameter is a map which defines which
     *              fields have to be included in the resulting map object. See {@code #parseFastaRecord}
     *              for more details
     *
     * @param closure
     *      For each sequence  produced this closure is invoked receiving it as argument.
     *      If the closure declares a second argument, it will contain the current sequence index.
     *      The closure return value is added to the container specified by the {@code into} parameter,
     *      otherwise is ignored.
     *
     * @return
     *      The container object specified by using the {@code into} parameter, or the last
     *      byte array produced if container object has been specified.
     *
     */

    @Deprecated
    static chopFasta( Object obj, Map options = [:] ) {
        chopFasta( obj, options, null )
    }

    /**
     * Chop a FASTA formatted text by splitting into chunks containing one or more *sequences*.
     * <p>
     *    The number of *sequences* in each chunk is defined by an optional {@code count} parameter (default: 1).
     * <p>
     *    Chunks can be sent to a {@code Collection} or a Dataflow channel by specifying the optional {@code into}
     *    named parameter.
     *    For example:
     *    <pre>
     *    def text = ..
     *    def list = text.chopFasta( into: [] )
     *
     *    def channel = text.chopFasta( count: 2, into: new DataflowQueue() )
     *    </pre>
     *
     *    <p>
     *    By default the sequence is returned as string object. When the *record* option is specified, the sequence
     *    is returned as map object. For example:
     *
     *    <pre>
     *    def text = ..
     *    def list = text.chopFasta( into: [], record: [id: true, seq: true] )
     *    </pre>
     *    The resulting {@code list} contains {@code Map} objects with two entries: the fasta id and the sequence string
     *
     *    <p>
     *    Option *count* cannot be greater than 1 when *record* is used.
     *
     *
     *
     * @param obj The target data object, it can be an instance of the following classes:
     *      <li>{@code CharSequence}
     *      <li>{@code Path}
     *      <li>{@code File}
     *      <li>{@code InputStream}
     *
     * @param options
     *      Optional named parameters:
     *      <li>{@code count}: The number of byte in each byte array piece. default: 1
     *      <li>{@code into}: The target container to which are sent the chopped byte array .
     *          It can be an instance of {@code Collection} or {@code DataflowQueue}
     *      <li>{@code charset}: The {@code Charset} to use when reading a binary file
     *      <li>{@code record}: Specifying the record option the sequence is returned
     *              as a map object instead of a string. This parameter is a map which defines which
     *              fields have to be included in the resulting map object. See {@code #parseFastaRecord}
     *              for more details
     *
     * @param closure
     *      For each sequence produced this closure is invoked receiving it as argument.
     *      If the closure declares a second argument, it will contain the current sequence index.
     *      The closure return value is added to the container specified by the {@code into} parameter,
     *      otherwise is ignored.
     *
     * @return
     *      The container object specified by using the {@code into} parameter, or the last
     *      byte array produced if container object has been specified.
     *
     * @see #parseFastaRecord(java.lang.CharSequence, java.util.Map)
     * @see DataflowExtensions#chopFasta(groovyx.gpars.dataflow.DataflowReadChannel)
     *
     */
    @Deprecated
    static chopFasta( Object obj, Map options = [:], Closure closure ) {
        log.warn "Method 'chopFasta' has been deprecated and will be removed in a future release -- consider using 'splitFasta' instead"

        def opt = options != null ? options : [:]
        if( closure ) opt.each = closure
        new FastaSplitter().options(options).target(obj).split()

    }



    /**
     * Splits a {@code CharSequence} in text chunks having the specified number of lines
     *
     * @param sequence A sequence of chars to which apply the chunkLines operation
     * @param n The number of lines in each 'chunk'
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */
    @Deprecated
    public static void chunkLines( CharSequence sequence, int n = 1, Closure block ) {
        assert sequence != null
        chunkLines( new StringReader(sequence.toString()), [size: n], block )
    }

    /**
     * Splits a {@code CharSequence} in text chunks having the specified number of lines
     *
     * @param sequence A sequence of chars to which apply the chunkLines operation
     * @param options Specifies the chunk operation options. The following attributes are supported:
     *   <li>{@code size}: The number of lines in each text chunk (default: 1)
     *   <li>{@code autoClose}: Close automatically any channel eventually specified in the chunk closure (default: true)
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */
    @Deprecated
    public static void chunkLines( CharSequence sequence, Map options, Closure block ) {
        assert sequence != null
        chunkLines( new StringReader(sequence.toString()), options, block )
    }

    /**
     * Splits a {@code Reader} in text chunks having the specified number of lines
     *
     * @param reader A reader to which apply the chunkLines operation
     * @param n The number of lines in each 'chunk'
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */
    @Deprecated
    static void chunkLines( Reader reader, int n = 1,  Closure block) {
        chunkLines( reader, [size: n], block)
    }


    /**
     * Splits a {@code Reader} in text chunks having the specified number of lines
     *
     * @param reader A reader to which apply the chunkLines operation
     * @param options Specifies the chunk operation options. The following attributes are supported:
     *   <li>{@code size}: The number of lines in each text chunk (default: 1)
     *   <li>{@code autoClose}: Close automatically any channel eventually specified in the chunk closure (default: true)
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */

    @Deprecated
    static void chunkLines( Reader reader, Map options,  Closure block) {
        log.warn "Method 'chunkLines' has been deprecated and will be removed in a future release -- consider using 'splitText' instead"

        assert reader != null
        assert options != null
        assert block != null

        log.trace "Chunk options: ${options}"

        int size = options?.size ?: 1
        log.trace "Chunk size: $size"

        BufferedReader reader0 = reader instanceof BufferedReader ? reader : new BufferedReader(reader)
        try {
            // -- wrap the owner to intercept any reference to an external dataflow instance
            final interceptor = new WritableChannelInterceptor(block)

            String line
            StringBuilder buffer = new StringBuilder()
            int c=0
            while( (line = reader0.readLine()) != null ) {
                if ( c ) buffer << '\n'
                buffer << line
                if ( ++c == size ) {
                    c = 0
                    block.call( buffer.toString() )

                    buffer.setLength(0)
                }
            }

            if ( buffer.size() ) {
                block.call( buffer.toString() )
            }

            // send a poison pill to any written channel
            def close = options.autoClose
            if( close == null || close == true ) {
                interceptor.closeChannels()
            }
            else {
                log.debug "Skipping channel autoClose"
            }
        }
        finally {
            reader0.closeQuietly()
        }

    }


    /**
     * Splits a {@code File} in text chunks having the specified number of lines
     *
     * @param file A file to which apply the chunkLines operation
     * @param n The number of lines in each 'chunk'
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */
    @Deprecated
    static void chunkLines( File file, int n = 1, Closure block ) {
        assert file
        chunkLines( new FileReader(file), [size: n], block )
    }

    /**
     * Splits a {@code File} in text chunks having the specified number of lines
     *
     * @param file A file to which apply the chunkLines operation
     * @param options Specifies the chunk operation options. The following attributes are supported:
     *   <li>{@code size}: The number of lines in each text chunk (default: 1)
     *   <li>{@code autoClose}: Close automatically any channel eventually specified in the chunk closure (default: true)
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */
    @Deprecated
    static void chunkLines( File file, Map options, Closure block ) {
        assert file
        chunkLines( new FileReader(file), options, block )
    }

    @Deprecated
    static void chunkLines( Path file, int n = 1, Closure block ) {
        assert file
        chunkLines( Files.newBufferedReader(file, Charset.defaultCharset()), n, block )
    }


    @Deprecated
    static void chunkLines( Path file, Map options, Closure block ) {
        assert file
        chunkLines( Files.newBufferedReader(file, Charset.defaultCharset()), options, block )
    }

    /**
     * Splits a {@code InputStream} in text chunks having the specified number of lines
     *
     * @param stream The stream of chars to which apply the chunk operation
     * @param n The number of lines in each 'chunk'
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */
    @Deprecated
    static void chunkLines( InputStream stream, int n = 1, Closure block ) {
        assert stream
        chunkLines( new InputStreamReader(stream), [size: n], block )
    }

    /**
     * Splits a {@code InputStream} in text chunks having the specified number of lines
     *
     * @param stream A text stream to which apply the chunkLines operation
     * @param options Specifies the chunk operation options. The following attributes are supported:
     *   <li>{@code size}: The number of lines in each text chunk (default: 1)
     *   <li>{@code autoClose}: Close automatically any channel eventually specified in the chunk closure (default: true)
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */

    @Deprecated
    static void chunkLines( InputStream stream, Map options, Closure block ) {
        assert stream
        chunkLines( new InputStreamReader(stream), options, block )
    }

    /**
     * Splits a text formatted in multi-FASTA format in chunks containing the specified number of 'sequences'
     * <p>
     * Read more about the FASTA format http://en.wikipedia.org/wiki/FASTA_format
     *
     * @param text A {@code CharSequence} holding the sequences to be splitted
     * @param n The number of 'sequences' in each text chunk (default: 1)
     * @param block A closure invoked for each chunk of sequences, it takes a single (optional) string parameter which represents the current chunk os sequences
     */
    @Deprecated
    public static void chunkFasta( CharSequence text, int n = 1, Closure block ) {
        assert text != null
        chunkFasta( new StringReader(text.toString()), [size: n], block )
    }

    /**
     * Splits a text formatted in multi-FASTA format in chunks containing the specified number of 'sequences'
     * <p>
     * Read more about the FASTA format http://en.wikipedia.org/wiki/FASTA_format
     *
     * @param text A {@code CharSequence} holding the sequences to be splitted
     * @param options Specifies the chunk operation options. The following attributes are supported:
     *   <li>{@code size}: The number of sequences in each chunk (default: 1)
     *   <li>{@code autoClose}: Close automatically any channel eventually specified in the chunk closure (default: true)
     * @param block A closure invoked for each chunk of sequences, it takes a single (optional) string parameter which represents the current chunk os sequences
     */

    @Deprecated
    public static void chunkFasta( CharSequence text, Map options, Closure block ) {
        assert text != null
        chunkFasta( new StringReader(text.toString()), options, block )
    }

    /**
     * Splits a text formatted in multi-FASTA format in chunks containing the specified number of 'sequences'
     * <p>
     * Read more about the FASTA format http://en.wikipedia.org/wiki/FASTA_format
     *
     * @param text A {@code InputStream} holding the sequences to be splitted
     * @param n The number of 'sequences' in each text chunk (default: 1)
     * @param block A closure invoked for each chunk of sequences, it takes a single (optional) string parameter which represents the current chunk os sequences
     */
    @Deprecated
    static void chunkFasta( InputStream text, int n = 1, Closure block ) {
        assert text != null
        chunkFasta( new InputStreamReader(text), [size: n], block)
    }

    /**
     * Splits a text formatted in multi-FASTA format in chunks containing the specified number of 'sequences'
     * <p>
     * Read more about the FASTA format http://en.wikipedia.org/wiki/FASTA_format
     *
     * @param text A {@code InputStream} holding the sequences to be splitted
     * @param options Specifies the chunk operation options. The following attributes are supported:
     *   <li>{@code size}: The number of sequences in each chunk (default: 1)
     *   <li>{@code autoClose}: Close automatically any channel eventually specified in the chunk closure (default: true)
     * @param block A closure invoked for each chunk of sequences, it takes a single (optional) string parameter which represents the current chunk os sequences
     */
    @Deprecated
    static void chunkFasta( InputStream text, Map options, Closure block ) {
        assert text != null
        chunkFasta( new InputStreamReader(text), options, block)
    }

    /**
     * Splits a text formatted in multi-FASTA format in chunks containing the specified number of 'sequences'
     * <p>
     * Read more about the FASTA format http://en.wikipedia.org/wiki/FASTA_format
     *
     * @param text A {@code File} holding the sequences to be splitted
     * @param n The number of 'sequences' in each text chunk (default: 1)
     * @param block A closure invoked for each chunk of sequences, it takes a single (optional) string parameter which represents the current chunk os sequences
     */
    @Deprecated
    static void chunkFasta( File file, int n = 1, Closure block ) {
        assert file
        chunkFasta( new FileReader(file), n, block )
    }

    @Deprecated
    static void chunkFasta( Path path, int n = 1, Closure block ) {
        assert path
        chunkFasta( Files.newBufferedReader(path, Charset.defaultCharset()), n, block )
    }


    /**
     * Splits a text formatted in multi-FASTA format in chunks containing the specified number of 'sequences'
     * <p>
     * Read more about the FASTA format http://en.wikipedia.org/wiki/FASTA_format
     *
     * @param text A {@code File} holding the sequences to be splitted
     * @param options Specifies the chunk operation options. The following attributes are supported:
     *   <li>{@code size}: The number of sequences in each chunk (default: 1)
     *   <li>{@code autoClose}: Close automatically any channel eventually specified in the chunk closure (default: true)
     * @param block A closure invoked for each chunk of sequences, it takes a single (optional) string parameter which represents the current chunk os sequences
     */
    @Deprecated
    static void chunkFasta( File file, Map options, Closure block ) {
        assert file
        chunkFasta( new FileReader(file), options, block )
    }

    @Deprecated
    static void chunkFasta( Path file, Map options, Closure block ) {
        assert file
        chunkFasta( Files.newBufferedReader(file, Charset.defaultCharset()), options, block )
    }


    /**
     * Splits a text formatted in multi-FASTA format in chunks containing the specified number of 'sequences'
     * <p>
     * Read more about the FASTA format http://en.wikipedia.org/wiki/FASTA_format
     *
     * @param text A {@code Reader} holding the sequences to be splitted
     * @param n The number of 'sequences' in each text chunk (default: 1)
     * @param block A closure invoked for each chunk of sequences, it takes a single (optional) string parameter which represents the current chunk os sequences
     */
    @Deprecated
    static void chunkFasta( Reader text, int n = 1,  Closure block ) {
        chunkFasta( text, [size: n], block)
    }

    /**
     * Splits a text formatted in multi-FASTA format in chunks containing the specified number of 'sequences'
     * <p>
     * Read more about the FASTA format http://en.wikipedia.org/wiki/FASTA_format
     *
     * @param text A {@code Reader} holding the sequences to be splitted
     * @param options Specifies the chunk operation options. The following attributes are supported:
     *   <li>{@code size}: The number of sequences in each chunk (default: 1)
     *   <li>{@code autoClose}: Close automatically any channel eventually specified in the chunk closure (default: true)
     * @param block A closure invoked for each chunk of sequences, it takes a single (optional) string parameter which represents the current chunk os sequences
     */
    @Deprecated
    static void chunkFasta( Reader text, Map options, Closure block ) {
        log.warn "Method 'chunkFasta' has been deprecated and will be removed in a future release -- consider using 'splitFasta' instead"

        assert text != null
        assert block != null
        assert options != null

        log.trace "Chunk options: $options"

        int size = options.size ?: 1
        log.trace "Chunk size: $size"

        BufferedReader reader0 = text instanceof BufferedReader ? text : new BufferedReader(text)
        try {

            // -- wrap the owner to intercept any reference to an external dataflow instance
            final interceptor = new WritableChannelInterceptor(block)

            String line
            StringBuilder buffer = new StringBuilder()
            int blockCount=0
            boolean openBlock = false
            while( (line = reader0.readLine()) != null ) {

                if ( line == '' ) {
                    buffer << '\n'
                }
                else if ( !openBlock && line.charAt(0)=='>' ) {
                    openBlock = true
                    buffer << line << '\n'
                }
                else if ( openBlock && line.charAt(0)=='>') {
                    // another block is started

                    if ( ++blockCount == size ) {
                        // invoke the closure, passing the read block as parameter
                        block.call(buffer.toString())

                        buffer.setLength(0)
                        blockCount=0
                    }

                    buffer << line << '\n'

                }
                else {
                    buffer << line << '\n'
                }

            }

            if ( buffer.size() ) {
                block.call(buffer.toString())
            }

            // send a poison pill to any written channel
            def close = options.autoClose
            if( close == null || close == true ) {
                interceptor.closeChannels()
            }
            else {
                log.trace "Skipping channel autoClose"
            }
        }
        finally {
            reader0.closeQuietly()
        }
    }



    // --==== Dataflow each operator extension ===--


    /**
     * Implements 'each' iterator for {@code DataflowQueue}.
     * <p>
     *     The key feature of this operator is to spread transparently the 'poison pill'
     *     from the source {@code queue} into the referenced queue in the code block
     * <p>
     *     For example:
     *     <code>
     *     source = new DataflowQueue()
     *     target1 = new DataflowQueue()
     *     target2 = new DataflowQueue()
     *
     *     source.each {
     *          target1 << it
     *          target2 << it * it
     *     }
     *     </code>
     *
     *
     *
     * @param queue
     * @param group
     * @param code
     */
    @Deprecated
    static void each( DataflowStreamReadAdapter queue, Closure closure ) {
        each(queue, Dataflow.retrieveCurrentDFPGroup(), closure)
    }

    @Deprecated
    static void each( DataflowStreamReadAdapter queue, final Pool pool, Closure closure ) {
        each( queue, new DefaultPGroup(pool), closure)
    }

    @Deprecated
    static void each( DataflowStreamReadAdapter queue, final PGroup group, Closure code ) {
        each0( queue, group, code )
    }

    @Deprecated
    static void each( DataflowQueue queue, Closure closure ) {
        each(queue, Dataflow.retrieveCurrentDFPGroup(), closure)
    }

    @Deprecated
    static void each( DataflowQueue queue, final Pool pool, Closure closure ) {
        each( queue, new DefaultPGroup(pool), closure)
    }

    @Deprecated
    static void each( DataflowQueue queue, final PGroup group, Closure code ) {
        each0( queue, group, code )
    }

    @Deprecated
    static void each ( WriteChannelWrap channel, Closure code ) {
        each0( channel.target as DataflowReadChannel, Dataflow.retrieveCurrentDFPGroup(), code )
    }

    /**
     * Implements the 'each' operator.
     * <p>
     * The goal is intercept any writes to 'dataflow' queues so that to
     * being able to spread the 'poison pill' from the source 'queue' to the target
     *
     * @param queue A a {@code DataflowQueue} or {@code DataflowBroadcast} object instance
     * @param group
     * @param code The closure code block
     */
    private static void each0( DataflowReadChannel queue, final PGroup group, Closure code ) {
        log.warn "The operator 'each' has been deprecated -- use 'subscribe' instead"

        // -- wrap the owner to intercept any reference to an external dataflow instance
        final interceptor = new WritableChannelInterceptor(code)

        // -- when a 'PoisonPill' is received, spread it over over any written channel
        def listenForPoisonPill = new DataflowEventAdapter() {
            public Object controlMessageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
                if ( message instanceof PoisonPill ) {
                    interceptor.getWrittenChannels() *. bind( message )
                }
                return message;
            }

            public boolean onException(final DataflowProcessor processor, final Throwable e) {
                NextflowExtensions.log.error("${e.getMessage()} -- See the file '.nextflow.log' for more error details", e)
                return true
            }
        }

        group.operator(inputs:[queue], outputs:[], listeners: [listenForPoisonPill], code )
    }

    /**
     * Implements a semantically equivalent 'each' iterator over a {@code DataflowStreamReadAdapter} channel

     * @param channel
     * @param closure
     */
    @Deprecated
    static void eachWithIndex( DataflowStreamReadAdapter channel, Closure closure ) {
        eachWithIndex(channel, Dataflow.retrieveCurrentDFPGroup(), closure)
    }

    @Deprecated
    static void eachWithIndex( DataflowStreamReadAdapter channel, final Pool pool, Closure closure ) {
        eachWithIndex( channel, new DefaultPGroup(pool), closure)
    }

    @Deprecated
    static void eachWithIndex( DataflowStreamReadAdapter channel, final PGroup group, Closure code ) {
        eachWithIndex0( channel, group, code )
    }

    @Deprecated
    static void eachWithIndex( DataflowQueue queue, Closure closure ) {
        eachWithIndex(queue, Dataflow.retrieveCurrentDFPGroup(), closure)
    }

    @Deprecated
    static void eachWithIndex( DataflowQueue channel, final Pool pool, Closure closure ) {
        eachWithIndex( channel, new DefaultPGroup(pool), closure)
    }

    @Deprecated
    static void eachWithIndex( DataflowQueue channel, final PGroup group, Closure code ) {
        eachWithIndex0( channel, group, code )
    }

    @Deprecated
    static void eachWithIndex( WriteChannelWrap channel, Closure code ) {
        eachWithIndex0( channel.target as DataflowReadChannel, Dataflow.retrieveCurrentDFPGroup(), code )
    }

    @Deprecated
    private static void eachWithIndex0( DataflowReadChannel channel, final PGroup group, Closure code ) {
        log.warn "The operator 'eachWithIndex' has been deprecated -- use 'subscribe' instead"

        // -- wrap the owner to intercept any reference to an external dataflow instance
        final interceptor = new WritableChannelInterceptor(code)

        // -- when a 'PoisonPill' is received, spread it over over any written channel
        def listenForPoisonPill = new DataflowEventAdapter() {
            public Object controlMessageArrived(final DataflowProcessor arg0, final DataflowReadChannel<Object> arg1, final int arg2, final Object message) {
                if ( message instanceof PoisonPill ) {
                    interceptor.getWrittenChannels() *. bind(message)
                }
                return message;
            }

            public boolean onException(final DataflowProcessor processor, final Throwable e) {
                NextflowExtensions.log.error("${e.getMessage()} -- See the file '.nextflow.log' for more error details", e)
                return true
            }
        }

        def index = new AtomicLong()
        group.operator(inputs:[channel], outputs:[], listeners: [listenForPoisonPill]) { entry ->
             // invoke the user 'code' passing the current index as an extra parameter
            code.call(entry,index.getAndIncrement())

        }
    }


    /**
     * Keep trace of any reference to a {@code DataflowQueue} or {@code DataflowBroadcast}
     */
    private static class WritableChannelInterceptor {

        /* The target owner object */
        def target

        /* Any reference to a {@code DataflowQueue} or {@code DataflowBroadcast} in the closure block */
        def List<WriteChannelWrap> channels = []

        def private added = []

        def private boolean closed = false

        WritableChannelInterceptor( Closure code ) {
            assert code
            // replace the closure 'owner' by 'this' instance
            target = code.owner
            code.@owner = this
        }

        /** All the channel for which a 'bind' operation has been invoked */
        List<WriteChannelWrap> getWrittenChannels() {
            channels.findAll{ WriteChannelWrap it -> it.receivedData }
        }

        def boolean isClosed() { closed }

        def void closeChannels( ) {
            if( !closed ) {
                channels *.unwrap() .each { channel -> channel << PoisonPill.instance }
                closed = true
            }
        }

        def getProperty(String name) {

            def result = target.getProperty(name)
            if( result instanceof DataflowQueue && !added.contains(result)) {
                added << result
                channels << ( result = new WriteChannelWrap(result) )
            }
            else if ( result instanceof DataflowBroadcast && !added.contains(result) ) {
                added << result
                channels << ( result = new WriteChannelWrap(result))
            }

            return result
        }

    }

    /**
     * Wrap a {@code WriteChannelWrap} a keep track of any bind operation
     */
    @TupleConstructor
    private static class WriteChannelWrap {

        @Delegate
        DataflowWriteChannel target

        /** Whenever a 'bind' operation has been invoked on the target channel */
        boolean receivedData

        /** The reference to the target channel */
        DataflowWriteChannel unwrap() { target }

        DataflowWriteChannel leftShift(final Object value) {
            receivedData = true
            target.leftShift(value)
        }

        void bind(final Object value) {
            receivedData = true
            target.bind(value)
        }

        DataflowWriteChannel leftShift(final DataflowReadChannel ref) {
            receivedData = true
            target.leftShift(ref)
        }

    }


    /**
     * INTERNAL ONLY API
     * <p>
     * Add the {@code update} method to an {@code Agent} so that it call implicitly
     * the {@code Agent#updateValue} method
     *
     */

    static void update( Agent self, Closure message ) {
        assert message != null

        self.send {
            message.call(it)
            updateValue(it)
        }

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


    private static Lock MAP_LOCK = new ReentrantLock()

    static def <T> T getOrCreate( Map self, key, factory ) {

        MAP_LOCK.withLock {
            if( self.containsKey(key) )
                return self.get(key)

            def result = factory instanceof Closure ? factory.call(key) : factory
            self.put(key,result)
            return result
        }

    }


}
