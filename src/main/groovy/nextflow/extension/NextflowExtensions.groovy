/*
 * Copyright (c) 2012, the authors.
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
import nextflow.util.CacheHelper
import nextflow.util.Duration
import org.codehaus.groovy.runtime.DefaultGroovyMethods
import org.codehaus.groovy.runtime.StringGroovyMethods

/**
 * Provides extension methods to chunk text and file
 *
 * See more about extension methods
 * http://docs.codehaus.org/display/GROOVY/Creating+an+extension+module
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class NextflowExtensions {

    static private PATTERN_RIGHT_TRIM = ~/\s+$/

    static private PATTERN_LEFT_TRIM = /^\s+/

    static private PATTERN_FASTA_DESC = ~/^\S+\s+(.*)/

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


    @groovy.transform.PackageScope
    static chopInvoke( Closure closure, Object obj, int index ) {
        if( !closure ) return obj
        def types = closure.getParameterTypes()
        if( types.size()>1 ) {
            return closure.call(obj, index)
        }
        else {
            return closure.call(obj)
        }
    }

    /**
     * Parse a {@code CharSequence} as a FASTA formatted text, retuning a {@code Map} object
     * containing the fields as specified by the @{code record} parameter.
     * <p>
     *  For example:
     *  <pre>
     *  def fasta = '''
     *      >5524211 cytochrome b [Elephas maximus maximus]
     *      LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV
     *      IENY/
     *      '''.stripIndent()
     *
     *   def record = fasta.parseFastaRecord( [ id: true, seq: true ]
     *   assert record.id == '5524211'
     *   assert record.seq = 'LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVIENY'
     *  </pre>
     *
     *
     * @param fasta
     *      The fasta formatted text to be parsed
     * @param record
     *      The map object that is used to specify which fields are required to be returned in the result map.
     *      The following field can be used:
     *      <li>{@code id} The fasta ID
     *      <li>{@code seq} The sequence string
     *      <li>{@code desc} The description in the fasta header
     *      <li>{@code head} The fasta header (first line w/o the '>' character)
     *      <li>{@code text} The complete fasta text block
     *      <li>{@code width} The width of the fasta formatted block. If 0 is specified the sequence is not broken into multiple lines
     *      <li>{@code hash} The hashCode of the entered FASTA sequence
     *      <li>{@code uuid} A random {@code UUID} id for this sequence
     *
     *
     * @return
     */
    static parseFastaRecord( CharSequence fasta, Map record ) {
        assert fasta != null
        if( !record ) return

        String head = null

        def body = new StringBuilder()
        fasta.eachLine { String line ->
            line = line.trim()
            if( !line ) return
            if( line.startsWith(';')) return
            if( !head ) {
                head = line.substring(1);
                return
            }

            if( body.size() )
                body.append('\n')

            body.append(line)
        }

        def result = [:]
        if( record.id && head ) {
            int p = head.indexOf(' ')
            result.id = p != -1 ? head.substring(0,p) : head
        }
        if( record.desc && head ) {
            def m = PATTERN_FASTA_DESC.matcher(head)
            result.desc = m.matches() ? m.group(1) : null
        }
        if( record.head ) {
            result.head = head
        }
        if( record.text ) {
            result.text = fasta.toString()
        }
        if( record.seq ) {
            if( record.width == null ) {
                result.seq = body.toString()
            }
            else if( record.width.toString().isInteger()) {
                def buff = new StringBuilder()
                int len = record.width.toInteger()

                if( len > 0 ) {
                    body.chopString(count:len, ignoreNewLine: true) { str ->
                        if( buff.size() ) buff.append('\n')
                        buff.append(str)
                    }
                }
                else {
                    body.eachLine { buff.append(it) }
                }

                result.seq = buff.toString()
            }
            else {
                throw new IllegalArgumentException("Invalid 'width' argument value: ${record.width}")
            }
        }
        if( record.uuid ) {
            result.uuid = UUID.randomUUID()
        }
        if( record.hash ) {
            result.hash = CacheHelper.hasher(fasta).hash()
        }

        return result
    }

    @groovy.transform.PackageScope
    static Charset getCharset( object ) {

        if( object instanceof Map ) {
            if( object.containsKey('charset') )
                object = object.charset
            else
                return Charset.defaultCharset()
        }

        if( object instanceof Charset )
            return object

        if( object instanceof String && Charset.isSupported(object) )
            return Charset.forName(object)

        if( object != null )
            log.warn "Invalid charset object: $object -- using defualt: ${Charset.defaultCharset()}"

        Charset.defaultCharset()
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
    static chopLines( Object obj, Map options = [:], Closure closure ) {

        if( obj instanceof CharSequence )
            return chopLinesImpl( new StringReader(obj.toString()), options, closure )

        if( obj instanceof Reader )
            return chopLinesImpl( obj, options, closure )

        if( obj instanceof Path ) {
            def charset = getCharset(options)
            return chopLinesImpl( Files.newBufferedReader(obj, charset), options, closure )
        }

        if( obj instanceof InputStream ) {
            def charset = getCharset(options)
            return chopLinesImpl( new InputStreamReader(obj,charset), options, closure )
        }

        if( obj instanceof File ) {
            return chopLinesImpl( new FileReader(obj), options, closure )
        }

        throw new IllegalAccessException("Object of class '${obj.class.name}' does not support 'chopText' method")
    }

    static private chopLinesImpl( Reader reader, Map options = null , Closure closure ) {
        assert reader != null
        assert options != null

        log.trace "Chop options: ${options}"
        final size = options.count ?: 1
        final into = options.into
        if( into && !(into instanceof Collection) && !(into instanceof DataflowQueue) )
            throw new IllegalArgumentException("Argument 'into' can be a subclass of Collection or a DataflowQueue type -- Entered value type: ${into.class.name}")

        BufferedReader reader0 = reader instanceof BufferedReader ? reader : new BufferedReader(reader)

        def result = null
        String line
        int index = 0
        StringBuilder buffer = new StringBuilder()
        int c=0

        try {

            while( (line = reader0.readLine()) != null ) {
                if ( c ) buffer << '\n'
                buffer << line
                if ( ++c == size ) {
                    c = 0
                    result = chopInvoke(closure, buffer.toString(), index++ )
                    if( into != null )
                        into << result

                    buffer.setLength(0)
                }
            }

        }
        finally {
            reader0.closeQuietly()
        }

        /*
         * if there's something remaining in the buffer it's supposed
         * to be the last entry
         */
        if ( buffer.size() ) {
            result = chopInvoke(closure, buffer.toString(), index )
            if( into != null )
                into << result
        }

        /*
         * now close and return the result
         * - when the target it's a channel, send stop message
         * - when it's a list return it
         * - otherwise return the last value
         */
        if( into instanceof DataflowWriteChannel ) {
            into << PoisonPill.instance
            return into
        }
        if( into != null )
            return into

        return result
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
    static chopString( Object obj, Map options, Closure closure ) {

        if( obj instanceof CharSequence )
            return chopStringImpl( new StringReader(obj.toString()), options, closure )

        if( obj instanceof Path ) {
            def charset = getCharset(options.charset)
            return chopStringImpl( Files.newBufferedReader(obj, charset), options, closure )
        }

        if( obj instanceof Reader )
            return chopStringImpl( obj, options, closure )

        if( obj instanceof InputStream ) {
            def charset = getCharset(options)
            return chopStringImpl( new InputStreamReader(obj,charset), options, closure )
        }

        if( obj instanceof File )
            return chopStringImpl( new FileReader(obj), options, closure )

        throw new IllegalAccessException("Object of class '${obj.class.name}' does not support 'chopString' method")
    }

    static private chopStringImpl( Reader reader, Map options = [:] , Closure closure  ) {
        assert reader != null
        assert options != null

        log.trace "Chop options: ${options}"
        final count = options.count ?: 1
        final into = options.into
        final ignoreNewLine = options.ignoreNewLine == true ?: false
        if( into && !(into instanceof Collection) && !(into instanceof DataflowQueue) )
            throw new IllegalArgumentException("Argument 'into' can be a subclass of Collection or a DataflowQueue type -- Entered value type: ${into.class.name}")

        def result = null
        def index = 0
        def buffer = new StringBuilder()
        int c = 0
        def ch

        try {

            while( (ch=reader.read()) != -1 ) {
                if( ignoreNewLine && ( ch == '\n' as char || ch == '\r' as char ))
                    continue
                buffer.append( (char)ch )
                if ( ++c == count ) {
                    c = 0
                    result = chopInvoke(closure, buffer.toString(), index++ )
                    if( into != null )
                        into << result

                    buffer.setLength(0)
                }
            }

        }
        finally {
            reader.closeQuietly()
        }

        /*
         * if there's something remaining in the buffer it's supposed
         * to be the last entry
         */
        if ( buffer.size() ) {
            result = chopInvoke(closure, buffer.toString(), index++ )
            if( into != null )
                into << result
        }

        /*
         * now close and return the result
         * - when the target it's a channel, send stop message
         * - when it's a list return it
         * - otherwise return the last value
         */
        if( into instanceof DataflowWriteChannel ) {
            into << PoisonPill.instance
            return into
        }
        if( into != null )
            return into

        return result
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
    static chopBytes( Object obj, Map options, Closure closure ) {

        if( obj instanceof CharSequence )
            return chopBytesImpl( new ByteArrayInputStream(obj.toString().bytes), options, closure )

        if( obj instanceof Path ) {
            return chopBytesImpl( Files.newInputStream(obj), options, closure )
        }

        if( obj instanceof InputStream ) {
            return chopBytesImpl( obj, options, closure )
        }

        if( obj instanceof File ) {
            return chopBytesImpl( new FileInputStream(obj), options, closure )
        }

        throw new IllegalAccessException("Object of class '${obj.class.name}' does not support 'chopBytes' method")

    }

    static private chopBytesImpl( InputStream stream, Map options = [:], Closure<byte[]> closure ) {
        assert stream != null
        assert options != null

        log.trace "Chop options: ${options}"
        final count = options.count ?: 1
        final into = options.into
        if( into && !(into instanceof Collection) && !(into instanceof DataflowQueue) )
            throw new IllegalArgumentException("Argument 'into' can be a subclass of Collection or a DataflowQueue type -- Entered value type: ${into.class.name}")

        def result = null

        int c=0
        int index = 0
        byte[] buffer = new byte[count]
        byte item

        try {

            while( (item=stream.read()) != -1 ) {
                buffer[c] = (byte)item

                if ( ++c == count ) {
                    c = 0
                    result = chopInvoke(closure, buffer, index++ )
                    if( into != null )
                        into << result

                    buffer = new byte[count]
                }
            }

        }
        finally {
            stream.closeQuietly()
        }


        /*
         * if there's something remaining in the buffer it's supposed
         * to be the last entry
         */

        if ( c ) {
            if( c != count ) {
                def copy = new byte[c]
                System.arraycopy(buffer,0,copy,0,c)
                buffer = copy
            }

            result = chopInvoke(closure, buffer, index++ )
            if( into != null )
                into << result
        }

        /*
         * now close and return the result
         * - when the target it's a channel, send stop message
         * - when it's a list return it
         * - otherwise return the last value
         */
        if( into instanceof DataflowWriteChannel ) {
            into << PoisonPill.instance
            return into
        }
        if( into != null )
            return into

        return result
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
    static chopFasta( Object obj, Map options = [:], Closure closure ) {

        if( obj instanceof CharSequence )
            return chopFastaImpl( new StringReader(obj.toString()), options, closure )

        if( obj instanceof Path ) {
            def charset = getCharset(options)
            return chopFastaImpl( Files.newBufferedReader(obj, charset), options, closure )
        }

        if( obj instanceof Reader ) {
            return chopFastaImpl( obj as Reader, options, closure)
        }

        if( obj instanceof File ) {
            return chopFastaImpl( new FileReader(obj), options, closure )
        }

        if( obj instanceof InputStream ) {
            def charset = getCharset(options)
            return chopFastaImpl( new InputStreamReader(obj,charset), options, closure)
        }

        throw new IllegalAccessException("Object of class '${obj.class.name}' does not support 'chopFasta' method")
    }

    static chopFastaImpl( Reader text, Map options = [:], Closure<String> closure ) {
        assert text != null
        assert options != null

        log.trace "Chunk options: $options"
        final count = options.count ?: 1
        final into = options.into
        final rec = options.record && options.record instanceof Map

        if( into && !(into instanceof Collection) && !(into instanceof DataflowQueue) )
            throw new IllegalArgumentException("Argument 'into' can be a subclass of Collection or a DataflowQueue type -- Entered value type: ${into.class.name}")

        if( rec && count>1 )
            throw new IllegalArgumentException("When using 'record' option 'count' cannot be greater than 1")

        BufferedReader reader0 = text instanceof BufferedReader ? text : new BufferedReader(text)

        def result = null
        String line
        StringBuilder buffer = new StringBuilder()
        int index = 0
        int blockCount=0
        boolean openBlock = false

        try {

            while( (line = reader0.readLine()) != null ) {

                if( line.startsWith(';')) continue

                if ( line == '' ) {
                    buffer << '\n'
                }
                else if ( !openBlock && line.charAt(0)=='>' ) {
                    openBlock = true
                    buffer << line << '\n'
                }
                else if ( openBlock && line.charAt(0)=='>') {
                    // another block is started

                    if ( ++blockCount == count ) {
                        // invoke the closure, passing the read block as parameter
                        def record = rec ? parseFastaRecord(buffer.toString(), (Map)options.record) : buffer.toString()
                        result = chopInvoke(closure, record, index++ )
                        if( into != null ) {
                            into << result
                        }

                        buffer.setLength(0)
                        blockCount=0
                    }

                    buffer << line << '\n'

                }
                else {
                    buffer << line << '\n'
                }
            }

        }
        finally {
            reader0.closeQuietly()
        }


        /*
         * if there's something remaining in the buffer it's supposed
         * to be the last entry
         */
        if ( buffer.size() ) {
            def record = rec ? parseFastaRecord(buffer.toString(), (Map)options.record) : buffer.toString()
            result = chopInvoke(closure, record, index )
            if( into != null )
                into << result
        }

        /*
         * now close and return the result
         * - when the target it's a channel, send stop message
         * - when it's a list return it
         * - otherwise return the last value
         */
        if( into instanceof DataflowWriteChannel ) {
            into << PoisonPill.instance
            return into
        }
        if( into != null )
            return into

        return result
    }


    /**
     * Splits a {@code CharSequence} in text chunks having the specified number of lines
     *
     * @param sequence A sequence of chars to which apply the chunkLines operation
     * @param n The number of lines in each 'chunk'
     * @param A closure invoked for each chunk of text, it takes a single (optional) string parameter which represents the current text chunk
     */
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
    static void each( DataflowStreamReadAdapter queue, Closure closure ) {
        each(queue, Dataflow.retrieveCurrentDFPGroup(), closure)
    }

    static void each( DataflowStreamReadAdapter queue, final Pool pool, Closure closure ) {
        each( queue, new DefaultPGroup(pool), closure)
    }

    static void each( DataflowStreamReadAdapter queue, final PGroup group, Closure code ) {
        each0( queue, group, code )
    }


    static void each( DataflowQueue queue, Closure closure ) {
        each(queue, Dataflow.retrieveCurrentDFPGroup(), closure)
    }

    static void each( DataflowQueue queue, final Pool pool, Closure closure ) {
        each( queue, new DefaultPGroup(pool), closure)
    }

    static void each( DataflowQueue queue, final PGroup group, Closure code ) {
        each0( queue, group, code )
    }


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
    static void eachWithIndex( DataflowStreamReadAdapter channel, Closure closure ) {
        eachWithIndex(channel, Dataflow.retrieveCurrentDFPGroup(), closure)
    }

    static void eachWithIndex( DataflowStreamReadAdapter channel, final Pool pool, Closure closure ) {
        eachWithIndex( channel, new DefaultPGroup(pool), closure)
    }

    static void eachWithIndex( DataflowStreamReadAdapter channel, final PGroup group, Closure code ) {
        eachWithIndex0( channel, group, code )
    }


    static void eachWithIndex( DataflowQueue queue, Closure closure ) {
        eachWithIndex(queue, Dataflow.retrieveCurrentDFPGroup(), closure)
    }

    static void eachWithIndex( DataflowQueue channel, final Pool pool, Closure closure ) {
        eachWithIndex( channel, new DefaultPGroup(pool), closure)
    }

    static void eachWithIndex( DataflowQueue channel, final PGroup group, Closure code ) {
        eachWithIndex0( channel, group, code )
    }


    static void eachWithIndex( WriteChannelWrap channel, Closure code ) {
        eachWithIndex0( channel.target as DataflowReadChannel, Dataflow.retrieveCurrentDFPGroup(), code )
    }


    private static void eachWithIndex0( DataflowReadChannel channel, final PGroup group, Closure code ) {

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


}
