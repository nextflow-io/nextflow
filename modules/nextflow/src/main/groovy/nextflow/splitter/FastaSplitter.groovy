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

package nextflow.splitter
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.util.CacheHelper
import nextflow.util.MemoryUnit
/**
 * Split FASTA formatted text content or files
 *
 * @link http://en.wikipedia.org/wiki/FASTA_format
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@InheritConstructors
class FastaSplitter extends AbstractTextSplitter {

    static private Pattern PATTERN_FASTA_DESC = ~/^\S+\s+(.*)/

    private boolean incrementBySize

    private StringBuilder buffer = new StringBuilder()

    private String line

    @Override
    protected Map<String,Object> validOptions() {
        def result = super.validOptions()
        result.record = [ Boolean, Map ]
        result.size = [String, MemoryUnit, Number]
        return result
    }

    FastaSplitter options( Map options ) {
        super.options(options)

        if( options.size && options.by )
            throw new AbortOperationException("Parameter `by` and `size` conflicts -- check operator `$operatorName`")

        if( options.size ) {
            final size = parseChunkSize(options.size)
            counter = new EntryCounter(size, true)
            incrementBySize = true
        }

        return this
    }


    protected long parseChunkSize( value ) {
        assert value != null
        if( value instanceof Number )
            return value.toLong()
        if( value instanceof CharSequence )
            return MemoryUnit.of(value.toString()).toBytes()
        if( value instanceof MemoryUnit )
            return value.toBytes()
        throw new IllegalArgumentException("Not a valid `size` value: $value [${value.getClass().getName()}] -- check operator `${getOperatorName()}`")
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
     *   assert record.sequence = 'LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVIENY'
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
     *      <li>{@code header} The fasta header (first line including the '>' character)
     *      <li>{@code text} The complete fasta text block
     *      <li>{@code width} The width of the fasta formatted block.
     *      <li>{@code string} The sequence is returned as single line string (w/o newline char)
     *      <li>{@code hash} The hashCode of the entered FASTA sequence
     *      <li>{@code uuid} A random {@code UUID} id for this sequence
     *
     *
     * @return
     */
    final static parseFastaRecord( CharSequence fasta, Map record ) {
        assert fasta != null
        if( !record ) return

        String header = null

        def body = new StringBuilder()
        fasta.eachLine { line ->
            line = (line as String).trim()
            if( !line ) return
            if( line.startsWith(';')) return
            if( !header ) {
                if( line.startsWith('>') )
                    header = line.substring(1)
                return
            }

            body.append(line).append('\n')
        }

        def result = [:]
        if( record.id && header ) {
            int p = header.indexOf(' ')
            result.id = p != -1 ? header.substring(0,p) : header
        }

        if( record.desc && header ) {
            def m = PATTERN_FASTA_DESC.matcher(header)
            result.desc = m.matches() ? m.group(1) : null
        }

        if( record.header ) {
            result.put('header', header)
        }

        if( record.text ) {
            result.text = fasta.toString()
        }

        if( record.seqString ) {
            def buff = new StringBuilder()
            body.eachLine { it -> buff.append(it as String) }
            result.seqString = buff.toString()
        }

        if( record.sequence ) {
            if( record.width == null ) {
                result.sequence = body.toString()
            }
            else {
                if( ! record.width.toString().isInteger() )
                    throw new IllegalArgumentException("Invalid 'width' argument value: ${record.width}")

                int len = record.width as int
                if( len <= 0 )
                    throw new IllegalArgumentException("Invalid 'width' argument value: ${record.width} -- it must be a value greater than 0")

                def buff = new StringBuilder()
                new StringSplitter()
                        .options(by:len, ignoreNewLine: true)
                        .target(body)
                        .each { str ->
                            buff.append(str).append('\n')
                        }
                result.sequence = buff.toString()
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


    @Override
    protected process( Reader targetObject ) {
        line = null
        super.process(targetObject)
    }

    @Override
    protected fetchRecord(BufferedReader reader) {

        buffer.setLength(0)
        boolean openBlock = false
        int size = 0

        // check if a line from a previous iteration is available
        if( line==null )
            line = reader.readLine()

        while( line != null ) {
            if( line.startsWith(';') ) {
                // ignore this line
            }
            else if ( line == '' ) {
                buffer << '\n'
            }
            else if ( !openBlock && line.charAt(0)=='>' ) {
                openBlock = true
                buffer << line << '\n'
            }
            else if ( openBlock && line.charAt(0)=='>' ) {   // another block is started
                break
            }
            else {
                buffer << line << '\n'
                size += line.size()
            }

            // next line
            line = reader.readLine()
        }

        if( incrementBySize )
            counter.setIncrement(size)

        String result = buffer.size() ? buffer.toString() : null
        if( result && recordMode )
            return parseFastaRecord(result, recordFields)

        return result
    }
}
