/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
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

package nextflow.splitter
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import nextflow.util.CacheHelper
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

    @Override
    protected Map<String,Object> validOptions() {
        def result = super.validOptions()
        result.record = [ Boolean, Map ]
        return result
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
        fasta.eachLine { String line ->
            line = line.trim()
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
            body.eachLine { buff.append(it) }
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


    private StringBuilder buffer = new StringBuilder()

    private String line

    @Override
    protected process( Reader targetObject ) {
        line = null
        super.process(targetObject)
    }

    @Override
    protected fetchRecord(BufferedReader reader) {

        buffer.setLength(0)
        boolean openBlock = false

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
            }

            // next line
            line = reader.readLine()
        }

        String result = buffer.size() ? buffer.toString() : null
        if( result && recordMode )
            return parseFastaRecord(result, recordFields)

        return result
    }
}
