/*
 * Copyright 2013-2023, Seqera Labs
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

import com.google.gson.stream.JsonReader
import com.google.gson.stream.JsonToken
import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j

/**
 * Split a JSON document into records
 *
 * @author Pierre Lindenbaum PhD Institut-du-Thorax Nantes France.
 */
@Slf4j
@CompileStatic
@InheritConstructors
class JsonSplitter extends AbstractTextSplitter {

    public static final String OBJECT_KEY = 'key'
    public static final String OBJECT_VALUE = 'value'

    /**
     * json path to find a specific region of the json
     */
    private String jsonPath

    /**
     * number of items processed so far
     */
    private long itemsCount = 0

    /**
     * Set the splitter options by specifying a map of named parameters.
     * Valid parameters are:
     * <li>{@code jsonPath}
     *
     * @param options
     * @return The splitter instance itself
     */
    @Override
    JsonSplitter options(Map opts) {
        super.options(opts)
        if (opts.jsonPath) {
            this.jsonPath = (opts.jsonPath as String)
        }
        return this
    }

    /**
     * @return A map representing the valid options for the splitter. The map keys define the
     * accepted parameter names, the values the valid values for each of them.
     */
    @Override
    protected Map<String,?> validOptions() {
        def result = super.validOptions()
        result.remove('charset')
        result.remove('compress')
        result.remove('decompress')
        result.remove('keepHeader')
        result.remove('by')
        result.remove('elem')
        result.jsonPath = [String]
        return result
    }

    @Override
    protected process(Reader reader) {
        def result = null
        counter.reset() // <-- make sure to start
        itemsCount = 0
        try {
            final jsonreader = new JsonReader(reader)
            result = scanJsonForPath(jsonreader, "", this.jsonPath ?: "")

            jsonreader.close()
            if( result != null ) {
                // make sure to process collected entries
                if( collector && collector.hasChunk() ) {
                    result = invokeEachClosure(closure, collector.nextChunk())
                }
            }
        }
        catch( Exception ex ) {
            log.warn "Error while parsing JSON: ${ex.message ?: ex}", ex
            throw ex
        }
        finally {
            reader.closeQuietly()
            if( collector instanceof Closeable )
                collector.closeQuietly()
        }
        return result
    }

    /**
     * Return a section of a JSON document specified by a path,
     * or null if not found.
     *
     * @param reader
     * @param currentPath
     * @param expectedPath
     */
    private scanJsonForPath(JsonReader reader, String currentPath, String expectedPath) {
        final foundPath = currentPath == expectedPath
        final token = reader.peek()

        if( token == JsonToken.END_DOCUMENT )
            return null

        else if( token == JsonToken.BEGIN_ARRAY ) {
            if( foundPath )
                return processArray(reader)

            // loop over the current array
            def itemIndex = 0
            reader.beginArray()
            while( reader.peek() != JsonToken.END_ARRAY ) {
                final ret = scanJsonForPath(reader, "${currentPath}[${itemIndex}]", expectedPath)
                if( ret != null )
                    return ret
                itemIndex++
            }
            reader.endArray()
        }

        else if( token == JsonToken.BEGIN_OBJECT ) {
            if( foundPath )
                return processObject(reader)

            // loop over the current object
            reader.beginObject()
            while( reader.peek() != JsonToken.END_OBJECT ) {
                final key = reader.nextName()
                final ret = scanJsonForPath(reader, "${currentPath}${(currentPath.endsWith('.') || currentPath.isEmpty()) ? '' : '.'}${key}", expectedPath)
                if( ret != null )
                    return ret
            }
            reader.endObject()
        }

        else {
            if( foundPath )
                throw new IOException("Expected JSON to be an object or an array but got ${token.name()} for json-path \"${currentPath}\".")

            reader.skipValue()
        }

        return null
    }

    /**
     * Process and emit each record of a JSON array
     *
     * @param reader
     */
    private processArray( JsonReader reader ) {
        def result = null

        reader.beginArray()
        while( reader.hasNext() ) {
            if( reader.peek() == JsonToken.END_ARRAY ) {
                reader.endArray()
                break
            }

            final value = fromJson( reader )
            result = processChunk( value )

            // -- check the limit of allowed records has been reached
            if( limit > 0 && ++itemsCount == limit )
                break
        }

        return result
    }

    /**
     * Process and emit each record of a JSON object
     *
     * @param reader
     */
    private processObject( JsonReader reader ) {
        def result = null

        reader.beginObject()
        while( reader.hasNext() ) {
            if( reader.peek() == JsonToken.END_OBJECT ) {
                reader.endObject()
                break
            }

            final map = [
                (OBJECT_KEY): reader.nextName(),
                (OBJECT_VALUE): fromJson(reader)
            ]

            // -- apply the splitting logic for the fetched record
            result = processChunk( map )

            // -- check the limit of allowed records has been reached
            if( limit > 0 && ++itemsCount == limit )
                break
        }

        return result
    }

    @Override
    protected fetchRecord(BufferedReader reader) {
        throw new IllegalStateException("should never be called")
    }

    /**
     * Convert a json stream to an Object
     *
     * @param reader
     */
    static Object fromJson( JsonReader reader ) {
        final token = reader.peek()
        switch( token ) {
            case JsonToken.NULL:
                reader.nextNull()
                return null

            case JsonToken.STRING:
                return reader.nextString()

            case JsonToken.BOOLEAN:
                return reader.nextBoolean()

            case JsonToken.NUMBER:
                final str = reader.nextString()
                if( str ==~ /\d+(\.\d+)?/ && str.isInteger() ) return str.toInteger()
                if( str ==~ /\d+(\.\d+)?/ && str.isLong() ) return str.toLong()
                if( str ==~ /\d+(\.\d+)?/ && str.isDouble() ) return str.toDouble()
                return str

            case JsonToken.BEGIN_ARRAY:
                final array = []
                reader.beginArray()
                while( reader.peek() != JsonToken.END_ARRAY )
                    array.add(fromJson(reader))
                reader.endArray()
                return array

            case JsonToken.BEGIN_OBJECT:
                final object = [:]
                reader.beginObject()
                while( reader.peek() != JsonToken.END_OBJECT ) {
                    final key = reader.nextName()
                    final value = fromJson(reader)
                    object[key] = value
                }
                reader.endObject()
                return object

            default:
                throw new IllegalStateException("Unexpected JSON token: ${token}")
        }
    }
}
