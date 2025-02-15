/*
 * Copyright 2013-2024, Seqera Labs
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
 * @author Pierre Lindenbaum PhD Institut-du-Thorax 44000 Nantes France.
 * Reviewed by bentsherman and pditommaso
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
    private String path

    /**
     * number of items processed so far
     */
    private long itemsCount = 0

    /**
     * Set the splitter options by specifying a map of named parameters.
     * Valid parameters are:
     * <li>{@code path}
     *
     * @param options
     * @return The splitter instance itself
     */
    @Override
    JsonSplitter options(Map opts) {
        super.options(opts)
        if (opts.path) {
            this.path = (opts.path as String)
        }
        return this
    }

    /**
     * @return A map representing the valid options for the splitter. The map keys define the
     * accepted parameter names, the values the valid values for each of them.
     */
    @Override
    protected Map<String,?> validOptions() {
        final baseOptions = super.validOptions()
        return [
            autoClose: baseOptions.autoClose,
            into: baseOptions.into,
            limit: baseOptions.limit,
            path: String,
        ]
    }

    @Override
    protected process(Reader reader) {
        def result = null
        counter.reset() // <-- make sure to start
        itemsCount = 0
        try {
            final jsonReader = new JsonReader(reader)
            result = scanJsonForPath(jsonReader, "", this.path ?: "")

            jsonReader.close()
            if( result && collector?.hasChunk() )
                result = invokeEachClosure(closure, collector.nextChunk())
        }
        catch( Exception ex ) {
            throw new IllegalStateException("Error while parsing JSON content - cause: ${ex.message ?: ex}", ex)
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
            def index = 0
            reader.beginArray()
            while( reader.peek() != JsonToken.END_ARRAY ) {
                final ret = scanJsonForPath(reader, getJsonArraySubPath(currentPath, index), expectedPath)
                if( ret != null )
                    return ret
                index++
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
                final ret = scanJsonForPath(reader, getJsonObjectSubPath(currentPath, key), expectedPath)
                if( ret != null )
                    return ret
            }
            reader.endObject()
        }

        else {
            if( foundPath )
                throw new IllegalStateException("Expected JSON to be an object or an array but got ${token.name()} for json-path \"${currentPath}\".")

            reader.skipValue()
        }

        return null
    }

    /**
     * Get an array index sub-path of a given JSON path.
     *
     * @param path
     * @param index
     */
    private String getJsonArraySubPath(String path, int index) {
        "${path}[${index}]"
    }

    /**
     * Get an object key sub-path of a given JSON path.
     *
     * @param path
     * @param key
     */
    private String getJsonObjectSubPath(String path, String key) {
        "${path}${(path.endsWith('.') || path.isEmpty()) ? '' : '.'}${key}"
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

            final map = new LinkedHashMap<>(1)
            map[OBJECT_KEY] = reader.nextName()
            map[OBJECT_VALUE] = fromJson(reader)

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
