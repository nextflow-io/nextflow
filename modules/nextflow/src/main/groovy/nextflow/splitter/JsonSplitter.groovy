/*
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

import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;

/**
 * Split a Json file in records
 *
 * @author Pierre Lindenbaum PhD Institut-du-Thorax Nantes France.
 */
@Slf4j
@CompileStatic
@InheritConstructors
class JsonSplitter extends AbstractTextSplitter {
    public static final String OBJECT_KEY = "key";
    public static final String OBJECT_VALUE = "value";

    /**
     * json path to find a specific region of the json
     */
    private String jsonPath;

    /** number of items processed so far */
    private long itemsCount = 0L;

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
        result.remove('compress')
        result.remove('keepHeader')
        result.remove('by')
        result.remove('elem')
        result.jsonPath= [String]
        return result
    }

    @Override
    protected process( Reader  reader) {
        def result = null
        counter.reset() // <-- make sure to start
        itemsCount = 0L;
        try {
            final JsonReader jsonreader = new JsonReader(reader);
            result = scanJsonForPath(jsonreader,
                    "",
                    (this.jsonPath==null?"":this.jsonPath)
                    );
            jsonreader.close();
            if (result!=null) {
                // make sure to process collected entries
                if ( collector && collector.hasChunk() ) {
                    result = invokeEachClosure(closure, collector.nextChunk())
                }
            }
        }
        catch (Exception ex) {
            log.warn "Something was wrong when parsing JSON. reason: ${ex.message ?: ex}", ex
            throw ex;
        }
        finally {
            reader.closeQuietly()
            if( collector instanceof Closeable )
                collector.closeQuietly()
        }
        return result
    }

    /** recursive loop over the json tree until  this.jsonPath was found
     * @param reader the JsonReader
     * @param currentPath json path. Root is an empty string
     * @param expectedPath expected Path
     * @returns the output of processArray/processObject or null
     */
    private scanJsonForPath(
            JsonReader reader,
            String currentPath,
            String expectedPath
    ) {
        // did we reach the expected path ?
        boolean got_path =  currentPath.equals(expectedPath)

        JsonToken token = reader.peek();
        if( token.equals(JsonToken.END_DOCUMENT)) {
            return null;
        }
        else if (token.equals(JsonToken.BEGIN_ARRAY)) {
            if(got_path) {
                return processArray(reader);
            }
            else /* loop over the current array */
            {
                int item_index = 0;
                reader.beginArray();
                while(reader.hasNext()) {
                    token = reader.peek();
                    /* end of array */
                    if ( token.equals(JsonToken.END_ARRAY) ) {
                        break;
                    }
                    def ret = scanJsonForPath(reader, currentPath + "[" + item_index + "]", expectedPath );
                    if (ret!=null) return ret;
                    item_index++;
                }
                reader.endArray();
            }
        }
        else if(token.equals(JsonToken.BEGIN_OBJECT)) {
            if(got_path) {
                return processObject(reader);
            }
            else /* loop over the current object */
            {
                reader.beginObject();
                while(reader.hasNext()) {
                    token = reader.peek();
                    /* end of array */
                    if ( token.equals(JsonToken.END_OBJECT) ) {
                        break;
                    }
                    final String key = reader.nextName();
                    def ret = scanJsonForPath(reader, currentPath + (currentPath.endsWith(".") || currentPath.isEmpty()? "": ".") + key, expectedPath);
                    if (ret!=null) return ret;
                }
                reader.endObject();
            }
        }
        else {
            if(got_path)  throw new IOException(
                "Expected JSON to be an object or an array but got " +
                token.name() + " for json-path \"" + currentPath +"\"."
                );
            reader.skipValue();
        }
        return null;
    }


    /** process and emit each record of a JSON array */
    private processArray( JsonReader reader ) {
        def result = null
        reader.beginArray();
        while(reader.hasNext()) {
            final JsonToken token = reader.peek();
            /* end of array */
            if ( token.equals(JsonToken.END_ARRAY) ) {
                reader.endArray();
                break;
            }
            final Object value = fromJson( reader );
            result = processChunk( value )

            // -- check the limit of allowed records has been reached
            if( limit>0 && ++itemsCount == limit ) {
                break
            }
        }

        return result;
    }

    /** process and emit each record of a JSON object */
    private processObject( JsonReader reader ) {
        def result = null
        reader.beginObject();
        while(reader.hasNext()) {
            final JsonToken token = reader.peek();
            /* end of array */
            if ( token.equals(JsonToken.END_OBJECT) ) {
                reader.endObject();
                break;
            }
            final Map<String,Object> map = new HashMap<>(2);
            final String key = reader.nextName();
            final Object value = fromJson(reader);
            map.put(OBJECT_KEY, key);
            map.put(OBJECT_VALUE, value);
            // -- apply the splitting logic for the fetched record
            result = processChunk( map )

            // -- check the limit of allowed records has been reached
            if( limit>0 && ++itemsCount == limit ) {
                break
            }
        }
        return result;
    }


    @Override
    protected fetchRecord(BufferedReader reader) {
        throw new IllegalStateException("should never be called");
    }

    /** convert a json stream to an Object */
    static Object fromJson(final JsonReader reader ) {
        final JsonToken token = reader.peek();
        switch(token) {
            case JsonToken.NULL : {
                reader.nextNull();
                return null;
            }
            case JsonToken.STRING : {
                return reader.nextString();
            }
            case JsonToken.BOOLEAN : {
                return reader.nextBoolean();
            }
            case JsonToken.NUMBER : {
                final String s = reader.nextString();
                // look like an integer ?
                if(!s.contains(".")) {
                    // first try as int
                    try {
                        int v = Integer.parseInt(s);
                        return v;
                    }
                    catch (NumberFormatException ex) {
                        //otherwise try as long
                        try {
                            long v = Long.parseLong(s);
                            return v;
                        }
                        catch (NumberFormatException ex2) {
                            //ignore
                        }
                    }
                }
                // try as double
                try {
                    double v = Double.parseDouble(s);
                    if(Double.isInfinite(v) || Double.isNaN(v)) return s;
                    return v;
                }
                catch (NumberFormatException ex) {
                    //ignore
                }
                //everything failed, return as string
                return s;
            }
            case JsonToken.BEGIN_ARRAY : {
                reader.beginArray();
                final List<Object> array = new ArrayList<>();
                for(;;) {
                    final JsonToken token2 = reader.peek();
                    if(token2.equals(JsonToken.END_ARRAY)) {
                        reader.endArray();
                        break;
                    }
                    array.add(fromJson(reader));
                }
                return array;
            }
            case JsonToken.BEGIN_OBJECT : {
                reader.beginObject();
                final Map<String,Object> object = new LinkedHashMap<>();
                for(;;) {
                    final JsonToken token2 = reader.peek();
                    if(token2.equals(JsonToken.END_OBJECT)) {
                        reader.endObject();
                        break;
                    }
                    final String key = reader.nextName();
                    final Object value = fromJson(reader);
                    object.put(key,value);
                }
                return object;
            }
            default: throw new IllegalStateException("Got json event: " + token);
        }
    }
}
