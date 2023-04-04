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
 * @author Pierre Lindenbaum PhD
 */
@Slf4j
@CompileStatic
@InheritConstructors
class JsonSplitter extends AbstractTextSplitter {
    /**
     * use lenient json parsing for JsonReader
     */
    private boolean lenient  = false;

    /**
     * Set the splitter options by specifying a map of named parameters.
     * Valid parameters are:
     * <li>{@code lenient}
     *
     * @param options
     * @return The splitter instance itself
     */
    @Override
    JsonSplitter options(Map opts) {
        super.options(opts)
        this.lenient = opts.lenient == true		
        return this
    }

    /**
     * @return A map representing the valid options for the splitter. The map keys define the
     * accepted parameter names, the values the valid values for each of them.
     */
    @Override
    protected Map<String,?> validOptions() {
        def result = super.validOptions()
	result.remove('file')       // <-- `file` mode not support by JsonSplitter (not sure about this)
        result.remove('compress')   // <-- `compress` mode not supported
        result.put("lenient",false);
        return result
    }

    @Override
    protected process( Reader  reader) {
        def result = null
        counter.reset() // <-- make sure to start
		itemsCount = 0
        try {
            final JsonReader jsonreader = new JsonReader(reader);
            jsonreader.setLenient(this.lenient);
        
            final JsonToken token = jsonreader.peek();
            if (token.equals(JsonToken.BEGIN_ARRAY)) {
            	processArray(jsonreader);
            	}
            else if(token.equals(JsonToken.BEGIN_OBJECT)) {
            	processObject(jsonreader);
            	}
            else {
            	throw new IOException("Expected JSON to be an object or an array but got " + token);
            	}
			/* extra content at the end ? */
			if(jsonreader.hasNext()) {
				final JsonToken last =  jsonreader.peek();
				if(!token.equals(JsonToken.END_DOCUMENT)) {
						log.warn("Extra content at the end of the json stream");
						}
				}
            jsonreader.close();
            }
        finally {
            reader.closeQuietly()
            if( collector instanceof Closeable )
                collector.closeQuietly()
        }

        return result
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
			 if( limit>0 && ++itemsCount == limit )
				 break
			}
	
	 // make sure to process collected entries
	    if ( collector && collector.hasChunk() ) {
		result = invokeEachClosure(closure, collector.nextChunk())
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
    	 final Map<String,Object> map = new HashMap<>();
    	 final String key = reader.nextName();
    	 final Object value = fromJson(reader);
    	 map.put(key,value);
    	  // -- apply the splitting logic for the fetched record
         result = processChunk( map )
		 
		 // -- check the limit of allowed records has been reached
		 if( limit>0 && ++itemsCount == limit )
			 break
    	 }

	 // make sure to process collected entries
	    if ( collector && collector.hasChunk() ) {
		result = invokeEachClosure(closure, collector.nextChunk())
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
