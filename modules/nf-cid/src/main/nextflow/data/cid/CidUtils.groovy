/*
 * Copyright 2013-2025, Seqera Labs
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
 *
 */

package nextflow.data.cid

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.data.cid.fs.CidPath
import nextflow.data.cid.model.WorkflowRun
import nextflow.data.cid.model.TaskRun
import nextflow.data.cid.serde.CidEncoder
import nextflow.data.cid.serde.CidSerializable
import nextflow.serde.gson.GsonEncoder

import java.nio.file.attribute.FileTime
import java.time.Instant

/**
 * Utils class for CID.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class CidUtils {

    private static final String[] EMPTY_ARRAY = new String[] {}

    /**
     * Query a CID store.
     * @param store CID store to query.
     * @param uri Query to perform in a URI-like format.
     *      Format 'cid://<key>[?QueryString][#fragment]' where:
     *       - Key: Element where the query will be applied. '/' indicates query will be applied in all the elements of the CID store.
     *       - QueryString: all param-value pairs that the CID element should fulfill in a URI's query string format.
     *       - Fragment: Element fragment to retrieve.
     * @return List of object fulfilling the query
     */
    static List query(CidStore store, URI uri) {
        String key = uri.authority ? uri.authority + uri.path : uri.path
        try {
            if (key == CidPath.SEPARATOR) {
                return globalSearch(store, uri)
            } else {
                final parameters = uri.query ? parseQuery(uri.query) : null
                final children = parseChildrenFormFragment(uri.fragment)
                return searchPath(store, key, parameters, children )
            }
        } catch(Throwable e){
            log.debug("Exception querying $uri. $e.message")
            return []
        }

    }

    private static List<CidSerializable> globalSearch(CidStore store, URI uri) {
        final results = store.search(uri.query)
        if (results && uri.fragment) {
            // If fragment is defined get the property of the object indicated by the fragment
            return filterResults(results, uri.fragment)
        }
        return results
    }

    private static List filterResults(List<CidSerializable> results, String fragment) {
        final filteredResults = []
        results.forEach {
            final output = navigate(it, fragment)
            if (output) {
                filteredResults.add(output)
            }
        }
        return filteredResults
    }

    /**
     * Get the array of the search path children elements from the fragment string
     * @param fragment String containing the elements separated by '.'
     * @return array with the parsed element
     */
    static String[] parseChildrenFormFragment(String fragment) {
        if( !fragment )
            return EMPTY_ARRAY
        return fragment.contains('.')
            ?  fragment.split("\\.")
            :  List.of(fragment) as String[]
    }

    /**
     * Search for objects inside a description
     * @param store CID store
     * @param key CID key where to perform the search
     * @param params Parameter-value pairs to be evaluated in the key
     * @param children  Sub-objects to evaluate and retrieve
     * @return List of object
     */
    protected static List<Object> searchPath(CidStore store, String key, Map<String, String> params, String[] children = []) {
        final object = store.load(key)
        if (!object) {
            throw new FileNotFoundException("Cid object $key not found.")
        }
        final results = new LinkedList<Object>()
        if (children && children.size() > 0) {
            treatSubObject(store, key, object, children, params, results)
        } else {
            treatObject(object, params, results)
        }

        return results
    }

    private static void treatSubObject(CidStore store, String key, CidSerializable object, String[] children, Map<String, String> params, LinkedList<Object> results) {
        final output = getSubObject(store, key, object, children)
        if (!output) {
            throw new FileNotFoundException("Cid object $key#${children.join('.')} not found.")
        }
        treatObject(output, params, results)
    }

    /**
     * Get a metadata sub-object.
     * If the requested sub-object is the workflow or task outputs, retrieves the outputs from the outputs description.
     *
     * @param store CidStore to retrieve metadata objects.
     * @param key Parent metadata key.
     * @param object Parent object.
     * @param children Array of string in indicating the properties to navigate to get the sub-object.
     * @return Sub-object or null in it does not exist.
     */
    static Object getSubObject(CidStore store, String key, CidSerializable object, String[] children) {
        if( isSearchingOutputs(object, children) ) {
            // When asking for a Workflow or task output retrieve the outputs description
            final outputs = store.load("${key}/outputs")
            if (!outputs)
                return null
            return navigate(outputs, children.join('.'))
        }
        return navigate(object, children.join('.'))
    }

    /**
     * Check if the Cid pseudo path or query is for Task or Workflow outputs.
     *
     * @param object Parent Cid metadata object
     * @param children Array of string in indicating the properties to navigate to get the sub-object.
     * @return return 'true' if the parent is a Task/Workflow run and the first element in children is 'outputs'. Otherwise 'false'
     */
    static boolean isSearchingOutputs(CidSerializable object, String[] children) {
        return (object instanceof WorkflowRun || object instanceof TaskRun) && children && children[0] == 'outputs'
    }

    /**
     * Evaluates object or the objects in a collection matches a set of parameter-value pairs. It includes in the results collection in case of match.
     * @param object Object or collection of objects to evaluate
     * @param params parameter-value pairs to evaluate in each object
     * @param results results collection to include the matching objects
     */
    protected static void treatObject(def object, Map<String, String> params, List<Object> results) {
        if (params) {
            if (object instanceof Collection) {
                (object as Collection).forEach { treatObject(it, params, results) }
            } else if (checkParams(object, params)) {
                results.add(object)
            }
        } else {
            results.add(object)
        }
    }
    /**
     * Parses a query string and store them in parameter-value Map.
     * @param queryString URI-like query string. (e.g. param1=value1&param2=value2).
     * @return Map containing the parameter-value pairs of the query string.
     */
    static Map<String, String> parseQuery(String queryString) {
        if (queryString) {
            return queryString.split('&').collectEntries {
                it.split('=').collect { URLDecoder.decode(it, 'UTF-8') }
            } as Map<String, String>
        }
        return [:]
    }

    /**
     * Check if an object fullfill the parameter-value
     * @param object Object to evaluate
     * @param params parameter-value pairs to evaluate
     * @return true if all object parameters exist and matches with the value, otherwise false.
     */
    static boolean checkParams(Object object, Map<String,String> params) {
        for (final entry : params.entrySet()) {
            final value = navigate(object, entry.key)
            if (!value || value.toString() != entry.value.toString() ) {
                return false
            }
        }
        return true
    }

    /**
     * Retrieves the sub-object or value indicated by a path.
     * @param obj Object to navigate
     * @param path Elements path separated by '.' e.g. field.subfield
     * @return sub-object / value
     */
    static Object navigate(Object obj, String path){
        if (!obj)
            return null
        // type has been replaced by class when evaluating CidSerializable objects
        if (obj instanceof CidSerializable && path == 'type')
                return obj.getClass()?.simpleName
        try{
            return path.tokenize('.').inject(obj) { current, key ->
                return getSubPath(current, key)
            }
        } catch (Throwable e) {
            log.debug("Error navigating to $path in object", e)
            return null
        }
    }

    private static Object getSubPath(current, String key) {
        if (current == null)
            return null

        if (current instanceof Map) {
            return current[key] // Navigate Map properties
        }
        if (current instanceof Collection) {
            return navigateCollection(current, key)
        }
        if (current.metaClass.hasProperty(current, key)) {
            return current.getAt(key) // Navigate Object properties
        }
        log.debug("No property found for $key")
        return null
    }

     private static Object navigateCollection(Collection collection, String key) {
         def results = []
         for (Object object: collection){
             final res = getSubPath(object, key)
             if (res) results.add(res)
         }
         if (results.isEmpty() ) {
             log.trace("No property found for $key")
             return null
         }
         // Return a single object if only ine results is found.
         return results.size() == 1 ? results[0] : results
     }

    /**
     * Helper function to convert from FileTime to ISO 8601.
     *
     * @param time File time to convert
     * @return Instant or null in case of not available (null)
     */
    static Instant toDate(FileTime time){
        if (time)
            return Instant.ofEpochMilli(time.toMillis())
        else
            return null
    }

    /**
     * Helper function to convert from String ISO 8601 to FileTime.
     *
     * @param date ISO formated time
     * @return Converted FileTime or null if date is not available (null or 'N/A')
     */
    static FileTime toFileTime(String date){
        if (!date)
            return null
        return FileTime.from(Instant.parse(date))
    }

    /**
     * Helper function to unify the encoding of outputs when querying and navigating the CID pseudoFS.
     * Outputs can include CidSerializable objects, collections or parts of these objects.
     * CidSerializable objects can be encoded with the CidEncoder, but collections or parts of
     * these objects require to extend the GsonEncoder.
     *
     * @param output Output to encode
     * @return Output encoded as a JSON string
     */
    static String encodeSearchOutputs(Object output, boolean prettyPrint) {
        if (output instanceof CidSerializable){
            return new CidEncoder().withPrettyPrint(prettyPrint).encode(output)
        } else {
            return new GsonEncoder<Object>() {}
                .withPrettyPrint(prettyPrint)
                .withSerializeNulls(true)
                .withTypeAdapterFactory(CidEncoder.newCidTypeAdapterFactory())
                .encode(output)
        }
    }
}
