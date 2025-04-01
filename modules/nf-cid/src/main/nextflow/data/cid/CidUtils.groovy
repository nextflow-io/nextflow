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
import nextflow.data.cid.serde.CidSerializable

import java.nio.file.Path
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
    public static List query(CidStore store, URI uri) {
        String key = uri.authority ? uri.authority + uri.path : uri.path
        try {
            if (key == CidPath.SEPARATOR) {
                final results = store.search(uri.query)
                if (results && uri.fragment){
                    // If fragment is defined get the property of the object indicated by the fragment
                    final filteredResults = []
                    results.forEach {
                        final output = navigate(it, uri.fragment)
                        if (output){
                            filteredResults.add(output)
                        }
                    }
                    return filteredResults
                }
                return results
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

    /**
     * Get the array of the search path children elements from the fragment string
     * @param fragment String containing the elements separated by '.'
     * @return array with the parsed element
     */
    public static String[] parseChildrenFormFragment(String fragment) {
        if (fragment) {
            if (fragment.contains('.')) {
                return fragment.split("\\.")
            } else {
                return [fragment] as String[]
            }
        } else {
            return [] as String[]
        }
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
        final results = new LinkedList<Object>()
        final object = store.load(key)
        if (object) {
            if (children && children.size() > 0) {
                final output = navigate(object, children.join('.'))
                if (output) {
                    treatObject(output, params, results)
                } else {
                    throw new FileNotFoundException("Cid object $key/${children ? children.join('/') : ''} not found.")
                }
            } else {
                treatObject(object, params, results)
            }
        } else {
            // If there isn't metadata check the parent to check if it is a subfolder of a task/workflow output
            final currentPath = Path.of(key)
            final parent = currentPath.getParent()
            if (parent) {
                ArrayList<String> newChildren = new ArrayList<String>()
                newChildren.add(currentPath.getFileName().toString())
                newChildren.addAll(children)
                return searchPath(store, parent.toString(), params, newChildren as String[])
            } else {
                throw new FileNotFoundException("Cid object $key/${children ? children.join('/') : ''} not found.")
            }
        }
        return results
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
    public static Map<String, String> parseQuery(String queryString) {
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
    public static boolean checkParams(Object object, Map<String,String> params) {
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
    public static Object navigate(Object obj, String path){
        if (!obj)
            return null
        try{
            // type has been replaced by class when evaluating CidSerializable objects
            if (obj instanceof CidSerializable && path == 'type')
                return obj.getClass()?.simpleName
            path.tokenize('.').inject(obj) { current, key ->
                if (current == null) return null

                if (current instanceof Map) {
                    return current[key] // Navigate Map properties
                }

                if (current.metaClass.hasProperty(current, key)) {
                    return current.getAt(key) // Navigate Object properties
                }
                log.trace("No property found for $key")
                return null // Property not found
            }
        } catch (Throwable e) {
            log.debug("Error navigating to $path in object", e)
            return null
        }
    }

    /**
     * Helper function to convert from FileTime to ISO 8601.
     *
     * @param time File time to convert
     * @return ISO Date format or 'N/A' in case of not available (null)
     */
    public static String toDate(FileTime time){
        if (time)
            return Instant.ofEpochMilli(time.toMillis()).toString()
        else
            return 'N/A'
    }

    /**
     * Helper function to convert from String ISO 8601 to FileTime.
     *
     * @param date ISO formated time
     * @return Converted FileTime or null if date is not available (null or 'N/A')
     */
    public static FileTime toFileTime(String date){
        if (!date || date == 'N/A')
            return null
        return FileTime.from(Instant.parse(date))
    }
}
