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

/**
 * Utils class for CID.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class CidUtils {

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
                final children = getChildrenFormFragment(uri.fragment)
                return searchPath(store, key, parameters, children )
            }
        } catch(Throwable e){
            log.debug("Exception querying $uri. $e.message")
            return []
        }

    }

    public static String[] getChildrenFormFragment(String fragment) {
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

    private static List<Object> searchPath(CidStore store, String path, Map<String, String> params, String[] children = []) {
        final results = new LinkedList<Object>()
        final object = store.load(path)
        if (object) {
            if (children && children.size() > 0) {
                final output = navigate(object, children.join('.'))
                if (output) {
                    treatObject(output, params, results)
                } else {
                    throw new FileNotFoundException("Cid object $path/${children ? children.join('/') : ''} not found.")
                }
            } else {
                treatObject(object, params, results)
            }
        } else {
            // If there isn't metadata check the parent to check if it is a subfolder of a task/workflow output
            final currentPath = Path.of(path)
            final parent = currentPath.getParent()
            if (parent) {
                ArrayList<String> newChildren = new ArrayList<String>()
                newChildren.add(currentPath.getFileName().toString())
                newChildren.addAll(children)
                return searchPath(store, parent.toString(), params, newChildren as String[])
            } else {
                throw new FileNotFoundException("Cid object $path/${children ? children.join('/') : ''} not found.")
            }
        }
        return results
    }

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

    public static Map<String, String> parseQuery(String queryString) {
        return queryString.split('&').collectEntries {
            it.split('=').collect { URLDecoder.decode(it, 'UTF-8') }
        } as Map<String, String>
    }

    public static boolean checkParams(Object object, Map<String,String> params) {
        for (final entry : params.entrySet()) {
            final value = navigate(object, entry.key)
            if (!value || value.toString() != entry.value.toString() ) {
                return false
            }
        }
        return true
    }

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
}
