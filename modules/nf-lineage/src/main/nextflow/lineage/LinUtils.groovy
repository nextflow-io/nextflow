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
 */

package nextflow.lineage

import static nextflow.lineage.fs.LinFileSystemProvider.*
import static nextflow.lineage.fs.LinPath.*

import java.nio.file.attribute.FileTime
import java.time.OffsetDateTime
import java.time.ZoneId

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.lineage.model.TaskRun
import nextflow.lineage.model.WorkflowRun
import nextflow.lineage.serde.LinEncoder
import nextflow.lineage.serde.LinSerializable
import nextflow.serde.gson.GsonEncoder
/**
 * Utils class for Lineage IDs.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class LinUtils {

    private static final String[] EMPTY_ARRAY = new String[] {}

    /**
     * Get a metadata lineage object or fragment from the Lineage store.
     *
     * @param store Lineage store.
     * @param uri Object or fragment to retrieve in URI-like format.
     *      Format 'lid://<key>[#fragment]' where:
     *       - Key: Metadata Element key
     *       - Fragment: Element fragment to retrieve.
     * @return Lineage metadata object or fragment.
     */
    static Object getMetadataObject(LinStore store, URI uri) {
        if( uri.scheme != SCHEME ) {
            throw new IllegalArgumentException("Invalid LID URI - scheme is different for $SCHEME")
        }
        final key = uri.authority ? uri.authority + uri.path : uri.path
        if( key == SEPARATOR ) {
            throw new IllegalArgumentException("Cannot get object from the root LID URI")
        }
        if ( uri.query )
            log.warn("Query string is not supported the Linage URI ($uri). It will be ignored.")

        final children = parseChildrenFromFragment(uri.fragment)
        return getMetadataObject0(store, key, children )
    }

    private static Object getMetadataObject0(LinStore store, String key, String[] children = []) {
        final object = store.load(key)
        if (!object) {
            throw new FileNotFoundException("Lineage object $key not found")
        }
        if (children && children.size() > 0) {
            return getSubObject(store, key, object, children)
        }
        return object
    }

    /**
     * Get the array of the search path children elements from the fragment string
     *
     * @param fragment String containing the elements separated by '.'
     * @return array with the parsed element
     */
    static String[] parseChildrenFromFragment(String fragment) {
        if( !fragment )
            return EMPTY_ARRAY
        final children = fragment.tokenize('.')
        new LinPropertyValidator().validate(children)
        return children as String[]
    }

    /**
     * Get a metadata sub-object.
     *
     * If the requested sub-object is the workflow or task outputs, retrieves the outputs from the outputs description.
     *
     * @param store Store to retrieve lineage metadata objects.
     * @param key Parent metadata key.
     * @param object Parent object.
     * @param children Array of string in indicating the properties to navigate to get the sub-object.
     * @return Sub-object or null in it does not exist.
     */
    static Object getSubObject(LinStore store, String key, LinSerializable object, String[] children) {
        if( isSearchingOutputs(object, children) ) {
            // When asking for a Workflow or task output retrieve the outputs description
            final outputs = store.load("${key}#output")
            if (!outputs)
                return null
            return navigate(outputs, children.join('.'))
        }
        return navigate(object, children.join('.'))
    }

    /**
     * Check if the Lid pseudo path or query is for Task or Workflow outputs.
     *
     * @param object Parent Lid metadata object
     * @param children Array of string in indicating the properties to navigate to get the sub-object.
     * @return return 'true' if the parent is a Task/Workflow run and the first element in children is 'outputs'. Otherwise 'false'
     */
    static boolean isSearchingOutputs(LinSerializable object, String[] children) {
        return (object instanceof WorkflowRun || object instanceof TaskRun) && children && children[0] == 'output'
    }

    /**
     * Evaluates object or the objects in a collection matches a set of parameter-value pairs. It includes in the results collection in case of match.
     *
     * @param object Object or collection of objects to evaluate
     * @param params parameter-value pairs to evaluate in each object
     * @param results results collection to include the matching objects
     */
    protected static void treatObject(def object, Map<String, List<String>> params, List<Object> results) {
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
     * Check if an object fulfill the parameter-value
     * 
     * @param object Object to evaluate
     * @param params parameter-value pairs to evaluate
     * @return true if all object parameters exist and matches with the value, otherwise false.
     */
    static boolean checkParams(Object object, Map<String, List<String>> params) {
        for( final entry : params.entrySet() ) {
            final value = navigate(object, entry.key)
            if( !checkParam(value, entry.value) ) {
                return false
            }
        }
        return true
    }

    private static boolean checkParam(Object value, List<String> expected) {
        if( !value )
            return false

        // If value collection, convert to String and check all expected values are in the value.
        if( value instanceof Collection ) {
            final colValue = value as Collection
            return colValue.collect { it.toString() }.containsAll(expected)
        }

        //Single object can't be compared with collection with one of more elements
        if( expected.size() > 1 ) {
            return false
        }

        return value.toString() == expected[0]
    }

    /**
     * Retrieves the sub-object or value indicated by a path.
     *
     * @param obj Object to navigate
     * @param path Elements path separated by '.' e.g. field.subfield
     * @return sub-object / value
     */
    static Object navigate(Object obj, String path) {
        if (!obj)
            return null
        // type has been replaced by class when evaluating LidSerializable objects
        if (obj instanceof LinSerializable && path == 'type')
            return obj.getClass()?.simpleName
        try {
            return path.tokenize('.').inject(obj) { current, key ->
                getSubPath(current, key)
            }
        }
        catch (Throwable e) {
            log.debug("Error navigating to $path in object", e)
            return null
        }
    }

    private static Object getSubPath(current, String key) {
        if (current == null) {
            return null
        }
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
        final results = []
        for (Object object : collection) {
            final res = getSubPath(object, key)
            if (res)
                results.add(res)
        }
        if (results.isEmpty() ) {
            log.trace("No property found for $key")
            return null
        }
        // Return a single object if only ine results is found.
        return results.size() == 1 ? results[0] : results
    }

    /**
     * Helper function to convert from FileTime to ISO 8601 with offset
     * of current timezone.
     *
     * @param time File time to convert
     * @return The {@link OffsetDateTime} for the corresponding file time or null in case of not available (null)
     */
    static OffsetDateTime toDate(FileTime time) {
        return time != null
            ? time.toInstant().atZone(ZoneId.systemDefault()).toOffsetDateTime()
            : null
    }

    /**
     * Helper function to convert from String ISO 8601 to FileTime.
     *
     * @param date ISO formated time
     * @return Converted FileTime or null if date is not available (null or 'N/A')
     */
    static FileTime toFileTime(OffsetDateTime date) {
        if (!date)
            return null
        return FileTime.from(date.toInstant())
    }

    /**
     * Helper function to unify the encoding of outputs when querying and navigating the lineage pseudoFS.
     * Outputs can include LinSerializable objects, collections or parts of these objects.
     * LinSerializable objects can be encoded with the LinEncoder, but collections or parts of
     * these objects require to extend the GsonEncoder.
     *
     * @param output Output to encode
     * @return Output encoded as a JSON string
     */
    static String encodeSearchOutputs(Object output, boolean prettyPrint = false) {
        if (output instanceof LinSerializable) {
            return new LinEncoder().withPrettyPrint(prettyPrint).encode(output)
        } else {
            return new GsonEncoder<Object>() {}
                .withPrettyPrint(prettyPrint)
                .withSerializeNulls(true)
                .withTypeAdapterFactory(LinEncoder.newLidTypeAdapterFactory())
                .encode(output)
        }
    }
}
