package nextflow.data.cid

import groovy.util.logging.Slf4j
import nextflow.data.cid.fs.CidPath
import nextflow.data.cid.serde.CidSerializable

import java.nio.file.Path

@Slf4j
class CidUtils {

    public static List query(CidStore store, URI uri) {
        String key = uri.authority ? uri.authority + uri.path : uri.path
        try {
            if (key == CidPath.SEPARATOR) {
                return store.search(uri.query)
            } else {
                return searchPath(store, key, uri.query ? parseQuery(uri.query) : null)
            }
        } catch(Throwable e){
            log.debug("Exception querying $uri. $e.message")
            return []
        }

    }

    private static List<Object> searchPath(CidStore store, String path, Map<String, String> params, String[] childs = []) {
        final results = new LinkedList<Object>()
        final object = store.load(path)
        if (object) {
            if (childs && childs.size() > 0) {
                final output = navigate(object, childs.join('.'))
                if (output) {
                    treatObject(output, params, results)
                } else {
                    throw new FileNotFoundException("Cid object $path/${childs ? childs.join('/') : ''} not found.")
                }
            } else {
                treatObject(object, params, results)
            }
        } else {
            // If there isn't metadata check the parent to check if it is a subfolder of a task/workflow output
            final currentPath = Path.of(path)
            final parent = currentPath.getParent()
            if (parent) {
                ArrayList<String> newChilds = new ArrayList<String>()
                newChilds.add(currentPath.getFileName().toString())
                newChilds.addAll(childs)
                return searchPath(store, parent.toString(), params, newChilds as String[])
            } else {
                throw new FileNotFoundException("Cid object $path/${childs ? childs.join('/') : ''} not found.")
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

    protected static Object navigate(Object obj, String path){
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
                log.debug("No property found for $key")
                return null // Property not found
            }
        } catch (Throwable e) {
            log.debug("Error navigating to $path in object", e)
            return null
        }
    }
}
