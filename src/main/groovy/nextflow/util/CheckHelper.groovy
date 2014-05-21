package nextflow.util

/**
 * Validation helper
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CheckHelper {


    static void checkParamsMap( String name, Map<String,?> params, Map<String,List<?>> valid )  {
        def keys = valid.keySet()
        params?.each {
            if( !keys.contains(it.key) )
                throw new IllegalArgumentException("Unknown argument '${it.key}' for operator '$name' -- Possible arguments: ${keys.join(', ')}")

            def values = valid.get(it.key)
            if( values && !values.contains(it.value))
                throw new IllegalArgumentException("Value '${it.value}' cannot be used in in parameter '${it.key}' for operator '$name' -- Possible values: ${values.join(', ')}")
        }
    }

}
