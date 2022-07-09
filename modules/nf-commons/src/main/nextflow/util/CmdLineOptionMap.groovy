package nextflow.util


import com.google.common.hash.Hasher
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Holder for parsed command line options.
 *
 * @author Manuele Simi <manuele.simi@gmail.com>
 */
@CompileStatic
@ToString(includes = 'options', includeFields = true)
@EqualsAndHashCode(includes = 'options', includeFields = true)
class CmdLineOptionMap implements CacheFunnel {

    final private Map<String, List<String>> options = new LinkedHashMap<String, List<String>>()
    final private static CmdLineOptionMap EMPTY = new CmdLineOptionMap()

    protected CmdLineOptionMap addOption(String key, String value) {
        if ( !options.containsKey(key) )
            options[key] = new ArrayList<String>(10)
        options[key].add(value)
        return this
    }

    boolean hasMultipleValues(String key) {
        options.containsKey(key) ? options[key].size() > 1 : false
    }

    boolean hasOptions() {
        options.size()
    }

    List<String> getValues(String key) {
        return options.containsKey(key) ? options[key] : new ArrayList<String>(10)
    }

    def getFirstValue(String key) {
        getFirstValueOrDefault(key, null)
    }

    boolean asBoolean() {
        return options.size()>0
    }

    boolean exists(String key) {
        options.containsKey(key)
    }

    def getFirstValueOrDefault(String key, String alternative) {
        options.containsKey(key) && options[key].get(0) ? options[key].get(0) : alternative
    }

    static CmdLineOptionMap fromMap(final Map map) {
        def optionMap = new CmdLineOptionMap()
        map.each {
            optionMap.addOption(it.key as String, it.value as String)
        }
        return optionMap
    }

    static CmdLineOptionMap emptyOption() {
        return EMPTY
    }

    @Override
    String toString() {
        def serialized = []
        options.each {
            serialized << "option{${it.key}: ${it.value.each {it}}}"
        }
        return "[${serialized.join(', ')}]"
    }

    @Override
    Hasher funnel(Hasher hasher, CacheHelper.HashMode mode) {
        return CacheHelper.hasher(hasher, options, mode)
    }
}
