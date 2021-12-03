package nextflow.util

/**
 * Holder for parsed command line options
 *
 * @author Manuele Simi <manuele.simi@gmail.com>
 */
class CmdLineOptionMap {

    final options = new LinkedHashMap<String, List<String>>()
    final static instance = new CmdLineOptionMap()

    protected addOption(String key, String value) {
        if ( !options.containsKey(key) )
            options[key] = new ArrayList<String>()
        options[key].add(value)
    }

    boolean hasMultipleValues(String key) {
        options.containsKey(key) ? options[key].size() > 1 : false
    }

    boolean hasOptions() {
        options.size()
    }

    def getValues(String key) {
        options.containsKey(key) ? options[key] : Collections.emptyList()
    }

    def getFirstValue(String key) {
        getFirstValueOrDefault(key, '')
    }

    boolean exists(String key) {
        options.containsKey(key)
    }

    def getFirstValueOrDefault(String key, String alternative) {
        options.containsKey(key) ? options[key].get(0) : alternative
    }

    static CmdLineOptionMap fromMap(final Map map) {
        def optionMap = new CmdLineOptionMap()
        map.each {
            optionMap.addOption(it.key as String, it.value as String)
        }
        return optionMap
    }

    static CmdLineOptionMap emptyOption() {
        return instance
    }

    @Override
    String toString() {
        def serialized = []

        options.each {
            serialized << "option{${it.key}: ${it.value.each {it}}}"
        }

        return "[${serialized.join(', ')}]"
    }
}
