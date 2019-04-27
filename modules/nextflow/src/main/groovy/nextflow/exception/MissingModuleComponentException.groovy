package nextflow.exception

import groovy.transform.PackageScope
import nextflow.script.ScriptMeta

/**
 * Exception thrown when a module component cannot be found
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class MissingModuleComponentException extends ProcessException {

    MissingModuleComponentException(ScriptMeta meta, String name ) {
        super(message(meta, name))
    }

    @PackageScope
    static String message(ScriptMeta meta, String name) {
        def result = "Cannot find a component with name '$name' in module: $meta.scriptPath"
        def names = meta.getDefinitions().collect { it.name }
        def matches = names.closest(name)
        result += "\n\nDid you mean one of these?\n" + matches.collect { "  $it"}.join('\n')
        return result
    }
}
