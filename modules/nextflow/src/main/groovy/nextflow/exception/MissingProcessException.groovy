package nextflow.exception

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.script.ScriptMeta

/**
 * Exception thrown when trying to invoke a process or function not existing
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class MissingProcessException extends ProcessException {

    MissingProcessException(ScriptMeta meta, MissingMethodException e) {
        super(message(meta, e.method), e)
    }

    @PackageScope
    static String message(ScriptMeta meta, String name) {
        def result = "Missing process or function with name '$name'"
        def matches = meta.getAllNames().closest(name)
        if( matches.size()==1 )
            result += " -- Did you mean '${matches[0]}' instead?"
        else if( matches.size()>1 )
            result += "\n\nDid you mean any of these instead?\n" + matches.collect { "  $it"}.join('\n') + '\n'
        return result
    }

}
