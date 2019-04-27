package nextflow.exception

import groovy.transform.InheritConstructors

/**
 * Exception thrown when is included in a script more than
 * a module component with the same name
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@InheritConstructors
class DuplicateModuleIncludeException extends ProcessException {
}
