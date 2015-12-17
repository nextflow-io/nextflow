package nextflow.exception

import groovy.transform.InheritConstructors

/**
 * Exception thrown when a process tries to output a file out located
 * out of it's working directory
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@InheritConstructors
class IllegalFileException extends ProcessException {

}
