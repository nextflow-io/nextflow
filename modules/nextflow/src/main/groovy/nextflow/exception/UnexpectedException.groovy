package nextflow.exception

import groovy.transform.InheritConstructors

/**
 * Unexpected unrecoverable runtime exception
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@InheritConstructors
class UnexpectedException extends Exception implements ShowOnlyExceptionMessage {
}
