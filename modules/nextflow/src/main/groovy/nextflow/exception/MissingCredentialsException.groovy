package nextflow.exception

import groovy.transform.InheritConstructors

/**
 * Exception thrown when credentials required by a service are not available
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@InheritConstructors
class MissingCredentialsException extends AbortOperationException {
}
