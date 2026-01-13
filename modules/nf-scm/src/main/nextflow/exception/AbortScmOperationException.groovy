package nextflow.exception

import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors

/**
 * Abort the SCM operation execution
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera,io>
 */
@InheritConstructors
@CompileStatic
class AbortScmOperationException extends RuntimeException{
}
