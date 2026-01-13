/*
 * Copyright 2021, Sage-Bionetworks
 *
 * All Rights reserved
 *
 */

package nextflow.secret

import groovy.transform.InheritConstructors
import nextflow.exception.AbortOperationException

/**
 * Thrown when a config secret name could not be resolved
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@InheritConstructors
class MissingSecretException extends AbortOperationException {
}
