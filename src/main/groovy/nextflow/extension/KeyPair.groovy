package nextflow.extension

import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Implements an helper key-value helper object used in dataflow operators
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString
@EqualsAndHashCode
class KeyPair {
    List keys
    List values
}