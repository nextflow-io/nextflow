package nextflow.exception

import groovy.transform.InheritConstructors

/**
 * Exception thrown when pod cannot be scheduled due to k8s `OutOfcpu` error reason
 */
@InheritConstructors
class K8sOutOfMemoryException extends RuntimeException implements ProcessRetryableException {
}
