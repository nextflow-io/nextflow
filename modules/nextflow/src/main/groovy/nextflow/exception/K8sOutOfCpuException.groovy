package nextflow.exception

import groovy.transform.InheritConstructors

/**
 * Exception thrown when pod cannot be scheduled due to k8s `OutOfmemory` error reason
 */
@InheritConstructors
class K8sOutOfCpuException extends RuntimeException implements ProcessRetryableException {
}
