package nextflow.exception

import groovy.transform.InheritConstructors

@InheritConstructors
class K8sOutOfMemoryException extends RuntimeException implements RetriableException {
}
