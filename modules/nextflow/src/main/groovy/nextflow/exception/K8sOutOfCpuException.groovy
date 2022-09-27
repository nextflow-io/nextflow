package nextflow.exception

import groovy.transform.InheritConstructors

@InheritConstructors
class K8sOutOfCpuException extends RuntimeException implements RetriableException {
}
