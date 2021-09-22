package nextflow.cloud.azure.batch

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzTaskKeyTest extends Specification {

    def 'should return the key pair' () {
        expect:
        new AzTaskKey('foo','bar').keyPair() == 'foo/bar'
    }

}
