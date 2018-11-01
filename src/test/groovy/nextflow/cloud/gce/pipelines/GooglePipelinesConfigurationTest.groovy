package nextflow.cloud.gce.pipelines

import spock.lang.Specification

class GooglePipelinesConfigurationTest extends Specification {

    def 'should construct correctly'() {
        given:
        def name = "testName"
        def testZone = "testZone"
        def instanceType = "testInstanceType"
        def preemp = true
        def remoteBinDir = null

        when:
        def config1 = new GooglePipelinesConfiguration(name,testZone,instanceType,remoteBinDir,preemp)
        def config2 = new GooglePipelinesConfiguration(name,testZone,instanceType)

        then:
        with(config1) {
            project == name
            zone == testZone
            vmInstanceType == instanceType
            preemptible == preemp
        }

        with(config2) {
            project == name
            zone == testZone
            vmInstanceType == instanceType
            !preemptible
        }
    }
}
