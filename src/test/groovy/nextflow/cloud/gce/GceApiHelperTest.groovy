package nextflow.cloud.gce

import com.google.api.services.compute.Compute
import com.google.api.services.compute.model.AttachedDisk
import com.google.api.services.compute.model.NetworkInterface
import nextflow.exception.AbortOperationException
import spock.lang.Shared
import spock.lang.Specification

//TODO: Implement real or stubbed tests depending on if we are running with google credentials or not
class GceApiHelperTest extends Specification {

    static String testProject = "testProject"
    static String testZone = "testZone"

    @Shared
    Compute computeStub = Stub()
    @Shared
    GceApiHelper sharedHelper = new GceApiHelper(testProject,testZone,computeStub)

    def 'should report error when region and zone are null'() {
        when:
        new GceApiHelper(null,null)
        then:
        thrown(AbortOperationException)
    }

    def 'should report error when project is missing in initialization'() {
        when:
        new GceApiHelper(null,testZone)
        then:
        thrown(AbortOperationException)
    }

    def 'should report error when zone is missing in initialization'() {
        when:
        new GceApiHelper(testProject,null)
        then:
        thrown(AbortOperationException)
    }

    def 'should read metadata if it is available'() {
        given:
        def helper = Spy(GceApiHelper, constructorArgs: [testProject,testZone,computeStub])
        helper.readGoogleMetadata(_) >> "metadata"
        when:
        def project = helper.readProject()
        def zone = helper.readZone()
        def instanceId = helper.readInstanceId()
        then:
        project == "metadata"
        zone == "metadata"
        instanceId == "metadata"
    }

    def 'should return a valid boot disk'() {
        when:
        AttachedDisk disk = sharedHelper.createBootDisk("testDisk","testimage")
        then:
        disk.getBoot()
        disk.getInitializeParams().getDiskName() == "testDisk"
        disk.getInitializeParams().getSourceImage() == GceApiHelper.imageName("testimage")
    }

    def 'should return a valid network interface'() {
        when:
        NetworkInterface netInt = sharedHelper.createNetworkInterface()
        then:
        netInt as NetworkInterface
    }

    //This is, of course, by no means a definite test that all names will be unique, but it will hopefully catch silly mistakes
    def 'should return random names'() {
        when:
        def randomNames = []
        1000.times {randomNames << sharedHelper.randomName()}
        def baseRandomNames = []
        1000.times {baseRandomNames << sharedHelper.randomName("test-")}
        then:
        randomNames.every {testString -> randomNames.count {it == testString} == 1}
        baseRandomNames.every {testString -> testString.startsWith("test") && baseRandomNames.count {it == testString} == 1}
        randomNames != baseRandomNames
    }

    def 'should validate label values correctly'() {
        given:
        def tooLong = ""
        64.times {tooLong += "x"}
        def justRight = ""
        63.times {justRight += "x"}

        when:
        def tlres = sharedHelper.validateLabelValue(tooLong)
        def jrres = sharedHelper.validateLabelValue(justRight)
        def illegalres = sharedHelper.validateLabelValue("æ !#%&")
        def legalres = sharedHelper.validateLabelValue("12345abcde")
        then:
        tlres != null
        jrres == null
        illegalres != null
        legalres == null
    }

    def 'should create correct scheduling'() {
        when:
        def nonPre = sharedHelper.createScheduling(false)
        def Pre = sharedHelper.createScheduling(true)
        then:
        !nonPre.getPreemptible()
        Pre.getPreemptible()
    }

    //TODO: Need to discover if there is "official" instructions how to do this correctly
    def 'should convert public ip to google dns name correctly'() {
        when:
        def ret = sharedHelper.publicIpToDns("192.168.1.254")
        then:
        ret == "254.1.168.192.bc.googleusercontent.com"
    }

    def 'should throw an error when converting an illegal public ip'() {
        when:
        sharedHelper.publicIpToDns("12.12")
        then:
        thrown(IllegalArgumentException)
    }

}
