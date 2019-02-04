/*
 * Copyright 2018, WuxiNextcode
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nextflow.cloud.google

import com.google.api.services.compute.Compute
import com.google.api.services.compute.model.*
import nextflow.exception.AbortOperationException
import spock.lang.Ignore
import spock.lang.IgnoreIf
import spock.lang.Shared
import spock.lang.Specification

@SuppressWarnings("UnnecessaryQualifiedReference")
class GceApiHelperTest extends Specification {

    static String testProject = "testProject"
    static String testZone = "testZone"

    static boolean runAgainstGce() {
        def path = System.getenv("GOOGLE_APPLICATION_CREDENTIALS")
        if( !path ) return false
        def exists = new File(path).exists()
        if( exists ) return true
        println "Google credentials file is missing: $path"
        return false
    }

    @Shared
    GceApiHelper sharedHelper

    def setupSpec() {

        sharedHelper = runAgainstGce() ?
                Spy(GceApiHelper, constructorArgs: [null,"us-central1-f"]) as GceApiHelper  :
                Spy(GceApiHelper, constructorArgs: [testProject,testZone,Stub(Compute)]) as GceApiHelper

        sharedHelper.readGoogleMetadata(_) >> "metadata"
    }

    @Ignore
    def 'should report error when region and zone are null'() {
        when:
        new GceApiHelper(null,null)
        then:
        thrown(AbortOperationException)
    }

    @IgnoreIf({GceApiHelperTest.runAgainstGce()})
    @Ignore
    //If we have a google credentials file, we can read the project name from it
    def 'should report error when project is missing in initialization'() {
        when:
        new GceApiHelper(null,testZone)
        then:
        thrown(AbortOperationException)
    }

    @Ignore
    def 'should report error when zone is missing in initialization'() {
        when:
        new GceApiHelper(testProject,null)
        then:
        thrown(AbortOperationException)
    }

    def 'should read metadata if it is available'() {
        when:
        def project = sharedHelper.readProject()

        then:
        if(runAgainstGce())
            assert project != "metadata" //should get a real project name here
        else
            assert project == "metadata" //should get the stubbed out value

        when:
        def zone = sharedHelper.readZone()
        then:
        zone == "metadata"

        when:
        def instanceId = sharedHelper.readInstanceId()
        then:
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

    //This is, of course, by no means a definite test that ensures that all names will be unique, but it will hopefully catch silly mistakes
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

    def 'should attach scripts as metadata to an instance'() {
        given:
        Instance gceInstance = new Instance()
        when:
        sharedHelper.setStartupScript(gceInstance,"startup")
        sharedHelper.setShutdownScript(gceInstance,"shutdown")
        def metadata = gceInstance.getMetadata()
        then:
        metadata.getItems().any{it.getKey() == "startup-script" && it.getValue() == "startup"}
        metadata.getItems().any{it.getKey() == "shutdown-script" && it.getValue() == "shutdown"}

    }

    def 'should block until a GCE operation returns status "DONE"'() {
        given:
        Compute.GlobalOperations globalOperations = Mock()
        Compute compute = Mock(Compute)

        def runningOp = new Operation().setStatus("RUNNING")
        def doneOp = new Operation().setStatus("DONE")

        Compute.GlobalOperations.Get computeGlobalOperations = Stub()
        computeGlobalOperations.execute() >>> [runningOp,doneOp]

        GceApiHelper helper = Spy(GceApiHelper,constructorArgs: [testProject,testZone,compute])

        when:
        def ret = helper.blockUntilComplete(runningOp,100,10)
        then:
        (1.._) * globalOperations.get(_,_) >>{
            computeGlobalOperations
        }
        (1.._) * compute.globalOperations() >> {
            globalOperations
        }
        !ret
    }

    def 'should block until multiple GCE operations returns status "DONE"'() {
        given:
        Compute.GlobalOperations globalOperations = Mock()
        Compute compute = Mock(Compute)

        def runningOp = new Operation().setStatus("RUNNING")
        def doneOp = new Operation().setStatus("DONE")

        Compute.GlobalOperations.Get computeGlobalOperations = Stub()
        computeGlobalOperations.execute() >>> [runningOp,doneOp,runningOp,doneOp]

        GceApiHelper helper = Spy(GceApiHelper,constructorArgs: [testProject,testZone,compute])

        when:
        def ret = helper.blockUntilComplete([runningOp,runningOp],100,10)
        then:
        (1.._) * globalOperations.get(_,_) >>{
            computeGlobalOperations
        }
        (1.._) * compute.globalOperations() >> {
            globalOperations
        }
        !ret
    }

    def 'should timeout while waiting too long for an operation to complete'() {
        given:
        Compute.GlobalOperations globalOperations = Mock()
        Compute compute = Mock(Compute)

        def runningOp = new Operation().setStatus("RUNNING")

        Compute.GlobalOperations.Get computeGlobalOperations = Stub()
        computeGlobalOperations.execute() >> {runningOp}

        GceApiHelper helper = Spy(GceApiHelper,constructorArgs: [testProject,testZone,compute])

        when:
        helper.blockUntilComplete(runningOp,100,10)
        then:
        (1.._) * globalOperations.get(_,_) >>{
            computeGlobalOperations
        }
        (1.._) * compute.globalOperations() >> {
            globalOperations
        }
        thrown(InterruptedException)
    }

    def 'should timeout while waiting too long for multiple operations to complete'() {
        given:
        Compute.GlobalOperations globalOperations = Mock()
        Compute compute = Mock(Compute)

        def runningOp = new Operation().setStatus("RUNNING")
        def doneOp = new Operation().setStatus("DONE")

        Compute.GlobalOperations.Get computeGlobalOperations = Stub()
        computeGlobalOperations.execute() >>> [runningOp,doneOp,runningOp,runningOp]

        GceApiHelper helper = Spy(GceApiHelper,constructorArgs: [testProject,testZone,compute])

        when:
        helper.blockUntilComplete([runningOp,runningOp],100,50)
        then:
        (1.._) * globalOperations.get(_,_) >>{
            computeGlobalOperations
        }
        (1.._) * compute.globalOperations() >> {
            globalOperations
        }
        thrown(InterruptedException)
    }


    @IgnoreIf({!GceApiHelperTest.runAgainstGce()})
    def 'should get the google credentials file'() {
        when:
        def fileContent = sharedHelper.getCredentialsFile()
        then:
        fileContent
    }

    @IgnoreIf({!GceApiHelperTest.runAgainstGce()})
    def 'should read project name from google credentials file'() {
        when:
        def project = sharedHelper.readProject()
        then:
        project != "metadata"
    }

    def 'should get a list of instances from gce api'() {
        given:
        Compute compute = Mock()
        Compute.Instances instances = Mock()
        Compute.Instances.List list = Mock()
        GceApiHelper helper = new GceApiHelper("testProject","testZone",compute)
        when:
        def instList = helper.getInstanceList("")
        then:
        1 * compute.instances() >> {instances}
        1 * instances.list(_,_) >> {list}
        1 * list.execute() >>{
            new InstanceList().setItems([new Instance().setName("instance")])
        }

        instList.size() == 1
        instList.get(0).getName() == "instance"
    }

    @IgnoreIf({!GceApiHelperTest.runAgainstGce()})
    def 'should return an Image details'() {
        given:
        def imageId = "centos-cloud/global/images/centos-7-v20180815"

        when:
        def image = sharedHelper.lookupImage(imageId)

        then:
        image.getSelfLink() == sharedHelper.imageName(imageId)

    }


}
