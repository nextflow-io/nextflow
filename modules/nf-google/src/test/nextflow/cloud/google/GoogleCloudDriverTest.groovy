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

import spock.lang.IgnoreIf
import spock.lang.Shared
import spock.lang.Specification

import com.google.api.services.compute.Compute
import com.google.api.services.compute.model.Image
import com.google.api.services.compute.model.Instance
import com.google.api.services.compute.model.MachineType
import com.google.api.services.compute.model.Operation
import nextflow.Global
import nextflow.cloud.CloudConfig
import nextflow.cloud.LaunchConfig
import nextflow.cloud.types.CloudInstanceStatus
import nextflow.exception.AbortOperationException
/**
 *
 * @author Ólafur Haukur Flygenring <olafurh@wuxinextcode.com>
 */
@SuppressWarnings("UnnecessaryQualifiedReference")
class GoogleCloudDriverTest extends Specification {

    final static String VALID_IMAGE_ID = "centos-cloud/global/images/centos-7-v20180815"
    final static String VALID_INSTANCETYPE = "n1-standard-1"
    final static String INVALID_INSTANCETYPE = "n1337-standard-1337"
    final static String ILLEGAL_NAME = "!\tINVALID\t!"
    final static String VALID_CLUSTER = "legal-label"
    final static String ZONE = "us-central1-f"

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
        //Use a real connection to GCE if credentials are defined
        if(runAgainstGce()) {
            print("Got google credentials. Running tests against GCE.")
            sharedHelper = new GceApiHelper(null,ZONE)
        }
        // else fake it so we make it
        else {
            print("No google credentials found.  Stubbing out GCE.")
            sharedHelper = Spy(GceApiHelper, constructorArgs: ["project", ZONE, Stub(Compute)]) as GceApiHelper
            sharedHelper.lookupMachineType(VALID_INSTANCETYPE) >> {
                new MachineType()
                        .setName(VALID_INSTANCETYPE)
                        .setGuestCpus(1)
                        .setMaximumPersistentDisks(65536)
                        .setMemoryMb(3875)
                        .setImageSpaceGb(10)
            }
            sharedHelper.lookupImage(VALID_IMAGE_ID) >> {
                new Image().setName(VALID_IMAGE_ID)
            }
        }
    }

    def 'should create bash profile properly'() {

        given:
        String bash
        LaunchConfig config
        def driver = new GoogleCloudDriver(Stub(GceApiHelper))

        when:
        config = Mock(LaunchConfig)
        config.getNextflow() >> new CloudConfig.Nextflow([version: '0.21.0'])
        config.getClusterName() >> 'cluster-x'

        bash = driver.scriptBashEnv(config)
        then:
        bash == '''
            export NXF_VER='0.21.0'
            export NXF_MODE='google'
            export NXF_EXECUTOR='ignite'
            export NXF_CLUSTER_JOIN='cloud:google:cluster-x'
            export GOOGLE_APPLICATION_CREDENTIALS=$HOME/.nextflow/gce_credentials.json
            '''
                .stripIndent().leftTrim()

        when:
        config = Mock(LaunchConfig)
        config.getClusterName() >> 'cluster-x'
        config.getNextflow() >> new CloudConfig.Nextflow([version: '0.21.0', trace: 'io.nextflow.TaskProcess'])
        bash = driver.scriptBashEnv(config)
        then:
        bash == '''
            export NXF_VER='0.21.0'
            export NXF_MODE='google'
            export NXF_EXECUTOR='ignite'
            export NXF_CLUSTER_JOIN='cloud:google:cluster-x'
            export NXF_TRACE='io.nextflow.TaskProcess'
            export GOOGLE_APPLICATION_CREDENTIALS=$HOME/.nextflow/gce_credentials.json
            '''
                .stripIndent().leftTrim()


        when:
        config = Mock(LaunchConfig)
        config.getClusterName() >> 'cluster-x'
        config.getNextflow() >> new CloudConfig.Nextflow([version: '0.21.0', options: '-XX:this-and-that'])
        bash = driver.scriptBashEnv(config)
        then:
        bash == '''
            export NXF_VER='0.21.0'
            export NXF_MODE='google'
            export NXF_EXECUTOR='ignite'
            export NXF_CLUSTER_JOIN='cloud:google:cluster-x'
            export NXF_OPTS='-XX:this-and-that'
            export GOOGLE_APPLICATION_CREDENTIALS=$HOME/.nextflow/gce_credentials.json
            '''
                .stripIndent().leftTrim()


        when:
        config = Mock(LaunchConfig)
        config.getNextflow() >> new CloudConfig.Nextflow([version: '0.21.0'])
        config.getClusterName() >> 'cluster-x'
        config.getUserName() >> 'illo'
        config.getSharedStorageId() >> 'ef-78289'
        config.getSharedStorageMount() >> '/mount/my/path'
        bash = driver.scriptBashEnv(config)
        then:
        bash == '''
            export NXF_VER='0.21.0'
            export NXF_MODE='google'
            export NXF_EXECUTOR='ignite'
            export NXF_CLUSTER_JOIN='cloud:google:cluster-x'
            export NXF_WORK='/mount/my/path/illo/work'
            export NXF_ASSETS='/mount/my/path/illo/projects'
            export GOOGLE_APPLICATION_CREDENTIALS=$HOME/.nextflow/gce_credentials.json
            '''
                .stripIndent().leftTrim()
    }

    def 'should fill up the cloud-boot template'() {

        given:
        def CONFIG = [
                gce:   [project: 'project', zone: 'zone'],
                cloud: [
                        imageId             : 'image',
                        instanceType        : 'instance',
                        preemptible         : true,
                        userName            : "testUser",
                        createUser          : true,
                        keyFile             : "keyFile",
                        keyHash             : "hash",
                        nextflow: [version: '0.23.0']
                ]]

        def driver = new GoogleCloudDriver(Stub(GceApiHelper))

        def config = CloudConfig.create(CONFIG)
                .setClusterName("my-cluster")
                .setRole('master')
                .build()

        when:
        def script = driver.cloudInitScript(config)
        then:
        script == '''
                    #!/bin/bash

                    su - testUser << 'EndOfScript'
                    (
                    set -e
                    set -x

                    mkdir -p $HOME/.nextflow
                    cat <<EOF >> $HOME/.nextflow/config
                    cloud {
                    \timageId='image'
                    \tinstanceType='instance'
                    \tpreemptible=true
                    \tuserName='testUser'
                    \tcreateUser=true
                    \tkeyFile='keyFile'
                    \tkeyHash='hash'
                    \tnextflow.version='0.23.0'
                    \tclusterName='my-cluster'
                    }

                    EOF

                    mkdir -p $HOME/bin
                    profile=$(mktemp)
                    cat <<EOF > $profile
                    export NXF_VER='0.23.0'
                    export NXF_MODE='google'
                    export NXF_EXECUTOR='ignite'
                    export NXF_CLUSTER_JOIN='cloud:google:my-cluster'
                    export GOOGLE_APPLICATION_CREDENTIALS=$HOME/.nextflow/gce_credentials.json

                    EOF

                    source $profile
                    cat $profile >> $HOME/.bash_profile
                    rm $profile

                    #
                    # Create credentials for Google cloud
                    #
                    if [[ -n "$HOME/.nextflow/gce_credentials.json" ]]; then
                    cat <<EOF >$HOME/.nextflow/gce_credentials.json

                    EOF
                    fi

                    #
                    # Launch docker and pull the container when DOCKER variable is defined
                    #
                    [[ '' ]] && for x in ''; do docker pull $x || true; done

                    #
                    # Install NEXTFLOW and launch it
                    #
                    version="v0.23.0"
                    
                    download_url=${NEXTFLOW_DOWNLOAD_URL:-http://www.nextflow.io/releases}
                    curl -fsSL ${download_url}/${version}/nextflow  > $HOME/nextflow
                    chmod +x $HOME/nextflow
                    $HOME/nextflow -download

                    # pull the nextflow pipeline repo
                    [[ '' ]] && $HOME/nextflow pull ''

                    # launch the nextflow daemon
                    if [[ 'master' == worker ]]; then
                      $HOME/nextflow node -cluster.join "$NXF_CLUSTER_JOIN" -cluster.interface eth0 -bg
                    fi

                    # save the environment for debugging
                    env | sort > boot.env

                    # just a marker file
                    touch READY

                    ) &> ~testUser/boot.log
                    EndOfScript
                    '''
                .stripIndent().leftTrim()

    }

    def 'should include docker pull line'() {
        given:
        def config = [
                gce:   [project: 'project', zone: 'zone'],
                cloud: [
                        dockerPull   : ['cbcrg/image1:tag', 'cbcrg/image2:tag']
                ]]

        def driver = new GoogleCloudDriver(Stub(GceApiHelper))

        when:
        def script = driver.cloudInitScript(CloudConfig.create(config).setClusterName('my-cluster'))
        then:
        script.contains("[[ 'cbcrg/image1:tag cbcrg/image2:tag' ]] && for x in 'cbcrg/image1:tag cbcrg/image2:tag'; do docker pull \$x || true; done")

    }

    def 'should return a valid Compute service'() {
        given:
        GoogleCloudDriver driver = new GoogleCloudDriver(sharedHelper)

        when:
        Compute client = driver.getClient()
        then:
        client in Compute

    }

    def 'should validate a config'() {
        given:
        GoogleCloudDriver driver = new GoogleCloudDriver(sharedHelper)
        def cfg

        when: 'ImageId is missing'
        cfg = CloudConfig.create(cloud: [instanceType: VALID_INSTANCETYPE])
        driver.validate(cfg)
        then:
        thrown AbortOperationException

        when: 'ImageId is invalid'
        cfg = CloudConfig.create(cloud: [instanceType: VALID_INSTANCETYPE, imageId: ILLEGAL_NAME])
        driver.validate(cfg)
        then:
        thrown AbortOperationException

        when: 'instanceType is missing'
        cfg = CloudConfig.create(cloud: [imageId: VALID_IMAGE_ID])
        driver.validate(cfg)
        then:
        thrown AbortOperationException

        when: 'illegal instance type name'
        cfg = CloudConfig.create(cloud: [imageId: VALID_IMAGE_ID, instanceType: ILLEGAL_NAME])
        driver.validate(cfg)
        then:
        thrown AbortOperationException

        when: 'invalid instance type'
        cfg = CloudConfig.create(cloud: [imageId: VALID_IMAGE_ID, instanceType: INVALID_INSTANCETYPE])
        driver.validate(cfg)
        then:
        thrown AbortOperationException

        when: 'invalid cluster name'
        cfg = CloudConfig.create(cloud: [imageId: VALID_IMAGE_ID, instanceType: VALID_INSTANCETYPE, clusterName: ILLEGAL_NAME])
        driver.validate(cfg)
        then:
        thrown AbortOperationException

        when: 'configured with a shared storage'
        cfg = CloudConfig.create(cloud: [imageId: VALID_IMAGE_ID, instanceType: VALID_INSTANCETYPE, clusterName: VALID_CLUSTER, sharedStorageId: "12345"])
        driver.validate(cfg)
        then:
        thrown AbortOperationException

        when: 'configured with a shared storage'
        cfg = CloudConfig.create(cloud: [imageId: VALID_IMAGE_ID, instanceType: VALID_INSTANCETYPE, clusterName: VALID_CLUSTER, sharedStorageId: "12345", sharedStorageMount: "12345"])
        driver.validate(cfg)
        then:
        thrown AbortOperationException

        when: 'has unsupported settings'
        cfg = CloudConfig.create(cloud: [imageId: VALID_IMAGE_ID, instanceType: VALID_INSTANCETYPE, clusterName: VALID_CLUSTER, spotPrice: 100])
        driver.validate(cfg)
        then:
        thrown AbortOperationException

        when: 'has a valid config'
        cfg = CloudConfig.create(cloud: [imageId: VALID_IMAGE_ID, instanceType: VALID_INSTANCETYPE, clusterName: VALID_CLUSTER])
        def result = driver.validate(cfg)
        then:
        result == null
    }

    def 'should describe an instanceType'() {
        given:
        GoogleCloudDriver driver = new GoogleCloudDriver(sharedHelper)

        when: 'valid instanceType'
        def validResult = driver.describeInstanceType(VALID_INSTANCETYPE)
        then:
        validResult?.id == VALID_INSTANCETYPE
        validResult?.cpus > 0

        when: 'invalid instanceType'
        def invalidResult = driver.describeInstanceType(INVALID_INSTANCETYPE)
        then:
        invalidResult == null
    }

    @IgnoreIf({!GoogleCloudDriverTest.runAgainstGce()})
    def 'should intitialize correctly when given google credentials'() {
        given:
        def config = [ project: "projectname", zone: "projectzone"]
        when:
        def driver = new GoogleCloudDriver(config)
        then:
        driver.helper.project == "projectname"
        driver.helper.zone == "projectzone"
    }

    @IgnoreIf({!GoogleCloudDriverTest.runAgainstGce()})
    def 'should initialize with values form the gobal config if given an empty map'() {
        given:
        def globalConfig = [ google: [project: "globalprojectname", zone: "globalprojectzone"]]
        def localCOnfig = [:]
        Global.config = globalConfig
        when:
        def driver = new GoogleCloudDriver(localCOnfig)
        then:
        driver.helper.project == "globalprojectname"
        driver.helper.zone == "globalprojectzone"
    }

    def 'should return the name of instances when asked to launch new instances'() {
        given:
        def helper = Spy(GceApiHelper, constructorArgs: ["project","zone",Stub(Compute)]) {
            if(!runAgainstGce()) {
                getCredentialsFile() >> "stubbed out"
            }
        }

        def driver = new GoogleCloudDriver(helper)

        def cloudConfig = CloudConfig.create([cloud: [keyHash: "hash"]]).build()
        cloudConfig.setClusterName("testcluster")
        cloudConfig.setPreemptible(true)
        when:
        def instances = driver.launchInstances(2,cloudConfig)
        then:
        1 * helper.blockUntilComplete(_ as Iterable<Operation>,_,_)
        2 * helper.blockUntilComplete(_ as Operation,_,_) >> null

        instances.size() == 2
        instances.each {
            assert it.startsWith("testcluster")
        }
    }

    def 'should wait correctly for instance to return a status'() {
        given:
        GceApiHelper helper = Mock()
        GoogleCloudDriver driver = Spy(GoogleCloudDriver, constructorArgs: [helper])

        when: 'wait correctly for READY status from instances'
        driver.waitInstanceStatus(["1"],CloudInstanceStatus.READY)
        then:
        1 * helper.getInstanceList('(name = "1")') >>> [[new Instance().setStatus("RUNNING").setName("1")]]

        when: 'wait correctly for STARTED status from instances'
        driver.waitInstanceStatus(["2"],CloudInstanceStatus.STARTED)
        then:
        1 * helper.getInstanceList('(name = "2")') >>> [[new Instance().setStatus("RUNNING").setName("2")]]

        when: 'wait correctly for TERMINATED status from instances'
        driver.waitInstanceStatus(["3"],CloudInstanceStatus.TERMINATED)
        then:
        1 * helper.getInstanceList('(name = "3")') >>> [[new Instance().setStatus("TERMINATED").setName("3")]]

        when: 'wait correctly for TERMINATED status from instances when instance list is empty'
        driver.waitInstanceStatus(["1","2"],CloudInstanceStatus.TERMINATED)
        then:
        1 * helper.getInstanceList('(name = "1") OR (name = "2")') >> {[] as List<Instance>}
    }

    def 'should run closure against instances'() {
        given:
        def inst1 = new Instance().setName("Inst1").setLabels(["a" : "b"])
        def inst2 = new Instance().setName("Inst2").setLabels(["a" : "b"])
        List<Instance> instances = [inst1,inst2]
        def allCount = 0
        GceApiHelper helper = GroovyMock()
        def driver = new GoogleCloudDriver(helper)


        when: 'should run against all instances'
        driver.eachInstance {allCount++}

        then:
        1 * helper.getInstanceList(_) >> {instances}
        allCount == instances.size()

        when: 'filter, tag and Id for each should also call getInstanceList'
        driver.eachInstanceWithTags(["a" : "a"]){}
        driver.eachInstanceWithFilter("filter"){}
        driver.eachInstanceWithIds(["1"]){}

        then:
        1 * helper.getInstanceList('(labels.a = "a")')
        1 * helper.getInstanceList("filter")
        1 * helper.getInstanceList('(name = "1")')

    }
}