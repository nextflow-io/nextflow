/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.cloud.gce

import com.google.api.services.compute.Compute
import com.google.api.services.compute.model.Image
import com.google.api.services.compute.model.Instance
import com.google.api.services.compute.model.MachineType
import com.google.api.services.compute.model.Operation
import nextflow.Global
import nextflow.cloud.CloudConfig
import nextflow.cloud.LaunchConfig
import nextflow.cloud.types.CloudInstance
import nextflow.cloud.types.CloudInstanceStatus
import nextflow.exception.AbortOperationException
import spock.lang.IgnoreIf
import spock.lang.Shared
import spock.lang.Specification

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
        System.getenv("GOOGLE_APPLICATION_CREDENTIALS")
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
            export NXF_MODE='ignite'
            export NXF_EXECUTOR='ignite'
            export NXF_CLUSTER_JOIN='cloud:gce:cluster-x'
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
            export NXF_MODE='ignite'
            export NXF_EXECUTOR='ignite'
            export NXF_CLUSTER_JOIN='cloud:gce:cluster-x'
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
            export NXF_MODE='ignite'
            export NXF_EXECUTOR='ignite'
            export NXF_CLUSTER_JOIN='cloud:gce:cluster-x'
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
            export NXF_MODE='ignite'
            export NXF_EXECUTOR='ignite'
            export NXF_CLUSTER_JOIN='cloud:gce:cluster-x'
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
                    

                    # Install java and docker.
                    # TODO: have this configured as an option in gce config block or remove altogether

                    if ! which java; then
                      if which apt-get; then
                        which add-apt-repository || apt-get -y install software-properties-common
                        which add-apt-repository || apt-get -y install python-software-properties
                        add-apt-repository -y ppa:openjdk-r/ppa
                        apt-get update
                        apt-get install -y openjdk-8-jdk
                      elif which yum; then
                        yum install -y java-1.8.0-openjdk
                      else
                        echo "Neither apt-get nor yum available, cannot install java"
                      fi
                    fi

                    if ! which docker; then
                    \tcurl -fsSL get.docker.com -o get-docker.sh
                    \tsh get-docker.sh

                    fi
                    usermod -aG docker testUser

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
                    export NXF_MODE='ignite'
                    export NXF_EXECUTOR='ignite'
                    export NXF_CLUSTER_JOIN='cloud:gce:my-cluster'
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
                    
                    download_url=${NEXFLOW_DOWNLOAD_URL:-http://www.nextflow.io/releases}
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
                    
                    # NFS share work folder
                    userhome=$(echo ~testUser)
                    sudo -u testUser mkdir -p $userhome/work

                    # TODO: Have automatic NFS support configurable or remove altogether
                    # TODO: apt-get/yum auto
                    yum -y install nfs-utils
                    echo "Userhome: ${userhome}"
                    if [[ 'master' == master ]]; then
                      systemctl enable nfs-server.service
                      systemctl start nfs-server.service
                      echo "${userhome}/work           *(rw,sync,no_root_squash,no_subtree_check)" >/etc/exports
                      exportfs -a
                    else
                      master=$(su - testUser ./nextflow cloud list my-cluster -driver gce | grep master | sed 's/ .*//')
                      mount $master:/$userhome/work $userhome/work
                    fi
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
    def 'should intitialize with values form the gobal config if given an empty map'() {
        given:
        def globalConfig = [ gce: [project: "globalprojectname", zone: "globalprojectzone"]]
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
        1 * driver.waitInstanceStatus(_ ,_)
        1 * helper.getInstanceList(_) >>> [[new Instance().setStatus("RUNNING").setName("1")]]

        when: 'wait correctly for STARTED status from instances'
        driver.waitInstanceStatus(["1","2","3"],CloudInstanceStatus.STARTED)
        then:
        1 * driver.waitInstanceStatus(_ ,_)
        3 * helper.getInstanceList(_) >>> [[new Instance().setStatus("PROVISIONING").setName("1")],[new Instance().setStatus("STAGING").setName("2")],[new Instance().setStatus("RUNNING").setName("3")]]

        when: 'wait correctly for TERMINATED status from instances'
        driver.waitInstanceStatus(["1"],CloudInstanceStatus.TERMINATED)
        then:
        1 * driver.waitInstanceStatus(_ ,_)
        1 * helper.getInstanceList(_) >>> [[new Instance().setStatus("TERMINATED").setName("1")]]

        when: 'wait correctly for TERMINATED status from instances when instance list is empty'
        driver.waitInstanceStatus(["1","2"],CloudInstanceStatus.TERMINATED)
        then:
        1 * driver.waitInstanceStatus(_ ,_)
        1 * helper.getInstanceList(_) >> {[] as List<Instance>}
    }
}