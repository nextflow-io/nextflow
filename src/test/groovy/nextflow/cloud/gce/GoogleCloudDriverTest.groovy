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

import com.amazonaws.services.ec2.AmazonEC2Client
import com.amazonaws.services.ec2.model.*
import com.amazonaws.services.ec2.waiters.AmazonEC2Waiters
import com.amazonaws.waiters.Waiter
import com.google.api.services.compute.Compute
import com.google.api.services.compute.model.MachineType
import nextflow.cloud.CloudConfig
import nextflow.cloud.LaunchConfig
import nextflow.cloud.aws.AmazonCloudDriver
import nextflow.cloud.types.CloudInstanceStatus
import nextflow.cloud.types.CloudInstanceType
import nextflow.exception.AbortOperationException
import nextflow.util.MemoryUnit
import spock.lang.Ignore
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Ólafur Haukur Flygenring <olafurh@wuxinextcode.com>
 */
class GoogleCloudDriverTest extends Specification {

    final static String VALID_IMAGEID = "VALID_IMAGEID"
    final static String VALID_INSTANCETYPE = "VALID_INSTANCETYPE"
    final static String INVALID = "!\tINVALID\t!"
    final static String VALID_CLUSTER = "legal-label"

    @Shared
    GceApiHelper stubbedHelper


    def setupSpec() {
        stubbedHelper = Spy(GceApiHelper,constructorArgs: ["","",Stub(Compute)])
        stubbedHelper.lookupMachineType(VALID_INSTANCETYPE) >> {new MachineType().setName(VALID_INSTANCETYPE)}
        stubbedHelper.lookupImage(VALID_IMAGEID) >> {new com.google.api.services.compute.model.Image().setName(VALID_IMAGEID)}
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
                      master=$(su - testUser ./nextflow cloud list my-cluster | grep master | sed 's/ .*//')
                      mount $master:/$userhome/work $userhome/work
                    fi
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
        GoogleCloudDriver driver = new GoogleCloudDriver(stubbedHelper)

        when:
        Compute client = driver.getClient()
        then:
        client in Compute

    }

    //TODO: Mock out the need to have a working google credentials for testing
    //TODO: Split up into multiple :and blocks
    def 'should validate a config'() {
        given:
        GoogleCloudDriver driver = new GoogleCloudDriver(stubbedHelper)
        def cfg

        when: 'ImageId is missing'
        cfg = CloudConfig.create(cloud: [instanceType: VALID_INSTANCETYPE])
        driver.validate(cfg)
        then:
        thrown AbortOperationException

        when: 'ImageId is invalid'
        cfg = CloudConfig.create(cloud: [instanceType: VALID_INSTANCETYPE, imageId: INVALID])
        driver.validate(cfg)
        then:
        thrown AbortOperationException

        when: 'instanceType is missing'
        cfg = CloudConfig.create(cloud: [imageId: VALID_IMAGEID])
        driver.validate(cfg)
        then:
        thrown AbortOperationException

        when: 'invalid instance type'
        cfg = CloudConfig.create(cloud: [imageId: VALID_IMAGEID, instanceType: INVALID])
        driver.validate(cfg)
        then:
        thrown AbortOperationException

        when: 'invalid cluster name'
        cfg = CloudConfig.create(cloud: [imageId: VALID_IMAGEID, instanceType: VALID_INSTANCETYPE, clusterName: INVALID])
        driver.validate(cfg)
        then:
        thrown AbortOperationException

        when: 'configured with a shared storage'
        cfg = CloudConfig.create(cloud: [imageId: VALID_IMAGEID, instanceType: VALID_INSTANCETYPE, clusterName: VALID_CLUSTER, sharedStorageId: "12345"])
        driver.validate(cfg)
        then:
        thrown AbortOperationException

        when: 'configured with a shared storage'
        cfg = CloudConfig.create(cloud: [imageId: VALID_IMAGEID, instanceType: VALID_INSTANCETYPE, clusterName: VALID_CLUSTER, sharedStorageId: "12345", sharedStorageMount: "12345"])
        driver.validate(cfg)
        then:
        thrown AbortOperationException

        when: 'has a valid config'
        cfg = CloudConfig.create(cloud: [imageId: VALID_IMAGEID, instanceType: VALID_INSTANCETYPE, clusterName: VALID_CLUSTER])
        def result = driver.validate(cfg)
        then:
        result == null
    }
}

