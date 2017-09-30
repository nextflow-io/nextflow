/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
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

package nextflow.cloud.aws
import java.nio.file.Files

import com.amazonaws.services.ec2.AmazonEC2Client
import com.amazonaws.services.ec2.model.BlockDeviceMapping
import com.amazonaws.services.ec2.model.EbsBlockDevice
import com.amazonaws.services.ec2.model.GroupIdentifier
import com.amazonaws.services.ec2.model.RequestSpotInstancesRequest
import com.amazonaws.services.ec2.model.RunInstancesRequest
import com.amazonaws.services.ec2.waiters.AmazonEC2Waiters
import com.amazonaws.waiters.Waiter
import nextflow.Global
import nextflow.cloud.CloudConfig
import nextflow.cloud.LaunchConfig
import nextflow.cloud.types.CloudInstanceStatus
import nextflow.cloud.types.CloudInstanceType
import nextflow.config.ConfigBuilder
import nextflow.util.MemoryUnit
import org.apache.commons.lang.SerializationUtils
import spock.lang.Specification
import spock.lang.Unroll
import spock.util.environment.RestoreSystemProperties
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AmazonCloudDriverTest extends Specification {

    def setup() {
        Global.config = [aws: [
            accessKey: 'xxx-alpha',
            secretKey: 'zzz-beta',
            region: 'eu-west-1'
        ]]
    }

    @RestoreSystemProperties
    def 'should create mime-multipart cloud-init' () {

        given:
        def folder = Files.createTempDirectory('test')
        folder.resolve('.ssh').mkdir()
        folder.resolve('.ssh/id_rsa.pub').text = 'ssh-rsa fake'
        System.setProperty('user.home', folder.toString())

        def file = folder.resolve('nextflow.config')
        file.text = '''
            aws {
                accessKey = 'abc'
                secretKey = 'xxx'
                region = 'eu-west-1'
            }

            cloud {
                imageId = 'ami-s8a9s8'
                instanceType = 'r3.xlarge'

                bootStorageSize = '10GB'

                instanceStorageDevice = '/dev/xyz'
                instanceStorageMount = '/scratch'

                sharedStorageId = '78xs8'
                sharedStorageMount = '/mnt/efs'
            }

            '''

        def config = new ConfigBuilder().buildConfig([file])
        def driver = new AmazonCloudDriver(Mock(AmazonEC2Client))

        when:
        def payload = driver.getUserDataAsString( CloudConfig.create(config).build() )
        then:
        payload.contains('MIME-Version: 1.0')
        payload.contains('Content-Type: multipart/mixed;')
        payload.contains('Content-Type: text/x-shellscript; charset=us-ascii')
        payload.contains('Content-Type: text/cloud-boothook; charset=us-ascii')

        cleanup:
        folder?.deleteDir()
    }

    def 'should create bash profile properly' () {

        given:
        def CONFIG = [
                cloud: [
                        nextflow: [version: '0.21.0'],
                ]
            ]
        def cloud = new AmazonCloudDriver(Mock(AmazonEC2Client))

        when:
        def cfg = (Map)SerializationUtils.clone(CONFIG)
        def bash1 = cloud.scriptBashEnv(CloudConfig.create(cfg).setClusterName('cluster-x'))
        then:
        bash1 == '''
            export NXF_VER='0.21.0'
            export NXF_MODE='ignite'
            export NXF_EXECUTOR='ignite'
            export NXF_CLUSTER_JOIN='cloud:aws:cluster-x'
            '''
            .stripIndent() .leftTrim()

        when:
        cfg = (Map)SerializationUtils.clone(CONFIG)
        cfg.cloud.nextflow.trace = 'io.nextflow.TaskProcess'
        def bash2 = cloud.scriptBashEnv(CloudConfig.create(cfg).setClusterName('cluster-x'))
        then:
        bash2 == '''
            export NXF_VER='0.21.0'
            export NXF_MODE='ignite'
            export NXF_EXECUTOR='ignite'
            export NXF_CLUSTER_JOIN='cloud:aws:cluster-x'
            export NXF_TRACE='io.nextflow.TaskProcess'
            '''
                .stripIndent() .leftTrim()

        when:
        cfg = (Map)SerializationUtils.clone(CONFIG)
        cfg.cloud.nextflow.options = '-XX:this-and-that'
        def bash3 = cloud.scriptBashEnv(CloudConfig.create(cfg).setClusterName('cluster-x'))
        then:
        bash3 == '''
            export NXF_VER='0.21.0'
            export NXF_MODE='ignite'
            export NXF_EXECUTOR='ignite'
            export NXF_CLUSTER_JOIN='cloud:aws:cluster-x'
            export NXF_OPTS='-XX:this-and-that'
            '''
                .stripIndent() .leftTrim()


        when:
        cfg = (Map)SerializationUtils.clone(CONFIG)
        cfg.cloud.sharedStorageId = 'ef-78289'
        cfg.cloud.sharedStorageMount = '/mount/my/path'
        cfg.cloud.userName = 'illo'
        def bash4 = cloud.scriptBashEnv(CloudConfig.create(cfg).setClusterName('cluster-x'))
        then:
        bash4 == '''
            export NXF_VER='0.21.0'
            export NXF_MODE='ignite'
            export NXF_EXECUTOR='ignite'
            export NXF_CLUSTER_JOIN='cloud:aws:cluster-x'
            export NXF_WORK='/mount/my/path/illo/work'
            export NXF_ASSETS='/mount/my/path/illo/projects'
            '''
                .stripIndent() .leftTrim()

        when:
        cfg = (Map)SerializationUtils.clone(CONFIG)
        cfg.cloud.instanceStorageDevice = '/dev/abc'
        cfg.cloud.instanceStorageMount = '/scratch'
        def bash5 = cloud.scriptBashEnv(CloudConfig.create(cfg).setClusterName('cluster-x'))
        then:
        bash5 == '''
            export NXF_VER='0.21.0'
            export NXF_MODE='ignite'
            export NXF_EXECUTOR='ignite'
            export NXF_CLUSTER_JOIN='cloud:aws:cluster-x'
            export NXF_TEMP='/scratch'
            '''
                .stripIndent() .leftTrim()


        when:
        cloud.accessKey = 'alpha'
        cloud.secretKey = 'beta'
        cloud.region = 'region-1'
        then:
        def bash6 = cloud.scriptBashEnv(CloudConfig.create(cfg).setClusterName('cluster-x'))
        then:
        bash6 == '''
            export NXF_VER='0.21.0'
            export NXF_MODE='ignite'
            export NXF_EXECUTOR='ignite'
            export NXF_CLUSTER_JOIN='cloud:aws:cluster-x'
            export NXF_TEMP='/scratch'
            export AWS_ACCESS_KEY_ID='alpha'
            export AWS_SECRET_ACCESS_KEY='beta'
            export AWS_DEFAULT_REGION='region-1'
            '''
                .stripIndent() .leftTrim()
    }



    def 'should fill up the cloud-boot template' () {

        given:
        def config = [
                aws: [accessKey: 'xxx', secretKey: 'yyy', region: 'eu-west-1'],
                cloud: [
                imageId: 'ami-123',
                instanceType: 'r3.abc',
                securityGroup: 'sg-df72b9ba',
                keyName: 'eurokey',
                subnetId: 'subnet-05222a43',
                bootStorageSize: '10GB',
                instanceStorageDevice: '/dev/abc',
                instanceStorageMount: '/scratch',
                sharedStorageId: 'fs-1803efd1',
                sharedStorageMount: '/mnt/efs',
                spotPrice: 0.55,
                nextflow: [version: '0.23.0'],
                autoscale: [ enabled: true, spotPrice: 0.33 ]

        ]]

        def driver = new AmazonCloudDriver(Mock(AmazonEC2Client))
        driver.accessKey = 'xxx'
        driver.secretKey = 'yyy'
        driver.region = 'eu-west-1'

        when:
        def cloud = CloudConfig.create(config)
                        .setRole('master')
                        .setClusterName('my-cluster')
                        .build()

        def script = driver.cloudInitScript(cloud)
        then:
        script ==   '''
                    #!/bin/bash
                    su - ec2-user << 'EndOfScript'
                    (
                    set -e
                    set -x

                    mkdir -p $HOME/.nextflow
                    cat <<EOF >> $HOME/.nextflow/config
                    cloud {
                    \timageId='ami-123'
                    \tinstanceType='r3.abc'
                    \tsecurityGroup='sg-df72b9ba'
                    \tkeyName='eurokey'
                    \tsubnetId='subnet-05222a43'
                    \tbootStorageSize='10GB'
                    \tinstanceStorageDevice='/dev/abc'
                    \tinstanceStorageMount='/scratch'
                    \tsharedStorageId='fs-1803efd1'
                    \tsharedStorageMount='/mnt/efs'
                    \tspotPrice=0.55
                    \tnextflow.version='0.23.0'
                    \tautoscale {
                    \t\tenabled=true
                    \t\tspotPrice=0.33
                    \t}
                    \tclusterName='my-cluster'
                    \tuserName='ec2-user'
                    }

                    EOF

                    mkdir -p $HOME/bin
                    profile=$(mktemp)
                    cat <<EOF > $profile
                    export NXF_VER='0.23.0'
                    export NXF_MODE='ignite'
                    export NXF_EXECUTOR='ignite'
                    export NXF_CLUSTER_JOIN='cloud:aws:my-cluster'
                    export NXF_WORK='/mnt/efs/ec2-user/work'
                    export NXF_ASSETS='/mnt/efs/ec2-user/projects'
                    export NXF_TEMP='/scratch'
                    export AWS_ACCESS_KEY_ID='xxx'
                    export AWS_SECRET_ACCESS_KEY='yyy'
                    export AWS_DEFAULT_REGION='eu-west-1'

                    EOF

                    source $profile
                    cat $profile >> $HOME/.bash_profile
                    rm $profile

                    #
                    # set instance name
                    #
                    instance="$(curl -s http://169.254.169.254/latest/meta-data/instance-id)"
                    zone="$(curl -s 169.254.169.254/latest/meta-data/placement/availability-zone)"
                    region="${zone::-1}"

                    #
                    # Launch docker and pull the container when DOCKER variable is defined
                    #
                    [[ '' ]] && for x in ''; do docker pull $x || true; done

                    #
                    # Install NEXTFLOW and launch it
                    #
                    version="v0.23.0"
                    curl -fsSL http://www.nextflow.io/releases/${version}/nextflow  > $HOME/nextflow
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

                    ) &> ~ec2-user/boot.log
                    EndOfScript
                    '''
                    .stripIndent().leftTrim()

    }

    def 'should include docker pull line' () {
        given:
        def config = [
                aws: [accessKey: 'xxx', secretKey: 'yyy', region: 'eu-west-1'],
                cloud: [
                        imageId: 'ami-123',
                        instanceType: 'r3.abc',
                        securityGroup: 'sg-df72b9ba',
                        keyName: 'eurokey',
                        subnetId: 'subnet-05222a43',
                        nextflow: [version: '0.23.0'],
                        dockerPull: ['cbcrg/image1:tag', 'cbcrg/image2:tag']

                ]]


        def driver = new AmazonCloudDriver(Mock(AmazonEC2Client))

        when:
        def script = driver.cloudInitScript(CloudConfig.create(config).setClusterName('my-cluster'))
        then:
        script.contains("[[ 'cbcrg/image1:tag cbcrg/image2:tag' ]] && for x in 'cbcrg/image1:tag cbcrg/image2:tag'; do docker pull \$x || true; done")

    }

    def 'should create script to mount instance storage' () {

        given:
        def driver = [:] as AmazonCloudDriver
        when:
        def snippet = driver.scriptMountInstanceStorage('/dev/xyz', '/mnt/scratch', 'ubuntu')
        then:
        snippet ==  """
                    mkfs.ext4 -E nodiscard /dev/xyz
                    mkdir -p /mnt/scratch
                    mount -o discard /dev/xyz /mnt/scratch
                    chown -R ubuntu:ubuntu /mnt/scratch
                    chmod 775 /mnt/scratch
                    """
                .stripIndent().leftTrim()
    }

    def 'should create script to mount EFS storage' () {

        given:
        def driver = [:] as AmazonCloudDriver
        when:
        def snippet = driver.scriptMountEFS('/dev/xyz', '/mnt/scratch', 'ubuntu')
        then:
        snippet ==  '''
                    zone="$(curl -s http://169.254.169.254/latest/meta-data/placement/availability-zone)"
                    region="${zone::-1}"
                    command -v nfsstat >/dev/null 2>&1 || yum install -y nfs-utils || apt-get -y install nfs-common
                    mkdir -p /mnt/scratch
                    mount -t nfs4 -o nfsvers=4.1 ${zone}./dev/xyz.efs.${region}.amazonaws.com:/ /mnt/scratch
                    chown ubuntu:ubuntu /mnt/scratch
                    chmod 775 /mnt/scratch
                    '''
                    .stripIndent().leftTrim()

    }

    def 'should return the script to create a user and install the key' () {

        given:
        def driver = [:] as AmazonCloudDriver
        when:
        def snippet = driver.scriptCreateUser('pditommaso', 'ssh-rsa AAAAB3NzaC1yc2EAAAABIwAAAQEA...')
        then:
        snippet ==  '''
                    useradd -m -s /bin/bash pditommaso
                    mkdir ~pditommaso/.ssh
                    echo "ssh-rsa AAAAB3NzaC1yc2EAAAABIwAAAQEA..." > ~pditommaso/.ssh/authorized_keys
                    chmod 700 ~pditommaso/.ssh
                    chmod 600 ~pditommaso/.ssh/authorized_keys
                    chown -R pditommaso:pditommaso ~pditommaso/.ssh
                    egrep -i "^wheel:" /etc/group > /dev/null && usermod -aG wheel pditommaso
                    egrep -i "^docker:" /etc/group > /dev/null && usermod -aG docker pditommaso
                    chmod +x /etc/sudoers
                    echo 'pditommaso ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
                    chmod -x /etc/sudoers
                    '''
                    .stripIndent().leftTrim()
    }



    def 'should check the right invocation for boothook script' () {

        given:
        def driver = Spy(AmazonCloudDriver)
        driver.scriptMountEFS(_,_,_) >> { "mount efs" }
        driver.scriptCreateUser(_,_) >> { "create user" }
        driver.scriptMountInstanceStorage(_,_,_) >> { "mount storage" }

        when:
        // no user to create
        // no storage to mount
        // no EFS to mount
        def config1 = [cloud: [
                keyName: 'my-key',
                imageId: 'ami-123',
                instanceType: 'r3.abc',
        ]]

        def result1 = driver.cloudBootHookScript(CloudConfig.create(config1).build())
        then:
        0 * driver.scriptCreateUser(_,_)
        0 * driver.scriptMountEFS(_, _, _)
        0 * driver.scriptMountInstanceStorage(_,_,_)
        result1 == ''


        when:
        // yes user to create
        // no storage to mount
        // no EFS to mount
        def config2 = [cloud: [ userName: 'ubuntu', keyHash: 'ssh-hash 89d98ds090f06da67da', ]]
        driver.cloudBootHookScript(CloudConfig.create(config2).build())
        then:
        1 * driver.scriptCreateUser("ubuntu","ssh-hash 89d98ds090f06da67da")
        0 * driver.scriptMountEFS(_, _, _)
        0 * driver.scriptMountInstanceStorage(_,_,_)


        when:
        // no user to create
        // no instance storage
        // yes EFS to mount
        def config3 = [cloud: [
                keyName:'eurokey',
                userName: 'the-user',
                sharedStorageId: 'fs-213',
                sharedStorageMount: '/mnt/efs']]
        driver.cloudBootHookScript(CloudConfig.create(config3).build())
        then:
        0 * driver.scriptCreateUser(_,_)
        0 * driver.scriptMountInstanceStorage(_,_,_)
        1 * driver.scriptMountEFS('fs-213', '/mnt/efs', 'the-user')


        when:
        // no user to create
        // instance storage to mount
        // no EFS to mount
        def config4 = [cloud: [
                keyName:'eurokey',
                userName: 'the-user',
                instanceType: 'r3.small',
                instanceStorageDevice: '/dev/xyz',
                instanceStorageMount: '/mnt/scratch']]
        driver.cloudBootHookScript(CloudConfig.create(config4).build())
        then:
        0 * driver.scriptCreateUser(_,_)
        1 * driver.scriptMountInstanceStorage('/dev/xyz','/mnt/scratch', 'the-user')
        0 * driver.scriptMountEFS(_, _, _)

    }


    @Unroll
    def 'should get instance type description #type' () {

        given:
        def driver = Spy(AmazonCloudDriver)

        when:
        def bean = driver.describeInstanceType(type)
        then:
        bean.cpus == cpus
        bean.memory.toString() == mem
        bean.disk.toString() == disk
        bean.numOfDisks == num

        where:
        type        | cpus      | mem       | disk      | num
        't2.nano'   | 1         | '512 MB'  | '0'       | 0
        'm3.2xlarge'| 8	        | '30 GB'	| '80 GB'   | 2
        'd2.8xlarge'| 36	    | '244 GB'  | '2 TB'	| 24


    }

    def 'should invoke the right waiters' () {

        given:
        def config = [aws: [accessKey: 'xxx', secretKey: 'yyy', region:"zzz"]]
        def instanceIds = ['i-111','i-222','i-333']
        def client = Mock(AmazonEC2Client)
        def waiters = Mock(AmazonEC2Waiters)
        def driver = new AmazonCloudDriver(client)
        client.waiters() >> waiters

        when:
        driver.waitInstanceStatus(instanceIds, CloudInstanceStatus.STARTED)
        then:
        1 * waiters.instanceRunning() >> { Mock(Waiter) }
        0 * waiters.instanceStatusOk()

        when:
        driver.waitInstanceStatus(instanceIds, CloudInstanceStatus.READY)
        then:
        1 * waiters.instanceRunning() >> { Mock(Waiter) }
        then:
        1 * waiters.instanceStatusOk() >> { Mock(Waiter) }

        when:
        driver.waitInstanceStatus(instanceIds, CloudInstanceStatus.TERMINATED)
        then:
        1 * waiters.instanceTerminated() >> { Mock(Waiter) }

    }


    def 'should make a plain instance request' () {

        given:
        final COUNT = 9
        final AMI = 'ami-123'
        final TYPE = 'mx2.123'
        final KEY = 'my-key'
        final SECURITY =  ['sg-123']
        final SUBNET = 'subnet-666'
        final BLOCK = new BlockDeviceMapping().withDeviceName('/dev/xvdc')
        final IAM_PROFILE = 'foo-role'
        def cfg = Mock(LaunchConfig)
        cfg.getNextflow() >> { new CloudConfig.Nextflow([version:'1.0', trace:'INFO']) }
        def driver = Spy(AmazonCloudDriver)

        when:
        def req = driver.makeRunRequest(COUNT, cfg)

        then:
        (1.._) * cfg.getImageId() >> AMI
        (1.._) * cfg.getInstanceType() >> TYPE
        1 * driver.getUserDataAsBase64(cfg)
        1 * driver.getBlockDeviceMappings(cfg) >> [ BLOCK ]
        (1.._) * cfg.getKeyName() >> KEY
        (1.._) * cfg.getSecurityGroups() >> SECURITY
        (1.._) * cfg.getSubnetId() >> SUBNET
        (1.._) * cfg.getIamProfile() >> IAM_PROFILE

        req instanceof RunInstancesRequest
        req.getMinCount() == COUNT
        req.getMaxCount() == COUNT
        req.getImageId() == AMI
        req.getInstanceType() == TYPE
        req.getKeyName() == KEY
        req.getSecurityGroupIds() == SECURITY
        req.getSubnetId() == SUBNET
        req.getBlockDeviceMappings().size() == 1
        req.getBlockDeviceMappings()[0] == BLOCK
        req.getIamInstanceProfile().getName() == IAM_PROFILE
    }

    def 'should create a spot instance request' () {

        given:
        final COUNT = 33
        final AMI = 'ami-996'
        final TYPE = 'lx2.663'
        final KEY = 'your-key'
        final SECURITY =  ['sg-775']
        final SUBNET = 'subnet-784'
        final PRICE = '1.52'
        final BLOCK = new BlockDeviceMapping().withDeviceName('/dev/xvdc')
        final IAM_PROFILE = 'bar-profile'
        def cfg = Mock(LaunchConfig)
        cfg.getNextflow() >> { new CloudConfig.Nextflow([version:'1.0', trace:'INFO']) }
        def driver = Spy(AmazonCloudDriver)

        when:
        def req = driver.makeSpotRequest(COUNT, cfg)

        then:
        (1.._) * cfg.getImageId() >> AMI
        (1.._) * cfg.getInstanceType() >> TYPE
        1 * driver.getUserDataAsBase64(cfg)
        1 * driver.getBlockDeviceMappings(cfg) >> [ BLOCK ]
        (1.._) * cfg.getKeyName() >> KEY
        (1.._) * cfg.getSecurityGroups() >> SECURITY
        (1.._) * cfg.getSubnetId() >> SUBNET
        (1.._) * cfg.getIamProfile() >> IAM_PROFILE
        1 * cfg.getSpotPrice() >> PRICE

        req instanceof RequestSpotInstancesRequest
        req.getInstanceCount() == COUNT
        req.getSpotPrice() == PRICE
        req.getLaunchSpecification().getImageId() == AMI
        req.getLaunchSpecification().getInstanceType() == TYPE
        req.getLaunchSpecification().getKeyName() == KEY
        req.getLaunchSpecification().getAllSecurityGroups() == [ new GroupIdentifier().withGroupId(SECURITY[0]) ]
        req.getLaunchSpecification().getSubnetId() == SUBNET
        req.getLaunchSpecification().getBlockDeviceMappings().size() == 1
        req.getLaunchSpecification().getBlockDeviceMappings()[0] == BLOCK
        req.getLaunchSpecification().getIamInstanceProfile().name == IAM_PROFILE

    }

    def 'should create block device mapping list' () {
        given:
        final SIZE = MemoryUnit.of('20 GB')
        final DEVICE = '/dev/xyz'
        final AMI = 'ami-3232'
        final SNAPSHOT = 'snap-2121'
        final BLOCK = new BlockDeviceMapping()
                    .withDeviceName('/root')
                    .withEbs( new EbsBlockDevice().withSnapshotId(SNAPSHOT) )

        def cfg = Mock(LaunchConfig)
        cfg.getImageId() >> AMI
        def driver = Spy(AmazonCloudDriver)

        when:
        def maps = driver.getBlockDeviceMappings(cfg)

        then:
        (1.._) * cfg.getInstanceStorageDevice() >> DEVICE
        (1.._) * cfg.getBootStorageSize() >> SIZE
        1 * driver.getRooDeviceMapping(AMI) >> BLOCK

        maps.size() == 2
        maps[0].deviceName == DEVICE
        maps[0].virtualName== 'ephemeral0'
        maps[1].deviceName == '/root'
        maps[1].ebs.snapshotId == SNAPSHOT
        maps[1].ebs.volumeSize == (int)SIZE.toGiga()

    }

    def 'should fetch aws region' () {
        given:
        def driver = Spy(AmazonCloudDriver, constructorArgs: [Mock(AmazonEC2Client)])

        when:
        def result = driver.fetchRegion()
        then:
        driver.getUrl('http://169.254.169.254/latest/meta-data/placement/availability-zone') >> 'eu-west-1a'
        result == 'eu-west-1'
    }

    def 'should fetch aws role' () {
        given:
        def driver = Spy(AmazonCloudDriver, constructorArgs: [Mock(AmazonEC2Client)])

        when:
        def result = driver.fetchIamProfile()
        then:
        driver.getUrl('http://169.254.169.254/latest/meta-data/iam/security-credentials/') >> 'iam-role-here'
        result == 'iam-role-here'
    }

    def 'should validate a config' () {

        given:
        def driver = Spy(AmazonCloudDriver, constructorArgs: [Mock(AmazonEC2Client)])
        def cfg

        when:
        cfg = CloudConfig.create(cloud: [imageId:'ami-123', instanceType: 't2.xxx'])
        driver.validate(cfg)
        then:
        1 * driver.describeInstanceType('t2.xxx') >> { [:] as CloudInstanceType }
        1 * driver.getAccessKey() >> 'xxx'
        1 * driver.getSecretKey() >> 'yyy'

        when:
        cfg = CloudConfig.create(cloud: [imageId:'ami-123', instanceType: 't2.xxx', iamProfile: 'foo-role'])
        driver.validate(cfg)
        then:
        1 * driver.describeInstanceType('t2.xxx') >> { [:] as CloudInstanceType }
        1 * driver.getAccessKey() >> null
        0 * driver.fetchIamProfile()

        when:
        cfg = CloudConfig.create(cloud: [imageId:'ami-123', instanceType: 't2.xxx'])
        driver.validate(cfg)
        then:
        1 * driver.describeInstanceType('t2.xxx') >> { [:] as CloudInstanceType }
        1 * driver.getAccessKey() >> null
        1 * driver.fetchIamProfile() >> 'secret'
        cfg.getIamProfile() == 'secret'

        when:
        cfg = CloudConfig.create(cloud: [imageId:'ami-123', instanceType: 't2.xxx'])
        driver.validate(cfg)
        then:
        1 * driver.describeInstanceType('t2.xxx') >> { [:] as CloudInstanceType }
        1 * driver.getAccessKey() >> null
        1 * driver.fetchIamProfile() >> null
        thrown(IllegalArgumentException)
    }

    def 'should validate config with profile' () {
        given:
        def driver = Spy(AmazonCloudDriver, constructorArgs: [Mock(AmazonEC2Client)])
        def cfg = Mock(LaunchConfig)


    }
}

