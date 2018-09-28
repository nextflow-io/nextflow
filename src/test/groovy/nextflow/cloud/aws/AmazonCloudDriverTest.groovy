/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.cloud.aws
import javax.mail.internet.MimeMessage
import javax.mail.internet.MimeMultipart

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
import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.lang.Unroll
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

    def 'should create mime-multipart cloud-init' () {

        given:
        def session = javax.mail.Session.getInstance(new Properties());
        def SCRIPT_INIT = 'cloud-init-script-text'
        def SCRIPT_BOOT = 'boot-hook-script'
        def config = Mock(LaunchConfig)
        def driver = Spy(AmazonCloudDriver)

        when:
        def payload = driver.getUserDataAsString(config)
        def message = new MimeMessage(session, new ByteArrayInputStream(payload.getBytes()))
        def parts = (MimeMultipart)message.getContent()
        then:
        1 * driver.cloudInitScript(config) >> SCRIPT_INIT
        1 * driver.cloudBootHookScript(config) >> SCRIPT_BOOT
        parts.getCount() == 2
        parts.getBodyPart(0).getContent().text == SCRIPT_INIT
        parts.getBodyPart(1).getContent().text == SCRIPT_BOOT
    }

    def 'should create bash profile properly' () {

        given:
        String bash
        LaunchConfig config
        def client = Mock(AmazonEC2Client)
        def driver = Spy(AmazonCloudDriver, constructorArgs:[client])

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
            export NXF_CLUSTER_JOIN='cloud:aws:cluster-x'
            '''
            .stripIndent() .leftTrim()

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
            export NXF_CLUSTER_JOIN='cloud:aws:cluster-x'
            export NXF_TRACE='io.nextflow.TaskProcess'
            '''
                .stripIndent() .leftTrim()


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
            export NXF_CLUSTER_JOIN='cloud:aws:cluster-x'
            export NXF_OPTS='-XX:this-and-that'
            '''
                .stripIndent() .leftTrim()


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
            export NXF_CLUSTER_JOIN='cloud:aws:cluster-x'
            export NXF_WORK='/mount/my/path/illo/work'
            export NXF_ASSETS='/mount/my/path/illo/projects'
            '''
                .stripIndent() .leftTrim()

        when:
        config = Mock(LaunchConfig)
        config.getNextflow() >> new CloudConfig.Nextflow([version: '0.21.0'])
        config.getClusterName() >> 'cluster-x'
        config.getInstanceStorageMount() >> '/scratch'
        bash = driver.scriptBashEnv(config)
        then:
        driver.hasInstanceStorage0(_) >> true
        bash == '''
            export NXF_VER='0.21.0'
            export NXF_MODE='ignite'
            export NXF_EXECUTOR='ignite'
            export NXF_CLUSTER_JOIN='cloud:aws:cluster-x'
            export NXF_TEMP='/scratch'
            '''
                .stripIndent() .leftTrim()


        when:
        driver.accessKey = 'alpha'
        driver.secretKey = 'beta'
        driver.region = 'region-1'

        config = Mock(LaunchConfig)
        config.getNextflow() >> new CloudConfig.Nextflow([version: '0.21.0'])
        config.getClusterName() >> 'cluster-x'

        bash = driver.scriptBashEnv(config)
        then:
        bash == '''
            export NXF_VER='0.21.0'
            export NXF_MODE='ignite'
            export NXF_EXECUTOR='ignite'
            export NXF_CLUSTER_JOIN='cloud:aws:cluster-x'
            export AWS_ACCESS_KEY_ID='alpha'
            export AWS_SECRET_ACCESS_KEY='beta'
            export AWS_DEFAULT_REGION='region-1'
            '''
                .stripIndent() .leftTrim()
    }



    def 'should fill up the cloud-boot template' () {

        given:
        def CONFIG = [
                aws: [accessKey: 'xxx', secretKey: 'yyy', region: 'eu-west-1'],
                cloud: [
                imageId: 'ami-123',
                instanceType: 'r3.abc',
                securityGroup: 'sg-df72b9ba',
                keyName: 'eurokey',
                subnetId: 'subnet-05222a43',
                bootStorageSize: '10GB',
                instanceStorageMount: '/scratch',
                sharedStorageId: 'fs-1803efd1',
                sharedStorageMount: '/mnt/efs',
                spotPrice: 0.55,
                nextflow: [version: '0.23.0'],
                autoscale: [ enabled: true, spotPrice: 0.33 ]

        ]]

        def driver = Spy(AmazonCloudDriver)
        driver.accessKey = 'xxx'
        driver.secretKey = 'yyy'
        driver.region = 'eu-west-1'

        def config = CloudConfig.create(CONFIG)
                        .setRole('master')
                        .setClusterName('my-cluster')
                        .build()

        when:
        def script = driver.cloudInitScript(config)
        then:
        1 * driver.hasInstanceStorage(config) >> true
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

    def 'should mount instance storage given mount path and device name' () {

        given:
        def script
        def driver = Spy(AmazonCloudDriver)
        def config = Mock(LaunchConfig)

        when: '''
            It is provided `mount path` and `device name`
            it should return a script to format and mount that device
        '''
        script = driver.scriptMountInstanceStorage(config)
        then:
        1 * config.getInstanceStorageMount() >> '/mnt/scratch'
        1 * config.getInstanceStorageDevice() >> '/dev/xyz'
        1 * config.getUserName() >> 'ubuntu'
        0 * driver.instanceStorageDeviceNames(config)
        script ==  """
                    # unmount any ephemeral volume
                    for x in \$(df | grep ephemeral | awk '{print \$1}'); do umount \$x; done
                    # format and mount storage volume
                    mkfs.ext4 -E nodiscard /dev/xyz
                    mkdir -p /mnt/scratch
                    mount -o discard /dev/xyz /mnt/scratch
                    chown -R ubuntu:ubuntu /mnt/scratch
                    chmod 775 /mnt/scratch
                    """
                .stripIndent().leftTrim()

    }


    def 'should mount lvm volume given mount path' () {

        given:
        def script
        def driver = Spy(AmazonCloudDriver)
        def config = Mock(LaunchConfig)

        when: '''
            the instanceType has more than one ephemeral instance storage volumes available
            it creates a LVM volume and mount that volume
        '''
        script = driver.scriptMountInstanceStorage(config)
        then:
        1 * config.getInstanceStorageMount() >> '/scratch'
        1 * config.getInstanceStorageDevice() >> null
        1 * config.getUserName() >> 'foo'
        1 * driver.instanceStorageDeviceNames(config) >> ['/dev/sdx', '/dev/sdz']
        script ==  """
                    # unmount any ephemeral volume
                    for x in \$(df | grep ephemeral | awk '{print \$1}'); do umount \$x; done
                    # install LVM2
                    command -v vgscan >/dev/null 2>&1 || yum install -y lvm2 || apt-get -y install lvm2
                    # create lvm volume
                    pvcreate -y /dev/sdx /dev/sdz
                    vgcreate -y eph /dev/sdx /dev/sdz
                    vgscan
                    lvcreate -y -n data -l 100%FREE eph
                    # format and mount storage volume
                    mkfs.ext4 -E nodiscard /dev/eph/data
                    mkdir -p /scratch
                    mount -o discard /dev/eph/data /scratch
                    chown -R foo:foo /scratch
                    chmod 775 /scratch
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
        0 * driver.scriptCreateUser(_,_) >> null
        0 * driver.scriptMountEFS(_, _, _) >> null
        0 * driver.scriptMountInstanceStorage(_) >> null
        result1 == ''


        when:
        // yes user to create
        // no storage to mount
        // no EFS to mount
        def config2 = [cloud: [ userName: 'ubuntu', keyHash: 'ssh-hash 89d98ds090f06da67da', ]]
        driver.cloudBootHookScript(CloudConfig.create(config2).build())
        then:
        1 * driver.scriptCreateUser("ubuntu","ssh-hash 89d98ds090f06da67da") >> null
        0 * driver.scriptMountEFS(_, _, _) >> null
        0 * driver.scriptMountInstanceStorage(_) >> null


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
        0 * driver.scriptMountInstanceStorage(_)
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
        0 * driver.scriptCreateUser(_,_) >> null
        1 * driver.describeInstanceType('r3.small') >> new CloudInstanceType(numOfDisks: 1)
        1 * driver.scriptMountInstanceStorage(_) >> null
        0 * driver.scriptMountEFS(_, _, _) >> null

    }

    def 'should return a list of device names' () {

        given:
        def driver = Spy(AmazonCloudDriver)
        expect:
        driver.instanceStorageDeviceNames(1) == ['/dev/sdb']
        driver.instanceStorageDeviceNames(2) == ['/dev/sdb', '/dev/sdc']
        driver.instanceStorageDeviceNames(4) == ['/dev/sdb', '/dev/sdc', '/dev/sdd','/dev/sde']

        when:
        driver.instanceStorageDeviceNames(0)
        then:
        thrown(IllegalArgumentException)
    }

    def 'should return list of NVMe device names' () {
        given:
        def driver = Spy(AmazonCloudDriver)
        expect:
        driver.instanceNVMeDeviceNames(1) == ['/dev/nvme0n1']
        driver.instanceNVMeDeviceNames(2) == ['/dev/nvme0n1', '/dev/nvme1n1']
        driver.instanceNVMeDeviceNames(4) == ['/dev/nvme0n1', '/dev/nvme1n1', '/dev/nvme2n1','/dev/nvme3n1']

        when:
        driver.instanceNVMeDeviceNames(0)
        then:
        thrown(IllegalArgumentException)
    }

    def 'should return list of device given configuration' () {

        given:
        def driver = Spy(AmazonCloudDriver)
        def config = Mock(LaunchConfig)
        List result

        when:
        result = driver.instanceStorageDeviceNames(config)
        then:
        1 * config.getInstanceStorageDevice() >> '/dev/sdx'
        result == ['/dev/sdx']

        when:
        result = driver.instanceStorageDeviceNames(config)
        then:
        1 * config.getInstanceStorageDevice() >> null
        1 * config.getInstanceType() >> 'xy.abc'
        1 * driver.describeInstanceType('xy.abc') >> new CloudInstanceType(numOfDisks: 3)
        result == ['/dev/sdb', '/dev/sdc', '/dev/sdd']

        when:
        result = driver.instanceStorageDeviceNames(config)
        then:
        1 * config.getInstanceStorageDevice() >> null
        1 * config.getInstanceType() >> 'f1.xxx'
        1 * driver.describeInstanceType('f1.xxx') >> new CloudInstanceType(numOfDisks: 4)
        result == ['/dev/nvme0n1', '/dev/nvme1n1', '/dev/nvme2n1', '/dev/nvme3n1']
    }

    def 'should check instance storage availability' () {
        given:
        def INSTANCE = 'm1.abc'
        def MOUNT = '/mnt/somewhere'
        def driver = Spy(AmazonCloudDriver)
        def config = Mock(LaunchConfig)
        def result

        when: '''
        the instance storage mount is not defined (ie. null) by definition => return FALSE
        '''
        result = driver.hasInstanceStorage0(config)
        then:
        1 * config.getInstanceStorageMount() >> null
        result == false

        when: '''
        the storage mount is defined BUT the instance has not ephemeral volumes => return FALSE
        '''
        result = driver.hasInstanceStorage0(config)
        then:
        1 * config.getInstanceType() >> INSTANCE
        1 * config.getInstanceStorageMount() >> MOUNT
        1 * driver.describeInstanceType(INSTANCE) >> new CloudInstanceType()
        result == false

        when: '''
        the mount path is defined and at least an ephemeral volume is provided by the instance ==> TRUE
        '''
        result = driver.hasInstanceStorage0(config)
        then:
        1 * config.getInstanceType() >> INSTANCE
        1 * config.getInstanceStorageMount() >> MOUNT
        1 * driver.describeInstanceType(INSTANCE) >> new CloudInstanceType(numOfDisks: 1)
        result == true

    }

    def 'should check NVMe storage' () {
        given:
        def driver = new AmazonCloudDriver()
        expect:
        driver.isNVMe('f1.xxx')
        driver.isNVMe('i3.xxx')
        !driver.isNVMe('any.other')
    }


    @Unroll
    def 'should get instance type description #type' () {

        given:
        def TYPE = 'm3.2xlarge'
        def CPUS = 8
        def MEM = MemoryUnit.of('30 GB')
        def DISK = MemoryUnit.of('80 GB')
        def NUM = 2
        def REGION = 'eu-west-1'

        def reader = Mock(AmazonPriceReader)
        def driver = Spy(AmazonCloudDriver)
        driver.region = REGION
        when:
        def bean = driver.describeInstanceType(TYPE)
        then:
        1 * driver.createReader(REGION) >> reader
        1 * reader.getInstanceType(TYPE) >> new CloudInstanceType(cpus: CPUS, memory: MEM, disk: DISK, numOfDisks: NUM)
        bean.cpus == CPUS
        bean.memory == MEM
        bean.disk == DISK
        bean.numOfDisks == NUM
    }

    def 'should invoke the right waiters' () {

        given:
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
        (1.._) * cfg.getInstanceRole() >> IAM_PROFILE

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
        (1.._) * cfg.getInstanceRole() >> IAM_PROFILE
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

    def 'should create an empty block device mapping' () {

        def config = Mock(LaunchConfig)
        def driver = Spy(AmazonCloudDriver)
        List<BlockDeviceMapping> mappings

        when:
        mappings = driver.getBlockDeviceMappings(config)
        then:
        1 * config.getBootStorageSize() >> null
        1 * driver.hasInstanceStorage(config) >> false
        mappings == null

    }

    def 'should create boot storage mapping' () {
        given:
        final SIZE = MemoryUnit.of('20 GB')
        final AMI = 'ami-3232'
        final SNAPSHOT = 'snap-2121'
        final BLOCK = new BlockDeviceMapping()
                    .withDeviceName('/root')
                    .withEbs( new EbsBlockDevice().withSnapshotId(SNAPSHOT) )

        def config = Mock(LaunchConfig)
        def driver = Spy(AmazonCloudDriver)
        List<BlockDeviceMapping> mappings

        when:
        mappings = driver.getBlockDeviceMappings(config)
        then:
        1 * driver.hasInstanceStorage(config) >> false
        1 * config.getBootStorageSize() >> SIZE
        1 * config.getImageId() >> AMI
        1 * driver.getRooDeviceMapping(AMI) >> BLOCK
        mappings.size() == 1
        mappings[0].deviceName == BLOCK.deviceName
        mappings[0].ebs.snapshotId == SNAPSHOT
        mappings[0].ebs.volumeSize == (int)SIZE.toGiga()

    }

    def 'should create instance storage mappings' () {
        given:
        def config = Mock(LaunchConfig)
        def driver = Spy(AmazonCloudDriver)
        List<BlockDeviceMapping> mappings

        when:
        mappings = driver.getBlockDeviceMappings(config)
        then:
        1 * config.getBootStorageSize() >> null
        1 * config.getInstanceType() >> 'r3.xxx'
        1 * driver.hasInstanceStorage(config) >> true
        1 * driver.instanceStorageDeviceNames(config) >> ['/dev/aa', '/dev/bb']
        then:
        mappings.size() ==2
        mappings[0].deviceName == '/dev/aa'
        mappings[0].virtualName == 'ephemeral0'
        mappings[1].deviceName == '/dev/bb'
        mappings[1].virtualName == 'ephemeral1'
    }

    def 'should not create block mapping for NVMe storage' () {

        given:
        def config = Mock(LaunchConfig)
        def driver = Spy(AmazonCloudDriver)
        List<BlockDeviceMapping> mappings

        when:
        mappings = driver.getBlockDeviceMappings(config)
        then:
        1 * config.getBootStorageSize() >> null
        1 * driver.hasInstanceStorage(config) >> true
        1 * driver.isNVMe(_) >> true
        then:
        mappings == null
    }

    def 'should fetch aws region' () {
        given:
        def driver = Spy(AmazonCloudDriver, constructorArgs: [Mock(AmazonEC2Client)])

        when:
        def result = driver.fetchRegion()
        then:
        1 * driver.getUrl('http://169.254.169.254/latest/meta-data/placement/availability-zone') >> 'eu-west-1a'
        result == 'eu-west-1'
    }

    def 'should fetch aws role' () {
        given:
        def driver = Spy(AmazonCloudDriver, constructorArgs: [Mock(AmazonEC2Client)])

        when:
        def result = driver.fetchIamRole()
        then:
        1 * driver.getUrl('http://169.254.169.254/latest/meta-data/iam/security-credentials/') >> 'iam-role-here'
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
        cfg = CloudConfig.create(cloud: [imageId:'ami-123', instanceType: 't2.xxx', instanceRole: 'foo-role'])
        driver.validate(cfg)
        then:
        1 * driver.describeInstanceType('t2.xxx') >> { [:] as CloudInstanceType }
        1 * driver.getAccessKey() >> null
        0 * driver.fetchIamRole()

        when:
        cfg = CloudConfig.create(cloud: [imageId:'ami-123', instanceType: 't2.xxx'])
        driver.validate(cfg)
        then:
        1 * driver.describeInstanceType('t2.xxx') >> { [:] as CloudInstanceType }
        1 * driver.getAccessKey() >> null
        1 * driver.fetchIamRole() >> 'secret'
        cfg.getInstanceRole() == 'secret'

        when:
        cfg = CloudConfig.create(cloud: [imageId:'ami-123', instanceType: 't2.xxx'])
        driver.validate(cfg)
        then:
        1 * driver.describeInstanceType('t2.xxx') >> { [:] as CloudInstanceType }
        1 * driver.getAccessKey() >> null
        1 * driver.fetchIamRole() >> null
        thrown(IllegalArgumentException)
    }

    def 'should validate config with profile' () {
        given:
        def driver = Spy(AmazonCloudDriver, constructorArgs: [Mock(AmazonEC2Client)])
        def cfg = Mock(LaunchConfig)


    }
}

