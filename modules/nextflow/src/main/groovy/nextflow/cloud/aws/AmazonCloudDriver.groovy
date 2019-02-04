/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import nextflow.cloud.CloudScripts

import static nextflow.cloud.CloudConst.TAG_CLUSTER_NAME
import static nextflow.cloud.CloudConst.TAG_CLUSTER_ROLE

import javax.mail.Session
import javax.mail.internet.MimeBodyPart
import javax.mail.internet.MimeMessage
import javax.mail.internet.MimeMultipart

import com.amazonaws.auth.BasicAWSCredentials
import com.amazonaws.regions.Region
import com.amazonaws.regions.RegionUtils
import com.amazonaws.services.batch.AWSBatchClient
import com.amazonaws.services.ec2.AmazonEC2Client
import com.amazonaws.services.ec2.model.BlockDeviceMapping
import com.amazonaws.services.ec2.model.CreateTagsRequest
import com.amazonaws.services.ec2.model.DescribeImagesRequest
import com.amazonaws.services.ec2.model.DescribeInstanceStatusRequest
import com.amazonaws.services.ec2.model.DescribeInstancesRequest
import com.amazonaws.services.ec2.model.DescribeSpotInstanceRequestsRequest
import com.amazonaws.services.ec2.model.DescribeSpotPriceHistoryRequest
import com.amazonaws.services.ec2.model.Filter
import com.amazonaws.services.ec2.model.GroupIdentifier
import com.amazonaws.services.ec2.model.IamInstanceProfileSpecification
import com.amazonaws.services.ec2.model.Image
import com.amazonaws.services.ec2.model.Instance
import com.amazonaws.services.ec2.model.LaunchSpecification
import com.amazonaws.services.ec2.model.RequestSpotInstancesRequest
import com.amazonaws.services.ec2.model.RequestSpotInstancesResult
import com.amazonaws.services.ec2.model.RunInstancesRequest
import com.amazonaws.services.ec2.model.RunInstancesResult
import com.amazonaws.services.ec2.model.SpotPrice
import com.amazonaws.services.ec2.model.Tag
import com.amazonaws.services.ec2.model.TerminateInstancesRequest
import com.amazonaws.util.Base64
import com.amazonaws.waiters.WaiterParameters
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.SimpleType
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.cloud.CloudDriver
import nextflow.cloud.LaunchConfig
import nextflow.cloud.types.CloudInstance
import nextflow.cloud.types.CloudInstanceStatus
import nextflow.cloud.types.CloudInstanceType
import nextflow.cloud.types.CloudSpotPrice
import nextflow.exception.AbortOperationException
import nextflow.processor.TaskTemplateEngine
import nextflow.util.ServiceName
/**
 * Implements the cloud driver of Amazon AWS
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@ServiceName('aws')
class AmazonCloudDriver implements CloudDriver {

    /**
     * Reference to {@link AmazonEC2Client} object
     */
    private AmazonEC2Client ec2Client

    /**
     * The AWS access key credentials (optional)
     */
    private String accessKey

    /**
     * The AWS secret key credentials (optional)
     */
    private String secretKey

    /**
     * The AWS region eg. {@code eu-west-1}. If it's not specified the current region is retrieved from
     * the EC2 instance metadata
     */
    private String region

    /**
     * Initialise the Amazon cloud driver with default (empty) parameters
     */
    AmazonCloudDriver() {
        this(Collections.emptyMap())
    }

    /**
     * Initialise the Amazon cloud driver with the specified parameters
     *
     * @param config
     *      A map holding the driver parameters:
     *      - accessKey: the access key credentials
     *      - secretKey: the secret key credentials
     *      - region: the AWS region
     */
    @CompileDynamic
    AmazonCloudDriver(Map config) {
        // -- get the aws credentials
        List credentials
        if( config.accessKey && config.secretKey ) {
            this.accessKey = config.accessKey
            this.secretKey = config.secretKey
        }
        else if( (credentials=Global.getAwsCredentials()) ) {
            this.accessKey = credentials[0]
            this.secretKey = credentials[1]
        }

        if( !accessKey && !fetchIamRole() )
            throw new AbortOperationException("Missing AWS security credentials -- Provide access/security keys pair or define a IAM instance profile (suggested)")

        // -- get the aws default region
        region = config.region ?: Global.getAwsRegion() ?: fetchRegion()
        if( !region )
            throw new AbortOperationException('Missing AWS region -- Make sure to define in your system environment the variable `AWS_DEFAULT_REGION`')

    }

    /**
     * only for testing purpose -- do not use
     */
    @CompileDynamic
    protected AmazonCloudDriver( AmazonEC2Client client ) {
        this.ec2Client = client
    }

    /**
     * @return The current set AWS access key
     */
    protected String getAccessKey() { accessKey }

    /**
     * @return The current set AWS secret key
     */
    protected String getSecretKey() { secretKey }

    /**
     * Retrieve the current IAM role eventually define for a EC2 instance.
     * See http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/iam-roles-for-amazon-ec2.html#instance-metadata-security-credentials
     *
     * @return
     *      The IAM role name associated to this instance or {@code null} if no role is defined or
     *      it's not a EC2 instance
     */
    protected String fetchIamRole() {
        try {
            def role = getUrl('http://169.254.169.254/latest/meta-data/iam/security-credentials/').readLines()
            if( role.size() != 1 )
                throw new IllegalArgumentException("Not a valid EC2 IAM role")
            return role.get(0)
        }
        catch( IOException e ) {
            log.trace "Unable to fetch IAM credentials -- Cause: ${e.message}"
            return null
        }
    }

    /**
     * Retrieve the AWS region from the EC2 instance metadata.
     * See http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-instance-metadata.html
     *
     * @return
     *      The AWS region of the current EC2 instance eg. {@code eu-west-1} or
     *      {@code null} if it's not an EC2 instance.
     */
    protected String fetchRegion() {
        try {
            def zone = getUrl('http://169.254.169.254/latest/meta-data/placement/availability-zone')
            zone ? zone.substring(0,zone.length()-1) : null
        }
        catch (IOException e) {
            log.debug "Cannot fetch AWS region", e
            return null
        }
    }

    /**
     * Helper method to map a region string to a {@link Region} object.
     *
     * @param region An AWS region string identifier eg. {@code eu-west-1}
     * @return A {@link Region} corresponding to the specified region string
     */
    private Region getRegionObj(String region) {
        final result = RegionUtils.getRegion(region)
        if( !result )
            throw new IllegalArgumentException("Not a valid AWS region name: $region");
        return result
    }

    /**
     * Gets or lazily creates an {@link AmazonEC2Client} instance given the current
     * configuration parameter
     *
     * @return
     *      An {@link AmazonEC2Client} instance
     */
    synchronized AmazonEC2Client getEc2Client() {

        if( ec2Client )
            return ec2Client

        def result = (accessKey && secretKey
                ? new AmazonEC2Client(new BasicAWSCredentials(accessKey, secretKey))
                : new AmazonEC2Client())

        if( region )
            result.setRegion(getRegionObj(region))

        return result
    }

    /**
     * Gets or lazily creates an {@link AWSBatchClient} instance given the current
     * configuration parameter
     *
     * @return
     *      An {@link AWSBatchClient} instance
     */
    @Memoized
    AWSBatchClient getBatchClient() {
        def result = (accessKey && secretKey
                ? new AWSBatchClient(new BasicAWSCredentials(accessKey, secretKey))
                : new AWSBatchClient())

        if( region )
            result.setRegion(getRegionObj(region))

        return result
    }

    /**
     * Checks if the instance type specified in the configuration provide one or more ephemeral instance storage.
     * See http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/InstanceStorage.html
     *
     * NOTE: this method is declared as memoized to avoid to warn multiple times for the same warning message.
     *
     * @param cfg
     *      A {@link LaunchConfig} object representing the current instance configuration.
     * @return
     *      {@code true} whenever the instance type provides one more instance storage volumes or
     *      {@code false} otherwise.
     */
    @Memoized
    @PackageScope
    boolean hasInstanceStorage(LaunchConfig cfg) {
        hasInstanceStorage0(cfg)
    }

    @PackageScope
    boolean hasInstanceStorage0(LaunchConfig cfg) {
        final mount = cfg.instanceStorageMount
        if( !mount )
            return false

        final type = cfg.instanceType
        final hasDisks = describeInstanceType(type).numOfDisks > 0
        if( !hasDisks ) {
            log.warn "EC2 instance type: $type does not provide any ephemeral storage -- skipping instance storage mount: $mount"
        }
        return hasDisks
    }

    /**
     * A BASH snippet used to initialise the EC2 instance user account.
     *
     * @param userName
     *          The Linux user name to be configured in the EC2 instance.
     * @param key
     *          The SSH public key to be configured in the EC2 instance.
     * @return
     *          A BASH script setting up the user account for the given parameters.
     */
    @PackageScope
    String scriptCreateUser(String userName, String key) {
        CloudScripts.scriptCreateUser(userName,key)
    }

    /**
     * A BASH snippet used to initialise the Elastic File System storage.
     * See https://aws.amazon.com/efs/
     *
     * @param fileSystemId The EFS storage unique ID
     * @param fileSystemMount The EFS posix mount path eg. {@code /mnt/data}
     * @param userName The Linux user which owns the file system permissions
     * @return A BASH script setting up the EFS storage for the given parameters
     */
    @PackageScope
    String scriptMountEFS(String fileSystemId, String fileSystemMount, String userName) {

        """\
        zone="\$(curl -s http://169.254.169.254/latest/meta-data/placement/availability-zone)"
        region="\${zone::-1}"
        command -v nfsstat >/dev/null 2>&1 || yum install -y nfs-utils || apt-get -y install nfs-common
        mkdir -p $fileSystemMount
        mount -t nfs4 -o nfsvers=4.1 \${zone}.${fileSystemId}.efs.\${region}.amazonaws.com:/ $fileSystemMount
        chown ${userName}:${userName} $fileSystemMount
        chmod 775 $fileSystemMount
        """
        .stripIndent()
    }

    /**
     * @return A BASH snippet used to un-mount any ephemeral instance storage volume
     */
    private String scriptUnmountEphemeralVolume() {
        """\
        # unmount any ephemeral volume
        for x in \$(df | grep ephemeral | awk '{print \$1}'); do umount \$x; done
        """
        .stripIndent()
    }

    /**
     * A BASH snippet that creates a LVM volume for the given device names
     *
     * @param devicePaths A list of one or more storage device names eg. {@code /dev/sdb}, {@code /dev/sdc}, etc
     * @param volumeGroup The LVM group name
     * @param volumeName The LVM volume name
     * @return A BASH script setting up a LVM volume for the given parameters
     */
    private String scriptCreateLogicalVolume(List<String> devicePaths, String volumeGroup, String volumeName) {
        assert devicePaths, "Device names list should contain at least one entry"

        """\
        # install LVM2
        command -v vgscan >/dev/null 2>&1 || yum install -y lvm2 || apt-get -y install lvm2
        # create lvm volume
        pvcreate -y ${devicePaths.join(' ')}
        vgcreate -y $volumeGroup ${devicePaths.join(' ')}
        vgscan
        lvcreate -y -n $volumeName -l 100%FREE $volumeGroup
        """
        .stripIndent()
    }

    /**
     * A BASH snippet formatting and mounting the specified device storage
     *
     * @param devicePath The posix device name of the storage to mount
     * @param mountPath The posix file path where mount the specified device
     * @param user The posix user name which own the file system permissions
     * @return A BASH script setting up the volume for the specified parameters
     */
    private String scriptMountStorageVolume(String devicePath, String mountPath, String user) {
        """\
        # format and mount storage volume
        mkfs.ext4 -E nodiscard $devicePath
        mkdir -p $mountPath
        mount -o discard $devicePath $mountPath
        chown -R $user:$user $mountPath
        chmod 775 $mountPath
        """
        .stripIndent()
    }

    /**
     * A BASH snippet formatting and mounting the one or more ephemeral instance storage volumes
     * given the current instance configuration.
     *
     * When the configuration object specifies the instance storage device name and mount path,
     * the volume is formatted and mounted accordingly.
     *
     * When the configuration object only specifies the mount path, all instance storage volumes
     * are grouped together as a sole logical volume and mounted in the specified path.
     *
     * See http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/add-instance-store-volumes.html
     *
     * @param cfg A {@link LaunchConfig} object holding the current user provided configuration
     * @return A BASH script setting up the one or more instance storage volumes given the current
     */
    @PackageScope
    String scriptMountInstanceStorage(LaunchConfig cfg) {
        final mountPath = cfg.getInstanceStorageMount()
        final devicePath = cfg.getInstanceStorageDevice()
        final user = cfg.getUserName()
        assert mountPath,   'Instance storage mount path parameter cannot be empty'
        assert user,        'Instance storage mount user parameter cannot be empty'

        def result = scriptUnmountEphemeralVolume()
        if( devicePath ) {
            result += scriptMountStorageVolume(devicePath, mountPath, user)
            return result
        }

        def names = instanceStorageDeviceNames(cfg)
        if( !names ) throw new IllegalStateException('Instance storage device names cannot be empty')

        result += scriptCreateLogicalVolume(names, 'eph', 'data')
        result += scriptMountStorageVolume('/dev/eph/data', mountPath, user)
        return result
    }

    /**
     * Retrieve the list of device names for the given number of volumes
     * accordingly AWS specification.
     *
     * See http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/device_naming.html#available-ec2-device-names
     *
     * @param num The number of devices, it must be greater that zero
     * @return A list of device names eg {@code ['/dev/sdb', '/dev/sdc', ..]}
     */
    @PackageScope
    List<String> instanceStorageDeviceNames(int num) {
        if( num<=0 )
            throw new IllegalArgumentException("Argument numOfDisks must be greater than zero")
        // instance storage device names always starts from `/dev/sdb
        // see http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/device_naming.html#device-name-limits`
        final int BASE = 'b' as char
        def names = []
        for( int i=0; i<num; i++ ) {
            names << "/dev/sd${(BASE+i) as char}"
        }

        return names
    }

    /**
     * Retrieve the list of device names for the given number of NVMe volumes
     * accordingly AWS specification.
     *
     * See http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/device_naming.html#available-ec2-device-names
     *
     * @param num The number of devices, it must be greater that zero
     * @return A list of device names eg {@code ['/dev/nvme0n1', '/dev/nvme1n1', ..]}
     */
    @PackageScope
    List<String> instanceNVMeDeviceNames(int num) {
        if( num<=0 )
            throw new IllegalArgumentException("Argument numOfDisks must be greater than zero")

        // see http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/device_naming.html#available-ec2-device-names
        def names = []
        for( int i=0; i<num; i++ ) {
            names << "/dev/nvme${i}n1"
        }

        return names
    }

    /**
     * Retrieve the list of ephemeral instance storage device names for the
     * instance type in the given launch configuration.
     *
     * http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/device_naming.html#available-ec2-device-names
     *
     */
    @PackageScope
    List<String> instanceStorageDeviceNames(LaunchConfig cfg) {
        final device = cfg.instanceStorageDevice
        if( device )
            return [device]

        final type = cfg.instanceType
        final instance = describeInstanceType(type)
        final num = instance.numOfDisks
        if( !num )
            return Collections.<String>emptyList()

        /*
         * certain instance family types have NVMe storage
         * see
         *   http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/InstanceStorage.html#instance-store-volumes
         *   http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/device_naming.html#available-ec2-device-names
         */
        isNVMe(type) ? instanceNVMeDeviceNames(num) : instanceStorageDeviceNames(num)
    }

    /**
     * Check if the specified instance type provide NVMe storage type
     *
     * @param type
     *      A instance type identifier eg. {@code c3.large}
     * @return
     *      {@code true} when the specified instance type provides NVMe storage type or
     *      {@code false} otherwise
     */
    @PackageScope
    boolean isNVMe(String type) {
        type.startsWith('i3.') || type.startsWith('f1.')
    }

    /**
     * The BASH environment to be defined in the target EC2 instance
     *
     * @param cfg
     *      The {@link LaunchConfig} object representing the instance configuration
     * @return
     *      A BASH snippet representing the environment to be defined in the target EC2 instance
     */
    @PackageScope
    String scriptBashEnv( LaunchConfig cfg ) {
        def profile = """\
        export NXF_VER='${cfg.nextflow.version}'
        export NXF_MODE='${cfg.nextflow.mode}'
        export NXF_EXECUTOR='ignite'
        export NXF_CLUSTER_JOIN='cloud:aws:${cfg.clusterName}'
        """
        .stripIndent()

        if( cfg.nextflow.trace )
            profile += "export NXF_TRACE='${cfg.nextflow.trace}'\n"

        if( cfg.nextflow.options )
            profile += "export NXF_OPTS='${cfg.nextflow.options}'\n"

        if( cfg.sharedStorageId && cfg.sharedStorageMount ) {
            profile += "export NXF_WORK='${cfg.sharedStorageMount}/${cfg.userName}/work'\n"
            profile += "export NXF_ASSETS='${cfg.sharedStorageMount}/${cfg.userName}/projects'\n"
        }

        if( hasInstanceStorage(cfg) )
            profile += "export NXF_TEMP='${cfg.instanceStorageMount}'\n"

        // access/secret keys are propagated only if IAM profile is not specified
        if( !cfg.instanceRole && accessKey && secretKey ) {
            profile += "export AWS_ACCESS_KEY_ID='$accessKey'\n"
            profile += "export AWS_SECRET_ACCESS_KEY='$secretKey'\n"
        }

        if( region )
            profile += "export AWS_DEFAULT_REGION='$region'\n"

        return profile
    }

    /**
     * A BASH script that download and initialise the nextflow daemon in the target EC2 instance.
     *
     * The script is injected by using the AWS user-data parameter and executed by the Cloud-init
     * mechanism as a shell script.
     *
     * When using Amazon Linux, this script is saved in the following path
     *  /var/lib/cloud/instance/scripts/part-001
     *
     * See
     *  http://cloudinit.readthedocs.io/en/latest/topics/format.html#user-data-script
     *  http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/amazon-linux-ami-basics.html#amazon-linux-cloud-init
     *  http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-instance-metadata.html#instancedata-add-user-data
     *
     * @param cfg
     *      The {@link LaunchConfig} object representing the instance configuration
     * @return
     *      A BASH snippet representing the nextflow initialisation in the target environment
     */
    @PackageScope
    String cloudInitScript(LaunchConfig cfg) {
        // load init script template
        def template =  this.class.getResourceAsStream('cloud-boot.txt')
        if( !template )
            throw new IllegalStateException("Missing `cloud-boot.txt` template resource")

        def binding = new CloudBootTemplateBinding(cfg)
        binding.nextflowConfig = cfg.renderCloudConfigObject()
        binding.bashProfile = scriptBashEnv(cfg)

        new TaskTemplateEngine()
                .setPlaceholder('!' as char)
                .setEnableShortNotation(false)
                .createTemplate(new InputStreamReader(template))
                .make(binding)
                .toString()

    }

    /**
     * A BASH script that initialises the Linux user profile and storage settings.
     *
     * The script is injected by using the AWS user-data parameter and executed by the Cloud-init
     * mechanism as a `boothook`.
     *
     * When using Amazon Linux, this script is saved in the following path
     *  /var/lib/cloud/instance/boothooks/part-002
     *
     * See
     *  http://cloudinit.readthedocs.io/en/latest/topics/format.html#cloud-boothook
     *  http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/amazon-linux-ami-basics.html#amazon-linux-cloud-init
     *  http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-instance-metadata.html#instancedata-add-user-data
     *
     * @param cfg
     *      The {@link LaunchConfig} object representing the instance configuration
     * @return
     *      A BASH snippet representing the EC2 instance initialisation in the target environment
     */
    @PackageScope
    String cloudBootHookScript( LaunchConfig cfg ) {
        def builder = []

        if( cfg.createUser ) {
            builder << scriptCreateUser(cfg.userName, cfg.keyHash)
        }

        if( cfg.sharedStorageId && cfg.sharedStorageMount ) {
            builder << scriptMountEFS(cfg.sharedStorageId, cfg.sharedStorageMount, cfg.userName)
        }

        if( hasInstanceStorage(cfg) ) {
            builder << scriptMountInstanceStorage(cfg)
        }

        if( builder ) {
            builder.add(0, '#!/bin/bash')
        }

        // note: `findAll` remove all empty strings
        builder.join('\n')
    }

    /**
     * Creates the user-data script that bootstrap and initialise the EC2 instance.
     *
     * The user data is encoded a MIME multipart message as requested by the cloud-init
     * mechanism. It is composed by two parts:
     * - the environment initialisation encoded as `text/x-shellscript`
     * - the instance initialisation encoded as  `text/cloud-boothook`
     *
     *
     * When executed in a Amazon Linux instance the relevant files are the following:
     * - configuration:
     *      /etc/cloud/cloud.cfg.d/00_defaults.cfg
     *      /etc/cloud/cloud.cfg
     * - logs:
     *      /var/log/cloud-init.log
     *      /var/log/cloud-init-output.log
     * - scripts:
     *      /var/lib/cloud/
     *      /var/lib/cloud/instance/scripts/part-001
     *      /var/lib/cloud/instance/boothooks/part-002
     *
     * See
     *  http://cloudinit.readthedocs.io/en/latest/topics/format.html#user-data-script
     *  http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/amazon-linux-ami-basics.html#amazon-linux-cloud-init
     *  http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-instance-metadata.html#instancedata-add-user-data
     *  http://cloudinit.readthedocs.io/en/latest/topics/format.html#mime-multi-part-archive
     *  https://forums.aws.amazon.com/thread.jspa?messageID=730555&#730555
     *
     * @param cfg
     *      The {@link LaunchConfig} object representing the instance configuration
     * @return
     *      A BASH script used to initilise and setup the target environment
     */
    @PackageScope
    byte[] getUserData0( LaunchConfig cfg ) {

        MimeMultipart multipart = new MimeMultipart()

        //
        // the first part contains the cloud-init script
        //
        def script1 = new MimeBodyPart()
        script1.setContent( cloudInitScript(cfg), 'text/x-shellscript' )
        multipart.addBodyPart(script1)

        //
        // Second part mount the instance and EFS storage
        // When using Amazon Linux this script is uploaded
        // at this path: /var/lib/cloud/instance/boothooks/part-002
        //
        def text = cloudBootHookScript(cfg)
        if( text ) {
            def script2 = new MimeBodyPart()
            script2.setContent(text, 'text/cloud-boothook')
            multipart.addBodyPart(script2)
        }

        //
        // put both in a mime-message and return a string
        //
        def message = new MimeMessage(Session.getInstance(new Properties()))
        message.setContent(multipart)

        def buffer = new ByteArrayOutputStream()
        message.writeTo(buffer)
        return buffer.toByteArray()
    }


    /**
     * @param cfg
     *      The {@link LaunchConfig} object representing the current instance configuration
     * @return
     *      The user-data content as a MIME-Multipart string
     */
    @PackageScope
    String getUserDataAsString(LaunchConfig cfg) {
        new String( getUserData0(cfg) )
    }

    /**
     * @param cfg
     *      The {@link LaunchConfig} object representing the current instance configuration
     * @return
     *      The user-data content as a MIME-Multipart string encoded as Base64
     */
    @PackageScope
    String getUserDataAsBase64( LaunchConfig cfg ) {
        Base64.encodeAsString(getUserData0(cfg))
    }

    /**
     * A list of {@link BlockDeviceMapping} to mount the required storage
     * volumes in the target instance.
     *
     * @param cfg
     *      The {@link LaunchConfig} object representing the current instance configuration
     * @return
     *      The list of required {@link BlockDeviceMapping}
     *      or {@code null} if no block device mapping is required.
     */
    @PackageScope
    List<BlockDeviceMapping> getBlockDeviceMappings(LaunchConfig cfg) {
        def result = new ArrayList<BlockDeviceMapping>()

        // -- map the ephemeral storage
        if( hasInstanceStorage(cfg) && !isNVMe(cfg.instanceType) ) {
            final deviceNames = instanceStorageDeviceNames(cfg)
            for( int i=0; i<deviceNames.size(); i++ ) {
                def mapping = new BlockDeviceMapping()
                mapping.deviceName = deviceNames[i]
                mapping.virtualName = "ephemeral$i"
                result << mapping
            }
        }

        // -- resize the boot storage
        final size = cfg.bootStorageSize
        if( size ) {
            def rootDevice = getRooDeviceMapping( cfg.imageId )
            rootDevice.ebs.volumeSize = (int)size.toGiga()
            rootDevice.ebs.encrypted = null // <-- note: reset to null the encrypted flag
            result << rootDevice
        }

        log.trace "AWS Block device mapping: $result"

        return result ?: null
    }

    /**
     * Retrieve the {@link BlockDeviceMapping} for the specified AMI
     *
     * @param imageId
     *      An virtual machine AMI identifier
     * @return
     *      The {@link BlockDeviceMapping} for the specified AMI
     */
    @PackageScope
    BlockDeviceMapping getRooDeviceMapping( String imageId ) {
        final req = new DescribeImagesRequest()
        req.setImageIds( [imageId] )
        final Image image = getEc2Client().describeImages(req).getImages().get(0)
        return image.blockDeviceMappings.find { it.deviceName == image.rootDeviceName }
    }

    /**
     * Creates an AWS {@link RunInstancesRequest} object for the given number of instances
     *
     * @param instanceCount
     *      The number of instance to launch
     * @param cfg
     *      The {@link LaunchConfig} object representing the requested instance configuration
     * @return
     *      The {@code RunInstancesRequest} object request to be submitted.
     */
    @PackageScope
    RunInstancesRequest makeRunRequest( int instanceCount, LaunchConfig cfg ) {
        def req = new RunInstancesRequest()
        req.minCount = instanceCount
        req.maxCount = instanceCount
        req.imageId = cfg.imageId
        req.instanceType = cfg.instanceType
        req.userData = getUserDataAsBase64(cfg)
        req.setBlockDeviceMappings( getBlockDeviceMappings(cfg) )
        if( cfg.instanceRole ) {
            def role = new IamInstanceProfileSpecification().withName(cfg.instanceRole)
            req.setIamInstanceProfile(role)
        }

        if( cfg.keyName )
            req.keyName = cfg.keyName

        if( cfg.securityGroups )
            req.securityGroupIds = cfg.securityGroups

        if( cfg.subnetId )
            req.subnetId = cfg.subnetId

        return req
    }

    /**
     * Creates an AWS {@link RequestSpotInstancesRequest} object for the given number of *spot* instances
     *
     * @param instanceCount
     *      The number of instance to launch
     * @param cfg
     *      The {@link LaunchConfig} object representing the requested instance configuration
     * @return
     *      The {@link RequestSpotInstancesRequest} object request to be submitted
     */
    @PackageScope
    RequestSpotInstancesRequest makeSpotRequest( int instanceCount, LaunchConfig cfg ) {

        if( !cfg.imageId ) throw new AbortOperationException("Missing `imageId` setting in cloud configuration")
        if( !cfg.instanceType ) throw new AbortOperationException('Missing `instanceType` setting in cloud configuration')

        def spec = new LaunchSpecification()
        cfg.with {
            spec.setImageId(imageId)
            spec.setInstanceType(instanceType)
            spec.setUserData( getUserDataAsBase64(cfg) )
            spec.setBlockDeviceMappings( getBlockDeviceMappings(cfg) )

            if( keyName )
                spec.keyName = keyName

            if( subnetId )
                spec.subnetId = subnetId

            if( securityGroups ) {
                def allGroups = securityGroups.collect { new GroupIdentifier().withGroupId(it) }
                spec.withAllSecurityGroups(allGroups)
            }

            if( instanceRole ) {
                def role = new IamInstanceProfileSpecification().withName(cfg.instanceRole)
                spec.setIamInstanceProfile(role)
            }
        }

        new RequestSpotInstancesRequest()
                .withInstanceCount(instanceCount)
                .withSpotPrice( cfg.getSpotPrice() )
                .withLaunchSpecification(spec)
    }

    @Override
    List<String> launchInstances(int instanceCount, LaunchConfig cfg) {

        def instanceIds = ( cfg.isSpot()
                ? launchSpotInstances0(instanceCount, cfg)
                : launchPlainInstances0(instanceCount, cfg))

        log.debug "AWS Ec2 instanceIds: $instanceIds"
        return instanceIds
    }

    private List<String> launchPlainInstances0(int instanceCount, LaunchConfig cfg) {
        log.debug "AWS submitting launch request for $instanceCount instances; config: $cfg"
        def req = makeRunRequest(instanceCount, cfg)
        RunInstancesResult response = getEc2Client().runInstances(req)
        return response.reservation.getInstances() *. instanceId
    }

    private List<String> launchSpotInstances0(int instanceCount, LaunchConfig cfg) {
        log.debug "AWS submitting launch request for $instanceCount SPOT instances; config: $cfg"

        def req = makeSpotRequest(instanceCount, cfg)
        RequestSpotInstancesResult response = getEc2Client().requestSpotInstances(req)
        def spotIds = response.spotInstanceRequests *. spotInstanceRequestId
        log.debug "AWS spot request IDs: ${spotIds.join(',')}"

        // -- wait for fulfill the spot request
        def waiter = getEc2Client().waiters().spotInstanceRequestFulfilled()
        def describeInstances = new DescribeSpotInstanceRequestsRequest().withSpotInstanceRequestIds(spotIds)
        waiter.run( new WaiterParameters<>().withRequest(describeInstances)  )

        def result = getEc2Client().describeSpotInstanceRequests(describeInstances)
        return result.spotInstanceRequests *. instanceId
    }

    /**
     * Convert a list of {@link Tag} to a {@link Map} object
     *
     * @param tags
     *      A list key=value {@link Tag}s
     * @return
     *      A map representing the specified key-value tags
     */
    private Map<String,String> tagListToMap( List<Tag> tags ) {
        def result = [:]
        tags.each { tag ->
            result.put( tag.key, tag.value )
        }
        return result
    }

    /**
     * Converts an AWS {@link Instance} model object to a {@link CloudInstance} object
     *
     * @param instance
     *      An {@link Instance} object representing a running instance
     * @return
     *      The {@link CloudInstance} object corresponding to the specified instance
     */
    private CloudInstance instanceToModel( Instance instance ) {

        def tags = tagListToMap(instance.tags)

        new CloudInstance(
                id: instance.instanceId,
                state: instance.state.name,
                publicDnsName: instance.publicDnsName,
                privateDnsName: instance.privateDnsName,
                publicIpAddress: instance.publicIpAddress,
                privateIpAddress: instance.privateIpAddress,
                role: tags.get(TAG_CLUSTER_ROLE),
                clusterName: tags.get(TAG_CLUSTER_NAME)
        )
    }

    /**
     * Iterate over the list of instances for the given tags
     *
     * @param tags
     *      One or more tag given as a {@link Map} object
     * @param callback
     *      A closure getting a single parameter of type {@link nextflow.cloud.types.CloudInstance}
     *      describing the properties for the current instance
     *
     * @see nextflow.cloud.types.CloudInstance
     */
    void eachInstanceWithTags(
            Map tags,
            @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudInstance']) Closure callback)
    {
        eachInstance0(null, tags, callback)
    }

    /**
     * Iterate over the list of instances for the given instance IDs
     *
     * @param instanceIds
     *      One or more instance IDs
     * @param callback
     *      A closure getting a single parameter of type {@link nextflow.cloud.types.CloudInstance}
     *      describing the properties for the current instance
     *
     * @see nextflow.cloud.types.CloudInstance
     */
    void eachInstanceWithIds(
            List<String> instanceIds,
            @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudInstance']) Closure callback)
    {
        eachInstance0(instanceIds, null, callback)
    }

    /**
     * Iterate over the list of available EC2 instances
     *
     * @param callback
     *      A closure getting a single parameter of type {@link nextflow.cloud.types.CloudInstance},
     *      representing the action to be applied on each instance
     *
     * @see nextflow.cloud.types.CloudInstance
     */
    void eachInstance( @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudInstance']) Closure callback) {
        eachInstance0(null, null, callback)
    }

    /**
     * Retrieve the private IPs of all nodes in the cluster with the specified name
     *
     * @param clusterName
     *      The name of the cluster for which retrieve the list of IPs
     * @return
     *      The list of node IPs in the specified cluster
     */
    @Override
    List<String> listPrivateIPs(String clusterName) {
        def result = []

        eachInstanceWithTags([(TAG_CLUSTER_NAME): clusterName]) { CloudInstance it ->
            result << it.privateIpAddress
        }

        return result
    }

    /**
     * Implements cluster instances iterator
     *
     * @param instanceIds
     *      An optional list of EC2 instance IDs
     * @param tags
     *      An optional map of tags filter the list of instances
     * @param callback
     *      A closure getting a single parameter of type {@link nextflow.cloud.types.CloudInstance},
     *      representing the action to be applied on each instance
     */
    private void eachInstance0(
            List<String> instanceIds,
            Map<String,String> tags,
            @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudInstance']) Closure callback)
    {
        def request = new DescribeInstancesRequest()

        // -- set filters
        def filters = []
        def running = new Filter('instance-state-name', ['running'])
        filters << running

        tags?.each { k, v ->
            filters << new Filter('tag-key', [k])
            if( v ) {
                def values = (v instanceof List<String> ? v as List : [v])
                filters << new Filter('tag-value', values)
            }
        }

        request.setFilters(filters)

        // -- set instance ids
        if( instanceIds )
            request.setInstanceIds(instanceIds)

        // -- submit the request
        def result = getEc2Client().describeInstances(request)
        def itr1 = result.reservations.iterator()
        while( itr1.hasNext() ) {
            def reservation = itr1.next()
            def itr2 = reservation.instances.iterator()
            while( itr2.hasNext() ) {
                def instance = itr2.next()
                callback.call(instanceToModel(instance))
            }
        }

    }

    /**
     * Wait for the specified instances reach the specified {@link CloudInstanceStatus}
     *
     * @param instanceIds One or more EC2 instance IDs
     */
    @Override
    void waitInstanceStatus(Collection<String> instanceIds, CloudInstanceStatus status) {

        switch( status ) {
            case CloudInstanceStatus.STARTED:
                waitRunning(instanceIds)
                break

            case CloudInstanceStatus.READY:
                waitRunning(instanceIds)
                waitStatusOk(instanceIds)
                break

            case CloudInstanceStatus.TERMINATED:
                waitTerminated(instanceIds)
                break

            default:
                throw new IllegalStateException("Unknown instance status: $status")

        }
    }

    /**
     * Wait for the specified instances reach RUNNING status
     *
     * @param instanceIds One or more EC2 instance IDs
     */
    @PackageScope
    void waitRunning( Collection<String> instancesId ) {
        def waiter = getEc2Client().waiters().instanceRunning()
        def describeInstances = new DescribeInstancesRequest().withInstanceIds(instancesId)
        waiter.run( new WaiterParameters<>().withRequest(describeInstances)  )
    }

    /**
     * Wait for the specified instances reach OK status
     *
     * @param instanceIds One or more EC2 instance IDs
     */
    @PackageScope
    void waitStatusOk( Collection<String> instanceIds ) {
        def waiter = getEc2Client().waiters().instanceStatusOk()
        def describeInstances = new DescribeInstanceStatusRequest().withInstanceIds(instanceIds)
        waiter.run( new WaiterParameters<>().withRequest(describeInstances) )
    }

    /**
     * Wait for the specified instances reach a termination status
     *
     * @param instanceIds One or more EC2 instance IDs
     */
    @PackageScope
    void waitTerminated( Collection<String> instanceIds ) {
        def waiter = getEc2Client().waiters().instanceTerminated()
        def describeInstances = new DescribeInstancesRequest().withInstanceIds(instanceIds)
        waiter.run( new WaiterParameters<>().withRequest(describeInstances) )
    }

    /**
     * Tag one or more instances with the specified key=value pairs
     *
     * @param instanceIds A list of instance IDs
     * @param tags A mpa of tags to be associated to the specified instances
     */
    @Override
    void tagInstances(Collection<String> instanceIds, Map<String,String> tags ) {

        if( !tags ) return

        def req = new CreateTagsRequest()
                .withResources(instanceIds)
                .withTags( tags.collect { k, v -> new Tag(k, v) } )

        getEc2Client().createTags(req)
    }

    /**
     * Iterate over list of spot prices for the given instance types
     *
     * @param instanceTypes
     *      A collection of instance types
     * @param callback
     *      A closure getting a single parameter of type {@link nextflow.cloud.types.CloudSpotPrice} describing
     *      the instance properties and price for the actual instance type
     */
    void eachSpotPrice(List<String> instanceTypes, Closure action) {

        def req = new DescribeSpotPriceHistoryRequest()
        if( instanceTypes ) {
            req.setInstanceTypes( instanceTypes as Collection<String> )
        }

        def history = getEc2Client().describeSpotPriceHistory(req).getSpotPriceHistory()
        for( SpotPrice entry : history ) {
            final price = new CloudSpotPrice(
                    type: entry.instanceType,
                    price: entry.spotPrice,
                    zone: entry.availabilityZone,
                    description: entry.productDescription,
                    timestamp: entry.timestamp
            )
            action.call(price)
        }
    }

    /**
     * Terminate one or more EC2 instances
     *
     * @param instanceIds A list of instance IDs
     */
    void terminateInstances( Collection<String> instanceIds ) {
        log.debug "Terminating EC2 instances: ids=${instanceIds?.join(',')}"
        if( !instanceIds ) return
        getEc2Client().terminateInstances(new TerminateInstancesRequest().withInstanceIds(instanceIds))
    }

    /**
     * Validate the cloud configuration object
     *
     * @param config The {@link LaunchConfig} configuration object
     */
    void validate(LaunchConfig config) {
        if( !config.imageId )
            throw new IllegalArgumentException("Missing mandatory cloud `imageId` setting")

        if( !config.instanceType )
            throw new IllegalStateException("Missing mandatory cloud `instanceType` setting")

        if( !describeInstanceType(config.instanceType) )
            throw new IllegalArgumentException("Unknown EC2 instance type: ${config.instanceType}")

        if( getAccessKey() && getSecretKey() ) {
            log.debug "AWS authentication based on access/secret credentials"
        }
        else {
            def role = config.instanceRole
            if( !role )
                config.instanceRole = role = fetchIamRole()
            if( !role )
                throw new IllegalArgumentException("Missing auth credentials -- Provide the accessKey/secretKey or IAM instance profile")
            log.debug "AWS authentication based IAM instance-role: `$role`"
        }
    }

    /**
     * @return The current instance identified
     */
    @Override
    String getLocalInstanceId() {
        try {
            return getUrl('http://169.254.169.254/latest/meta-data/instance-id')
        }
        catch( Exception e ) {
            log.debug "Oops.. Can't discover EC2 instance-id -- Cause: ${e.message ?: e}"
            return null
        }
    }

    /**
     * Retrieve a termination notice for spot instance notifying the instance is
     * going to be retired.
     *
     *  https://aws.amazon.com/blogs/aws/new-ec2-spot-instance-termination-notices/
     *  http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/spot-interruptions.html
     *
     * @return
     *  The termination notice string as returned by the instance metadata
     */
    @Override
    String getLocalTerminationNotice() {
        try {
            return getUrl('http://169.254.169.254/latest/meta-data/spot/termination-time')
        }
        catch (Exception e) {
            return null
        }
    }

    /**
     * Fetch a remote URL resource text content
     *
     * @param path
     *      A valid http/https resource URL
     * @param timeout
     *      Max connection timeout in millis
     * @return
     *      The resource URL content
     */
    protected String getUrl(String path, int timeout=150) {
        final url = new URL(path)
        final con = url.openConnection()
        con.setConnectTimeout(timeout)
        con.setReadTimeout(timeout)
        return con.getInputStream().text.trim()
    }

    /**
     * Describe an instance type by the given type
     *
     * Note: it's declared as Memoized to avoid to make a remote API call
     * for each invocation
     *
     * @param instanceType The instance type ID e.g. m3.2xlarge
     * @return
     *      The {@link CloudInstanceType} instance for the given
     *      instance or {@code null} if no record is found for the given instance type
     */
    @Memoized
    CloudInstanceType describeInstanceType( String instanceType ) {
        createReader(region).getInstanceType(instanceType)
    }

    /**
     * Creates a {@link AmazonPriceReader} for the specified region
     *
     * @param region An AWS region identifier eg {@code eu-west-1}
     * @return An {@link AmazonPriceReader} instance
     */
    @PackageScope
    AmazonPriceReader createReader(String region) {
        new AmazonPriceReader(region)
    }

}
