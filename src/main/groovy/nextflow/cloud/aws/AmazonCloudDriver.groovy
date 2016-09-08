/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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
import static nextflow.cloud.CloudConst.TAG_CLUSTER_NAME
import static nextflow.cloud.CloudConst.TAG_CLUSTER_ROLE

import javax.mail.Session
import javax.mail.internet.MimeBodyPart
import javax.mail.internet.MimeMessage
import javax.mail.internet.MimeMultipart

import com.amazonaws.auth.BasicAWSCredentials
import com.amazonaws.regions.RegionUtils
import com.amazonaws.services.ec2.AmazonEC2Client
import com.amazonaws.services.ec2.model.BlockDeviceMapping
import com.amazonaws.services.ec2.model.CreateTagsRequest
import com.amazonaws.services.ec2.model.DescribeImagesRequest
import com.amazonaws.services.ec2.model.DescribeInstanceStatusRequest
import com.amazonaws.services.ec2.model.DescribeInstancesRequest
import com.amazonaws.services.ec2.model.DescribeSpotInstanceRequestsRequest
import com.amazonaws.services.ec2.model.DescribeSpotPriceHistoryRequest
import com.amazonaws.services.ec2.model.Filter
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

    AmazonEC2Client ec2client

    private Map config

    private String accessKey

    private String secretKey

    private String region

    @CompileDynamic
    AmazonCloudDriver() {
        this.config = Global.config
        // -- get the aws credentials
        def credentials = Global.getAwsCredentials()
        if( !credentials )
            throw new AbortOperationException('Missing AWS access and secret keys -- Make sure to define in your system environment the following variables `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY`')
        this.accessKey = credentials[0]
        this.secretKey = credentials[1]
        // -- get the aws default region
        this.region = Global.getAwsRegion()
        if( !region )
            throw new AbortOperationException('Missing AWS region -- Make sure to deifne in your system environment the variable `AWS_DEFAULT_REGION`')
    }

    @CompileDynamic
    protected AmazonCloudDriver( Map config, AmazonEC2Client client ) {
        this.config = config
        (accessKey, secretKey) = Global.getAwsCredentials(null, config)
        this.region = Global.getAwsRegion(null, config)
        this.ec2client = client
    }

    static private AmazonEC2Client createEc2Client( String accessKey, String secretKey, String region ) {
        def result = new AmazonEC2Client(new BasicAWSCredentials(accessKey, secretKey))
        result.setRegion( RegionUtils.getRegion(region) )
        return result
    }

    synchronized private AmazonEC2Client getClient() {
        if( ec2client == null ) {
            if( !accessKey ) throw new IllegalArgumentException("Missing AWS access key config property")
            if( !secretKey ) throw new IllegalArgumentException("Missing AWS secret key config property")
            if( !region ) throw new IllegalArgumentException("Missing AWS region config property")

            ec2client = createEc2Client(accessKey, secretKey, region)
        }
        return ec2client
    }

    @PackageScope
    String scriptCreateUser(String userName, String key) {
        """\
        useradd $userName
        mkdir ~$userName/.ssh
        echo "${key.trim()}" > ~$userName/.ssh/authorized_keys
        chmod 700 ~$userName/.ssh
        chmod 600 ~$userName/.ssh/authorized_keys
        chown -R $userName:$userName ~$userName/.ssh
        egrep -i "^wheel:" /etc/group > /dev/null && usermod -aG wheel $userName
        egrep -i "^docker:" /etc/group > /dev/null && usermod -aG docker $userName
        chmod +x /etc/sudoers
        echo '${userName} ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
        chmod -x /etc/sudoers
        """
        .stripIndent()
    }

    @PackageScope
    String scriptMountEFS(String fileSystemId, String fileSystemMount, String userName) {

        """\
        zone="\$(curl -s http://169.254.169.254/latest/meta-data/placement/availability-zone)"
        region="\${zone::-1}"
        yum install -y nfs-utils
        mkdir -p $fileSystemMount
        mount -t nfs4 -o nfsvers=4.1 \${zone}.${fileSystemId}.efs.\${region}.amazonaws.com:/ $fileSystemMount
        chown ${userName}:${userName} $fileSystemMount
        chmod 775 $fileSystemMount
        """
        .stripIndent()
    }

    @PackageScope
    String scriptMountInstanceStorage(String devicePath, String mountPath, String user) {

        """\
        mkfs.ext4 -E nodiscard $devicePath
        mkdir -p $mountPath
        mount -o discard $devicePath $mountPath
        chown -R $user:$user $mountPath
        chmod 775 $mountPath
        """
        .stripIndent()
    }

    @PackageScope
    String scriptBashEnv( LaunchConfig cfg ) {
        def profile = """\
        export AWS_ACCESS_KEY_ID='$accessKey'
        export AWS_SECRET_ACCESS_KEY='$secretKey'
        export AWS_DEFAULT_REGION='$region'
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

        if( cfg.instanceStorageDevice && cfg.instanceStorageMount )
            profile += "export NXF_TEMP='${cfg.instanceStorageMount}'\n"

        return profile
    }


    @PackageScope
    String cloudInitScript(LaunchConfig cfg) {
        if( !accessKey ) throw new AbortOperationException('Missing AWS access id')
        if( !secretKey ) throw new AbortOperationException('Missing AWS secret key')
        if( !region ) throw new AbortOperationException('Missing AWS region')

        // load init script template
        def template =  this.class.getResourceAsStream('cloud-boot.txt')
        if( !template )
            throw new IllegalStateException("Missing `cloud-boot.txt` template resource")

        def binding = new CloudBootTemplateBinding(cfg)
        binding.awsAccessKey = accessKey
        binding.awsSecretKey = secretKey
        binding.awsRegion = region
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
     * Creates the user-data script that bootstrap the Ec2 instance.
     *
     * Note: the script is encoded as mime-multipart message
     *
     * See
     *      http://cloudinit.readthedocs.io/en/latest/topics/format.html#mime-multi-part-archive
            https://forums.aws.amazon.com/thread.jspa?messageID=730555&#730555
     * @param cfg
     * @return
     */
    @PackageScope
    String getUserDataAsString(LaunchConfig cfg) {
        new String( getUserData0(cfg) )
    }

    @PackageScope
    String getUserDataAsBase64( LaunchConfig cfg ) {
        Base64.encodeAsString(getUserData0(cfg))
    }

    @PackageScope
    String cloudBootHookScript( LaunchConfig cfg ) {
        def builder = []

        if( cfg.createUser ) {
            builder << scriptCreateUser(cfg.userName, cfg.keyHash)
        }

        if( cfg.sharedStorageId && cfg.sharedStorageMount ) {
            builder << scriptMountEFS(cfg.sharedStorageId, cfg.sharedStorageMount, cfg.userName)
        }

        if( cfg.instanceType?.startsWith('r3.') && cfg.instanceStorageDevice && cfg.instanceStorageMount ) {
            builder << scriptMountInstanceStorage(cfg.instanceStorageDevice, cfg.instanceStorageMount, cfg.userName)
        }

        builder.join('\n')
    }

    private byte[] getUserData0( LaunchConfig cfg ) {

        MimeMultipart multipart = new MimeMultipart()

        //
        // the first part contains the cloud-init script
        //
        def script1 = new MimeBodyPart()
        script1.setContent( cloudInitScript(cfg), 'text/x-shellscript' )
        multipart.addBodyPart(script1)

        //
        // second part mount the instance and EFS storage
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

    @PackageScope
    List<BlockDeviceMapping> getBlockDeviceMappings(LaunchConfig cfg) {
        def result = new ArrayList<BlockDeviceMapping>()

        // -- map the ephemeral storage
        if( cfg.instanceStorageDevice ) {
            def mapping = new BlockDeviceMapping()
            mapping.deviceName = cfg.instanceStorageDevice
            mapping.virtualName = 'ephemeral0'
            result << mapping
        }

        // -- resize the boot storage
        if( cfg.bootStorageSize ) {
            def rootDevice = getRooDeviceMapping( cfg.imageId )
            rootDevice.ebs.volumeSize = (int)cfg.bootStorageSize.toGiga()
            rootDevice.ebs.encrypted = null // <-- note: reset to null the encrypted flag
            log.debug "Root EBS device: ${rootDevice.ebs}"
            result << rootDevice
        }

        return result ?: null
    }

    BlockDeviceMapping getRooDeviceMapping( String imageId ) {
        final req = new DescribeImagesRequest()
        req.setImageIds( [imageId] )
        final Image image = getClient().describeImages(req).getImages().get(0)
        return image.blockDeviceMappings.find { it.deviceName == image.rootDeviceName }
    }

    @PackageScope
    RunInstancesRequest makeRunRequest( int instanceCount, LaunchConfig cfg ) {
        def req = new RunInstancesRequest()
        req.minCount = instanceCount
        req.maxCount = instanceCount
        req.imageId = cfg.imageId
        req.instanceType = cfg.instanceType
        req.userData = getUserDataAsBase64(cfg)
        req.setBlockDeviceMappings( getBlockDeviceMappings(cfg) )

        if( cfg.keyName )
            req.keyName = cfg.keyName

        if( cfg.securityGroups )
            req.securityGroupIds = cfg.securityGroups

        if( cfg.subnetId )
            req.subnetId = cfg.subnetId

        return req
    }

    @PackageScope
    RequestSpotInstancesRequest makeSpotRequest( int instanceCount, LaunchConfig cfg ) {

        if( !cfg.imageId ) throw new AbortOperationException("Missing `imageId` setting in cloud configuration")
        if( !cfg.instanceType ) throw new AbortOperationException('Missing `instanceType` setting in cloud configuration')

        def spec = new LaunchSpecification()
        cfg.with {
            spec.setImageId(imageId)
            spec.setInstanceType(instanceType)
            spec.setUserData( getUserDataAsBase64(cfg) )

            if( keyName )
                spec.keyName = keyName

            if( subnetId )
                spec.subnetId = subnetId

            if( securityGroups )
                spec.securityGroups = securityGroups
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
        RunInstancesResult response = getClient().runInstances(req)
        return response.reservation.getInstances() *. instanceId
    }

    private List<String> launchSpotInstances0(int instanceCount, LaunchConfig cfg) {
        log.debug "AWS submitting launch request for $instanceCount SPOT instances; config: $cfg"

        def req = makeSpotRequest(instanceCount, cfg)
        RequestSpotInstancesResult response = getClient().requestSpotInstances(req)
        def spotIds = response.spotInstanceRequests *. spotInstanceRequestId
        log.debug "AWS spot request IDs: ${spotIds.join(',')}"

        // -- wait for fulfill the spot request
        def waiter = getClient().waiters().spotInstanceRequestFulfilled()
        def describeInstances = new DescribeSpotInstanceRequestsRequest().withSpotInstanceRequestIds(spotIds)
        waiter.run( new WaiterParameters<>().withRequest(describeInstances)  )

        def result = getClient().describeSpotInstanceRequests(describeInstances)
        return result.spotInstanceRequests *. instanceId
    }

    private Map<String,String> tagListToMap( List<Tag> tags ) {
        def result = [:]
        tags.each { tag ->
            result.put( tag.key, tag.value )
        }
        return result
    }

    private CloudInstance instanceToModel( Instance it ) {

        def tags = tagListToMap(it.tags)

        new CloudInstance(
                id: it.instanceId,
                state: it.state.name,
                publicDnsName: it.publicDnsName,
                privateDnsName: it.privateDnsName,
                publicIpAddress: it.publicIpAddress,
                privateIpAddress: it.privateIpAddress,
                role: tags.get(TAG_CLUSTER_ROLE),
                clusterName: tags.get(TAG_CLUSTER_NAME)
        )
    }

    void eachInstanceWithTags(Map tags, Closure callback) {
        eachInstance0(null, tags, callback)
    }

    void eachInstanceWithIds(List<String> instanceIds, Closure callback) {
        eachInstance0(instanceIds, null, callback)
    }


    void eachInstance(Closure callback) {
        eachInstance0(null, null, callback)
    }

    @Override
    List<String> listPrivateIPs(String clusterName) {
        def result = []

        eachInstanceWithTags([(TAG_CLUSTER_NAME): clusterName]) { CloudInstance it ->
            result << it.privateIpAddress
        }

        return result
    }

    private void eachInstance0(List<String> instanceIds, Map<String,String> tags, Closure callback) {
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
        def result = getClient().describeInstances(request)
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

    @PackageScope
    void waitRunning( Collection<String> instancesId ) {
        def waiter = getClient().waiters().instanceRunning()
        def describeInstances = new DescribeInstancesRequest().withInstanceIds(instancesId)
        waiter.run( new WaiterParameters<>().withRequest(describeInstances)  )
    }

    @PackageScope
    void waitStatusOk( Collection<String> instanceIds ) {
        def waiter = getClient().waiters().instanceStatusOk()
        def describeInstances = new DescribeInstanceStatusRequest().withInstanceIds(instanceIds)
        waiter.run( new WaiterParameters<>().withRequest(describeInstances) )
    }

    @PackageScope
    void waitTerminated( Collection<String> instanceIds ) {
        def waiter = getClient().waiters().instanceTerminated()
        def describeInstances = new DescribeInstancesRequest().withInstanceIds(instanceIds)
        waiter.run( new WaiterParameters<>().withRequest(describeInstances) )
    }


    @Override
    void tagInstances(Collection<String> instanceIds, Map<String,String> tags ) {

        if( !tags ) return

        def req = new CreateTagsRequest()
                .withResources(instanceIds)
                .withTags( tags.collect { k, v -> new Tag(k, v) } )

        getClient().createTags(req)
    }



    void eachSpotPrice(List<String> instanceTypes, Closure action) {

        def req = new DescribeSpotPriceHistoryRequest()
        if( instanceTypes ) {
            req.setInstanceTypes( instanceTypes as Collection<String> )
        }

        def history = getClient().describeSpotPriceHistory(req).getSpotPriceHistory()
        history.each { SpotPrice entry ->

            def price = new CloudSpotPrice(
                    type: entry.instanceType,
                    price: entry.spotPrice,
                    zone: entry.availabilityZone,
                    description: entry.productDescription,
                    timestamp: entry.timestamp
            )
            action.call(price)

        }
    }

    void terminateInstances( Collection<String> instanceIds ) {
        log.debug "Terminating EC2 instances: ids=${instanceIds?.join(',')}"
        if( !instanceIds ) return
        getClient().terminateInstances(new TerminateInstancesRequest().withInstanceIds(instanceIds))
    }

    void validate(LaunchConfig config) {
        if( !config.imageId ) throw new IllegalArgumentException("Missing mandatory cloud `imageId` setting")
        if( !config.instanceType ) throw new IllegalStateException("Missing mandatory cloud `instanceType` setting")
    }

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
     * see
     *   https://aws.amazon.com/blogs/aws/new-ec2-spot-instance-termination-notices/
     *   http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/spot-interruptions.html
     *
     * @return
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

    private String getUrl(String path, int timeout=150) {
        final url = new URL(path)
        final con = url.openConnection()
        con.setConnectTimeout(timeout)
        con.setReadTimeout(timeout)
        return con.getInputStream().text.trim()
    }

    /**
     * Describe an instance type by the given type
     *
     * @param instanceType The instance type ID e.g. m3.2xlarge
     * @return The {@link CloudInstanceType} instance for the given instance
     * @throws IllegalArgumentException If the the specified instance type is unknown
     */
    @Memoized
    CloudInstanceType describeInstanceType( String instanceType ) {
        def result = new AmazonPriceReader(region).instanceTypeTable.get(instanceType)
        if( !result )
            throw new IllegalArgumentException("Unknown EC2 instance type: `$instanceType`")
        return result
    }

}
