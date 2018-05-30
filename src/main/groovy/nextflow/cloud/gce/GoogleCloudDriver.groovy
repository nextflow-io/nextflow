package nextflow.cloud.gce

import com.google.api.services.compute.Compute
import com.google.api.services.compute.model.AccessConfig
import com.google.api.services.compute.model.Instance
import com.google.api.services.compute.model.InstancesSetLabelsRequest
import com.google.api.services.compute.model.NetworkInterface
import com.google.api.services.compute.model.Operation
import com.google.common.collect.Maps
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.SimpleType
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.cloud.CloudDriver
import nextflow.cloud.LaunchConfig
import nextflow.cloud.types.CloudInstance
import nextflow.cloud.types.CloudInstanceStatus
import nextflow.cloud.types.CloudInstanceType
import nextflow.exception.AbortOperationException
import nextflow.util.ServiceName

import java.security.GeneralSecurityException

import static nextflow.cloud.CloudConst.TAG_CLUSTER_NAME
import static nextflow.cloud.CloudConst.TAG_CLUSTER_ROLE

@Slf4j
@CompileStatic
@ServiceName('gce')
/**
 * Cloud driver implementation for Google Compute Engine
 *
 * @author Vilmundur Pálmason <vilmundur@wuxinextcode.com>
 */
class GoogleCloudDriver implements CloudDriver {

    /**
     * The GCE zone eg. {@code us-central1-f}. If it's not specified the current region is retrieved from
     * the GCE instance metadata
     */
    private String zone
    private String project
    private String group = "inst"

    private GceApiHelper helper

    static long OPS_WAIT_TIMEOUT_MS = 5*60*1000

    /**
     * Initialise the Google cloud driver with default (empty) parameters
     */
    GoogleCloudDriver() {
        this(Collections.emptyMap())
    }

    /**
     * Initialise the Google cloud driver with the specified parameters
     *
     * @param config
     *      A map holding the driver parameters:
     *      - zone: the GCE zone
     *      - project: GCE project id
     */
    @CompileDynamic
    GoogleCloudDriver(Map config) {
        log.debug("Config: {}",config)
        log.debug("Global config: {}",Global.getConfig())
        this.zone = config.zone ?: Global.getConfig().gce.zone
        this.project = config.project ?: Global.getConfig().gce.project
        if (!(this.zone && this.project)) {
            throw new AbortOperationException("Need GCE project and region")
        }
        log.debug("Starting GoogleCloudDriver in project {} and zone {}",this.project,this.zone)
        this.helper = new GceApiHelper(project,zone)
    }

    /**
     * Gets {@link Compute} instance given the current
     * configuration parameter
     *
     * @return
     *      An {@link Compute} instance
     */
    synchronized Compute getClient() {
        return helper.compute()
    }

    @Override
    void validate(LaunchConfig config) {
        if( !config.imageId )
            throw new AbortOperationException("Missing mandatory cloud `imageId` setting")

        if( !config.instanceType )
            throw new AbortOperationException("Missing mandatory cloud `instanceType` setting")

        if( !helper.lookupMachineType(config.instanceType) )
            throw new AbortOperationException("Unknown GCE machine type: ${config.instanceType}")

        String validationError = helper.validateLabelValue(config.getClusterName());
        if (validationError != null) {
            throw new AbortOperationException("Invalid cluster name '"+config.getClusterName()+"': "+validationError);
        }
    }

    @Override
    List<String> launchInstances(int instanceCount, LaunchConfig config) {
        List<String> result = []
        List<Operation> ops = []
        instanceCount.times {
            def inst = new Instance();
            inst.setName(helper.randomName(config.getClusterName() + "-"));
            inst.setMachineType(helper.instanceType(config.getInstanceType()));

            helper.setBootDisk(inst,config.getImageId());
            helper.setNetworkInterface(inst);
            helper.setStartupScript(inst,config)
            def insert = helper.compute().instances().insert(project, zone, inst);

            result << inst.getName()
            ops << insert.execute()
        }
        helper.blockUntilComplete(ops,OPS_WAIT_TIMEOUT_MS);
        return result
    }

    @Override
    void waitInstanceStatus(Collection<String> instanceIds, CloudInstanceStatus status) {
        unsupported("waitInstanceStatus")
    }

    @Override
    void tagInstances(Collection<String> instanceIds, Map<String, String> tags) {
        Map<String,String> labels = [:]
        tags.each {k,v -> labels[k.toLowerCase()] = v}

        List<Operation> ops = []
        for (String instanceId: instanceIds) {
            def instance = helper.compute().instances().get(project,zone,instanceId).execute()

            // Preserve existing labels
            if (instance.getLabels() != null) {
                labels = labels + instance.getLabels()
            }

            def request = new InstancesSetLabelsRequest()
            request.setLabelFingerprint(instance.getLabelFingerprint())
            request.setLabels(labels)

            ops << helper.compute().instances().setLabels(project,zone,instanceId,request).execute()
        }
        helper.blockUntilComplete(ops,OPS_WAIT_TIMEOUT_MS);
    }

    @Override
    void eachSpotPrice(List<String> instanceTypes,
                       @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudSpotPrice']) Closure callback) {
        unsupported("eachSpotPrice")

    }

    @Override
    void eachInstanceWithTags(Map tags,
                              @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudInstance']) Closure callback) {
        def listRequest = helper.compute().instances().list(project, zone)

        if (tags != null && !tags.isEmpty()) {
            listRequest.setFilter(tags.collect(this.&tagToFilterExpression).join(" "))
        }
        listRequest.execute().getItems()?.each {inst -> callback.call(toNextflow(inst))}
    }

    @Override
    void eachInstanceWithIds(List<String> instanceIds,
                             @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudInstance']) Closure callback) {
        if (instanceIds.size() > 0) {
            def listRequest = helper.compute().instances().list(project, zone)
            listRequest.setFilter(instanceIds.collect(this.&instanceIdToFilterExpression).join(" or "))
            listRequest.execute().getItems()?.each { inst -> callback.call(toNextflow(inst)) }
        }
    }

    @Override
    void eachInstance(
            @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudInstance']) Closure callback) {
        unsupported("eachInstance")
    }

    @Override
    List<String> listPrivateIPs(String clusterName) {
        unsupported("listPrivateIPs")
        return null
    }

    @Override
    void terminateInstances(Collection<String> instanceIds) {
        unsupported("terminateInstances")
    }

    @Override
    String getLocalInstanceId() {
        unsupported("getLocalInstanceId")
    }

    @Override
    String getLocalTerminationNotice() {
        unsupported("getLocalTerminationNotice")
    }

    @Override
    CloudInstanceType describeInstanceType(String instanceType) {
        unsupported("describeInstanceType")
        return null
    }

    /**
     * @TODO: This method will be removed once all methods are implemented
     */
    def unsupported(String msg) {
        log.warn("UNSUPPORTED: "+msg)
    }

    def tagToFilterExpression(String k,v) {
        '(labels.' + k.toLowerCase() + '= "' + (v?:'*') + '")'
    }

    def instanceIdToFilterExpression(instanceId) {
        '(name = "' + instanceId + '")'
    }

    def toNextflow(Instance instance) {
        NetworkInterface iface
        AccessConfig accessConfig
        def labels = instance.getLabels() ?: {}

        if (instance.getNetworkInterfaces() != null && !instance.getNetworkInterfaces().isEmpty()) {
            iface = instance.getNetworkInterfaces()[0]
            if (iface.getAccessConfigs() != null && !iface.getAccessConfigs()?.isEmpty()) {
                accessConfig = iface.getAccessConfigs()[0]
            }
        }

        new CloudInstance(
                id: instance.getName(),
                privateIpAddress: iface?.getNetworkIP(),
                privateDnsName: helper.instanceIdToPrivateDNS(instance.getName()),
                publicIpAddress: accessConfig?.getNatIP(),
                publicDnsName: accessConfig?.getPublicPtrDomainName() ?: helper.publicIpToDns(accessConfig?.getNatIP()),
                state: instance.getStatus(),
                role: labels[TAG_CLUSTER_ROLE.toLowerCase()],
                clusterName: labels[TAG_CLUSTER_NAME.toLowerCase()]
        )
    }


}
