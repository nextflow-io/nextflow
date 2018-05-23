package nextflow.cloud.gce

import com.google.api.services.compute.Compute
import com.google.api.services.compute.model.Instance
import com.google.api.services.compute.model.Operation
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.SimpleType
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.cloud.CloudDriver
import nextflow.cloud.LaunchConfig
import nextflow.cloud.types.CloudInstanceStatus
import nextflow.cloud.types.CloudInstanceType
import nextflow.exception.AbortOperationException
import nextflow.util.ServiceName

import java.security.GeneralSecurityException

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
            throw new IllegalArgumentException("Missing mandatory cloud `imageId` setting")

        if( !config.instanceType )
            throw new IllegalStateException("Missing mandatory cloud `instanceType` setting")

        if( !helper.lookupMachineType(config.instanceType) )
            throw new IllegalArgumentException("Unknown GCE machine type: ${config.instanceType}")

    }

    @Override
    List<String> launchInstances(int instanceCount, LaunchConfig config) {
        def result = new ArrayList<>()
        def ops = new ArrayList<Operation>()
        instanceCount.times {
            Instance inst = new Instance();
            inst.setName(helper.randomName(config.getClusterName() + "-"));
            inst.setMachineType(helper.instanceType(config.getInstanceType()));
            helper.setBootDisk(inst,config.getImageId());
            helper.setNetworkInterface(inst);
            Compute.Instances.Insert insert = helper.compute().instances().insert(project, zone, inst);

            result << inst.getName()
            ops << insert.execute()
            System.out.println("Launching: "+result)
        }
        helper.blockUntilComplete(ops,5*60*1000);
        System.out.println("Launched: "+result)
        return result
    }

    @Override
    void waitInstanceStatus(Collection<String> instanceIds, CloudInstanceStatus status) {
        unsupported("waitInstanceStatus")
    }

    @Override
    void tagInstances(Collection<String> instanceIds, Map<String, String> tags) {
        unsupported("tagInstances")
    }

    @Override
    void eachSpotPrice(List<String> instanceTypes,
                       @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudSpotPrice']) Closure callback) {
        unsupported("eachSpotPrice")

    }

    @Override
    void eachInstanceWithTags(Map tags,
                              @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudInstance']) Closure callback) {
        unsupported("eachInstanceWithTags")

    }

    @Override
    void eachInstanceWithIds(List<String> instanceIds,
                             @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudInstance']) Closure callback) {
        unsupported("eachInstanceWithIds")
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
    private void unsupported(String msg) {
        log.warn("UNSUPPORTED: "+msg)
    }
}
