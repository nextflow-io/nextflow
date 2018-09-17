package nextflow.cloud.gce

import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.google.api.services.compute.Compute
import com.google.api.services.compute.model.AccessConfig
import com.google.api.services.compute.model.Instance
import com.google.api.services.compute.model.InstancesSetLabelsRequest
import com.google.api.services.compute.model.NetworkInterface
import com.google.api.services.compute.model.Operation
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.SimpleType
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.cloud.CloudDriver
import nextflow.cloud.CloudScripts
import nextflow.cloud.LaunchConfig
import nextflow.cloud.aws.CloudBootTemplateBinding
import nextflow.cloud.types.CloudInstance
import nextflow.cloud.types.CloudInstanceStatus
import nextflow.cloud.types.CloudInstanceType
import nextflow.exception.AbortOperationException
import nextflow.processor.TaskTemplateEngine
import nextflow.util.ServiceName

import java.nio.file.Path

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
    private String group = "inst"

    private GceApiHelper helper

    static long OPS_WAIT_TIMEOUT_MS = 5 * 60 * 1000

    static long POLL_WAIT = 1000

    /**
     * Initialise the Google cloud driver with default (empty) parameters
     */
    GoogleCloudDriver(GoogleCredential credential = null) {
        this(Collections.emptyMap(), credential)
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
    GoogleCloudDriver(Map config, GoogleCredential credential = null) {
        log.debug("Config: {}", config)
        log.debug("Global config: {}", Global.getConfig())
        String zone = config.zone ?: Global.getConfig()?.gce?.zone
        String project = config.project ?: Global.getConfig()?.gce?.project
        this.helper = credential ? new GceApiHelper(project, zone, credential) : new GceApiHelper(project, zone)
        log.debug("Starting GoogleCloudDriver in project {} and zone {}", helper.project, helper.zone)
    }

    /**
     * Gets {@link Compute} instance given the current
     * configuration parameter
     *
     * @return
     * An {@link Compute} instance
     */
    synchronized Compute getClient() {
        helper.compute
    }

    @Override
    void validate(LaunchConfig config) {
        if (!config.imageId)
            throw new AbortOperationException("Missing mandatory cloud `imageId` setting")

        if (!config.instanceType)
            throw new AbortOperationException("Missing mandatory cloud `instanceType` setting")

        if (!helper.lookupMachineType(config.instanceType))
            throw new AbortOperationException("Unknown GCE machine type: ${config.instanceType}")

        String validationError = helper.validateLabelValue(config.clusterName)
        if (validationError != null) {
            throw new AbortOperationException("Invalid cluster name '" + config.clusterName + "': " + validationError)
        }

        if (config.sharedStorageId && config.sharedStorageMount) {
            throw new AbortOperationException("Shared storage not supported in Google Cloud")
        }
    }

    @Override
    List<String> launchInstances(int instanceCount, LaunchConfig config) {
        List<String> result = []
        List<Operation> ops = []
        instanceCount.times {
            Instance inst = new Instance();
            inst.setName(helper.randomName(config.getClusterName() + "-"))
            inst.setMachineType(helper.instanceType(config.getInstanceType()))
            inst.setScheduling(helper.createScheduling(config))


            //TODO make helper return the instances instead of setting them
            helper.setBootDisk(inst, config.getImageId())
            helper.setNetworkInterface(inst)
            helper.setStartupScript(inst, gceStartupScript(config))
            def insert = helper.compute.instances().insert(helper.project, helper.zone, inst)

            result << inst.getName()
            ops << insert.execute()
        }
        helper.blockUntilComplete(ops, OPS_WAIT_TIMEOUT_MS);
        return result
    }

    @Override
    void waitInstanceStatus(Collection<String> instanceIds, CloudInstanceStatus status) {
        def instanceStatusList

        switch (status) {
            case CloudInstanceStatus.STARTED:
                instanceStatusList = ['PROVISIONING', 'STAGING', 'RUNNING']
                waitStatus(instanceIds, instanceStatusList)
                break

            case CloudInstanceStatus.READY:
                instanceStatusList = ['RUNNING']
                waitStatus(instanceIds, instanceStatusList)
                break

            case CloudInstanceStatus.TERMINATED:
                instanceStatusList = ['TERMINATED']
                waitStatus(instanceIds, instanceStatusList)
                break

            default:
                throw new IllegalStateException("Unknown instance status: $status")
        }
    }

    @PackageScope
    void waitStatus(Collection<String> instanceIds, List<String> instanceStatusList) {

        Set<String> remaining = new HashSet<>(instanceIds)
        while (!remaining.isEmpty()) {
            def filter = instanceIds.collect(this.&instanceIdToFilterExpression).join(" OR ")
            def listRequest = helper.compute.instances().list(helper.project, helper.zone)
            listRequest.setFilter(filter)
            List<Instance> instances = listRequest.execute().getItems()
            if (instances != null && !instances.isEmpty()) {
                for (Instance instance : instances) {
                    if (instanceStatusList.contains(instance.status)) {
                        remaining.remove(instance.getName())
                    }
                }
            } else if (instanceStatusList.contains('TERMINATED') && (instances == null || instances.isEmpty())) {
                break
            }
            if (!remaining.isEmpty()) Thread.sleep(POLL_WAIT)
        }
    }

    @Override
    void tagInstances(Collection<String> instanceIds, Map<String, String> tags) {
        Map<String, String> labels = [:]
        tags.each { k, v -> labels[k.toLowerCase()] = v }

        List<Operation> ops = []
        for (String instanceId : instanceIds) {
            def instance = helper.compute.instances().get(helper.project, helper.zone, instanceId).execute()

            // Preserve existing labels
            if (instance.getLabels() != null) {
                labels = labels + instance.getLabels()
            }

            def request = new InstancesSetLabelsRequest()
            request.setLabelFingerprint(instance.getLabelFingerprint())
            request.setLabels(labels)

            ops << helper.compute.instances().setLabels(helper.project, helper.zone, instanceId, request).execute()
        }
        helper.blockUntilComplete(ops, OPS_WAIT_TIMEOUT_MS)
    }

    @Override
    void eachSpotPrice(List<String> instanceTypes,
                       @ClosureParams(value = SimpleType, options = ['nextflow.cloud.types.CloudSpotPrice']) Closure callback) {
        unsupported("eachSpotPrice")

    }

    @Override
    void eachInstanceWithTags(Map tags,
                              @ClosureParams(value = SimpleType, options = ['nextflow.cloud.types.CloudInstance']) Closure callback) {
        String filter = null

        if (tags != null && !tags.isEmpty()) {
            filter = tags.collect(this.&tagToFilterExpression).join(" ")
        }
        eachInstanceWithFilter(filter, callback)
    }

    @Override
    void eachInstanceWithIds(List<String> instanceIds,
                             @ClosureParams(value = SimpleType, options = ['nextflow.cloud.types.CloudInstance']) Closure callback) {
        if (instanceIds.size() > 0) {
            eachInstanceWithFilter(instanceIds.collect(this.&instanceIdToFilterExpression).join(" OR "), callback)
        }
    }

    @Override
    void eachInstance(
            @ClosureParams(value = SimpleType, options = ['nextflow.cloud.types.CloudInstance']) Closure callback) {
        eachInstanceWithFilter(null as String, callback)
    }

    void eachInstanceWithFilter(String filter, @ClosureParams(value = SimpleType, options = ['nextflow.cloud.types.CloudInstance']) Closure callback) {
        def listRequest = helper.compute.instances().list(helper.project, helper.zone)
        listRequest.setFilter(filter)
        listRequest.execute().getItems()?.each { inst -> callback.call(toNextflow(inst)) }
    }

    @Override
    List<String> listPrivateIPs(String clusterName) {
        def result = []

        eachInstanceWithTags([(TAG_CLUSTER_NAME): clusterName]) { CloudInstance it ->
            result << it.privateIpAddress
        }

        return result
    }

    @Override
    void terminateInstances(Collection<String> instanceIds) {
        for (String idInstance : instanceIds) {
            Compute.Instances.Delete delete = helper.compute.instances().delete(helper.project, helper.zone, idInstance)
            delete.execute()
        }
    }

    @Override
    String getLocalInstanceId() {
        helper.readInstanceId()
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
        log.warn("UNSUPPORTED: " + msg)
    }

    def tagToFilterExpression(String k, v) {
        '(labels.' + k.toLowerCase() + '= "' + (v ?: '*') + '")'
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

    /**
     * @TODO: This is mostly a copy paste from AmazonCloudDriver
     */
    @PackageScope
    String gceStartupScript(LaunchConfig cfg) {
        def builder = []

        if (cfg.createUser) {
            builder << CloudScripts.scriptCreateUser(cfg.userName, cfg.keyHash)
        }

        builder << cloudInitScript(cfg)

        if (builder) {
            builder.add(0, '#!/bin/bash')
        }
        // note: `findAll` remove all empty strings
        builder.join('\n')
    }

    def GCE_CREDENTIAL_FILE = '$HOME/.nextflow/gce_credentials.json'

    /**
     * @TODO: This is mostly a copy paste from AmazonCloudDriver
     */
    @PackageScope
    @CompileDynamic
    String cloudInitScript(LaunchConfig cfg) {
        // load init script template
        def template = this.class.getResourceAsStream('cloud-boot.txt')
        if (!template)
            throw new IllegalStateException("Missing `cloud-boot.txt` template resource")

        def binding = new CloudBootTemplateBinding(cfg)
        binding.nextflowConfig = cfg.renderCloudConfigObject()
        binding.bashProfile = scriptBashEnv(cfg)
        def credFile = System.getenv("GOOGLE_APPLICATION_CREDENTIALS") as Path
        if (credFile) {
            binding.gceCredentialsFile = GCE_CREDENTIAL_FILE
            binding.gceCredentials = credFile.text
        }
        new TaskTemplateEngine()
                .setPlaceholder('!' as char)
                .setEnableShortNotation(false)
                .createTemplate(new InputStreamReader(template))
                .make(binding)
                .toString()

    }

    /**
     * @TODO: This is mostly a copy paste from AmazonCloudDriver
     */
    String scriptBashEnv(LaunchConfig cfg) {
        def profile = """\
        export NXF_VER='${cfg.nextflow.version}'
        export NXF_MODE='${cfg.nextflow.mode}'
        export NXF_EXECUTOR='ignite'
        export NXF_CLUSTER_JOIN='cloud:gce:${cfg.clusterName}'
        """
                .stripIndent()

        if (cfg.nextflow.trace)
            profile += "export NXF_TRACE='${cfg.nextflow.trace}'\n"

        if (cfg.nextflow.options)
            profile += "export NXF_OPTS='${cfg.nextflow.options}'\n"

        if (cfg.sharedStorageId && cfg.sharedStorageMount) {
            profile += "export NXF_WORK='${cfg.sharedStorageMount}/${cfg.userName}/work'\n"
            profile += "export NXF_ASSETS='${cfg.sharedStorageMount}/${cfg.userName}/projects'\n"
        }

        if (System.getenv("GOOGLE_APPLICATION_CREDENTIALS")) {
            profile += 'export GOOGLE_APPLICATION_CREDENTIALS='+GCE_CREDENTIAL_FILE+"\n"
        }
        if (System.getenv("NEXFLOW_DOWNLOAD_URL")) {
            profile += 'export NEXFLOW_DOWNLOAD_URL='+System.getenv("NEXFLOW_DOWNLOAD_URL")+"\n"
        }

        return profile
    }
}