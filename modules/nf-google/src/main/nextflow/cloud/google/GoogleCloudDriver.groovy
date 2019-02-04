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

import java.nio.file.Files
import java.nio.file.Path

import com.google.api.services.compute.Compute
import com.google.api.services.compute.model.AccessConfig
import com.google.api.services.compute.model.Instance
import com.google.api.services.compute.model.InstancesSetLabelsRequest
import com.google.api.services.compute.model.MachineType
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
import nextflow.file.FileHelper
import nextflow.processor.TaskTemplateEngine
import nextflow.util.MemoryUnit
import nextflow.util.ServiceName
import static nextflow.cloud.CloudConst.TAG_CLUSTER_NAME
import static nextflow.cloud.CloudConst.TAG_CLUSTER_ROLE
/**
 * Cloud driver implementation for Google Compute Engine.
 *
 * @author Vilmundur Pálmason <vilmundur@wuxinextcode.com>
 */
@Slf4j
@CompileStatic
@ServiceName('google')
class GoogleCloudDriver implements CloudDriver {

    /**
     * The GCE zone eg. {@code us-central1-f}. If it's not specified the current region is retrieved from
     * the GCE instance metadata
     */
    static final String TERMINATION_FILENAME = "/tmp/shutdown.begin"
    static final String GCE_CREDENTIAL_FILE = '$HOME/.nextflow/gce_credentials.json'

    static long OPS_WAIT_TIMEOUT_MS = 300000 //5 minutes
    static long POLL_WAIT = 1000

    private GceApiHelper helper

    /**
     * Only use for testing
     */
    GoogleCloudDriver(GceApiHelper helper) {
        this.helper = helper
    }

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
    @Deprecated
    GoogleCloudDriver(Map config) {
        String zone = config.zone ?: Global.config?.navigate('google.zone')
        String project = config.project ?: Global.config?.navigate('google.project')
        this.helper = new GceApiHelper(project, zone)
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

        try {
            if (!helper.lookupImage(config.imageId))
                throw new AbortOperationException("Unknown GCE image ID: ${config.imageId}")
        }
        catch( Exception e ) {
            throw new AbortOperationException("Invalid GCE image ID: ${config.imageId}", e)
        }

        try {
            if (!helper.lookupMachineType(config.instanceType))
                throw new AbortOperationException("Unknown GCE machine type: ${config.instanceType}")
        }
        catch ( Exception e ) {
            throw new AbortOperationException("Invalid GCE machine type: ${config.instanceType}", e)
        }

        String validationError = helper.validateLabelValue(config.clusterName)
        if (validationError != null) {
            throw new AbortOperationException("Invalid cluster name '" + config.clusterName + "': " + validationError)
        }

        checkUnsupportedSettings(config)
    }

    protected void checkUnsupportedSettings(LaunchConfig config){
        List<String> UNSUPPORTED = [
                'bootStorageSize',
                'instanceRole',
                'instanceStorageMount',
                'instanceStorageDevice',
                'securityGroup',
                'sharedStorageId',
                'sharedStorageMount',
                'subnetId',
                'spotPrice'
        ]

        def invalid = config.getAttributeNames().intersect((Collection)UNSUPPORTED)
        if( invalid )
            throw new AbortOperationException("The following cloud settings are not supported by GCE driver: ${invalid.join(', ')}")
    }

    @Override
    List<String> launchInstances(int instanceCount, LaunchConfig config) {
        List<String> result = []
        List<Operation> ops = []
        instanceCount.times {
            Instance inst = new Instance()
            inst.setName(helper.randomName(config.getClusterName() + "-"))
            inst.setMachineType(helper.instanceType(config.getInstanceType()))
            if(config.preemptible) {
                inst.setScheduling(helper.createScheduling(config.preemptible))
                helper.setShutdownScript(inst, gceShutdownScript())
            }

            inst.setDisks([helper.createBootDisk(inst.getName(), config.getImageId())])
            inst.setNetworkInterfaces([helper.createNetworkInterface()])
            helper.setStartupScript(inst, gceStartupScript(config))

            def insert = helper.compute.instances().insert(helper.project, helper.zone, inst)

            result << inst.getName()
            ops << insert.execute()
        }
        helper.blockUntilComplete(ops, OPS_WAIT_TIMEOUT_MS)
        return result
    }

    @Override
    void waitInstanceStatus(Collection<String> instanceIds, CloudInstanceStatus status) {
        def instanceStatusList

        // check google cloud instance status
        // https://cloud.google.com/compute/docs/instances/checking-instance-status
        switch (status) {
            case CloudInstanceStatus.STARTED:
                instanceStatusList = ['RUNNING']
                waitStatus(instanceIds, instanceStatusList)
                break

            // ideally 'READY' status should detected when the instance
            // has terminated the boot/initialisation phase, however the GCE
            // does not provide an API for that, therefore it's the same
            // as 'STARTED' status (see above)
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
            List<Instance> instances = helper.getInstanceList(filter)
            if (instances) {
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
        throw new UnsupportedOperationException("Spot prices feature is not supported by the GCE driver")
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
        helper.getInstanceList(filter)?.each { inst -> callback.call(toNextflow(inst)) }
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
        Path shutdownFile = FileHelper.asPath(TERMINATION_FILENAME)

        if(Files.exists(shutdownFile)) {
            shutdownFile.toFile().text
        }
        else {
            return null
        }
    }

    @Override
    CloudInstanceType describeInstanceType(String instanceType) {
        MachineType type = helper.lookupMachineType(instanceType)
        if(!type) {
            return null
        } else {
            new CloudInstanceType(type.getName(),type.getGuestCpus(),new MemoryUnit("${type.getMemoryMb()} MB"),new MemoryUnit("${type.getImageSpaceGb()} GB"),type.getMaximumPersistentDisks())
        }
    }

    static String tagToFilterExpression(String k, v) {
        '(labels.' + k.toLowerCase() + ' = "' + (v ?: '*') + '")'
    }

    static String instanceIdToFilterExpression(String instanceId) {
        '(name = "' + instanceId + '")'
    }

    def toNextflow(Instance instance) {
        NetworkInterface networkInterface
        AccessConfig accessConfig
        def labels = instance.getLabels() ?: {}

        if (instance.getNetworkInterfaces() != null && !instance.getNetworkInterfaces().isEmpty()) {
            networkInterface = instance.getNetworkInterfaces()[0]
            if (networkInterface.getAccessConfigs() != null && !networkInterface.getAccessConfigs()?.isEmpty()) {
                accessConfig = networkInterface.getAccessConfigs()[0]
            }
        }

        new CloudInstance(
                id: instance.getName(),
                privateIpAddress: networkInterface?.getNetworkIP(),
                privateDnsName: helper.instanceIdToPrivateDNS(instance.getName()),
                publicIpAddress: accessConfig?.getNatIP(),
                publicDnsName: accessConfig?.getPublicPtrDomainName() ?: helper.publicIpToDns(accessConfig?.getNatIP()),
                state: instance.getStatus(),
                role: labels[TAG_CLUSTER_ROLE.toLowerCase()],
                clusterName: labels[TAG_CLUSTER_NAME.toLowerCase()]
        )
    }

    @PackageScope
    String gceShutdownScript() {
        """
            #!/bin/bash            
            date >> $TERMINATION_FILENAME                                               
        """.stripIndent().leftTrim()
    }

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
        binding.gceCredentialsFile = GCE_CREDENTIAL_FILE
        binding.gceCredentials = helper.getCredentialsFile()

        new TaskTemplateEngine()
                .setPlaceholder('!' as char)
                .setEnableShortNotation(false)
                .createTemplate(new InputStreamReader(template))
                .make(binding)
                .toString()

    }

    String scriptBashEnv(LaunchConfig cfg) {
        def profile = """\
        export NXF_VER='${cfg.nextflow.version}'
        export NXF_MODE='google'
        export NXF_EXECUTOR='ignite'
        export NXF_CLUSTER_JOIN='cloud:google:${cfg.clusterName}'
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

        profile += "export GOOGLE_APPLICATION_CREDENTIALS=${GCE_CREDENTIAL_FILE}\n"

        if (System.getenv("NEXTFLOW_DOWNLOAD_URL")) {
            profile += 'export NEXTFLOW_DOWNLOAD_URL='+System.getenv("NEXTFLOW_DOWNLOAD_URL")+"\n"
        }

        return profile
    }
}