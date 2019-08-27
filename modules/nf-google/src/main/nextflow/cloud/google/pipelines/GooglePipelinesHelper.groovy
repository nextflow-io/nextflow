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

package nextflow.cloud.google.pipelines


import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.services.genomics.v2alpha1.Genomics
import com.google.api.services.genomics.v2alpha1.model.Action
import com.google.api.services.genomics.v2alpha1.model.CancelOperationRequest
import com.google.api.services.genomics.v2alpha1.model.Disk
import com.google.api.services.genomics.v2alpha1.model.Mount
import com.google.api.services.genomics.v2alpha1.model.Operation
import com.google.api.services.genomics.v2alpha1.model.Pipeline
import com.google.api.services.genomics.v2alpha1.model.Resources
import com.google.api.services.genomics.v2alpha1.model.RunPipelineRequest
import com.google.api.services.genomics.v2alpha1.model.ServiceAccount
import com.google.api.services.genomics.v2alpha1.model.VirtualMachine
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Helper class for Google Pipelines.
 *
 * @author Ólafur Haukur Flygenring <olafurh@wuxinextcode.com>
 * @author  Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@SuppressWarnings("GrMethodMayBeStatic")
class GooglePipelinesHelper {

    static final String SCOPE_CLOUD_PLATFORM = "https://www.googleapis.com/auth/cloud-platform"
    static final List<String> ENV_VAR_TO_INCLUDE = ["NXF_DEBUG"]

    Genomics genomicsClient
    GoogleCredential credential
    final String applicationName


    /**
     * As defined by Google pipeline API
     * https://cloud.google.com/genomics/reference/rest/Shared.Types/Flag
     */
    enum ActionFlags {
        FLAG_UNSPECIFIED,
        IGNORE_EXIT_STATUS,
        RUN_IN_BACKGROUND,
        ALWAYS_RUN,
        ENABLE_FUSE,
        PUBLISH_EXPOSED_PORTS,
        DISABLE_IMAGE_PREFETCH,
        DISABLE_STANDARD_ERROR_CAPTURE
    }

    GooglePipelinesHelper(GoogleCredential credential = null, String name = "Nextflow GooglePipelinesExecutor") {
        this.credential = credential
        this.applicationName = name
    }

    static String sanitizeName(String name) {
        name.replaceAll(/[^a-zA-Z0-9\-_]+/, '-').take(63)
    }

    void init() {
        if (!genomicsClient) {

            if (!credential)
                credential = GoogleCredential.applicationDefault

            if (credential.createScopedRequired()) {
                credential = credential.createScoped([SCOPE_CLOUD_PLATFORM])
            }

            genomicsClient = new Genomics.Builder(GoogleNetHttpTransport.newTrustedTransport(), JacksonFactory.defaultInstance, credential)
                    .setApplicationName(applicationName)
                    .build()
        }
    }

    Map<String, String> getEnvironment() {
        System.getenv().findAll { it ->
            ENV_VAR_TO_INCLUDE.contains(it.key)
        }
    }

    Action createAction(String name, String imageUri, List<String> commands, List<Mount> mounts, List<ActionFlags> flags = [], String entrypoint = null) {
        new Action()
                .setName(name)
                .setImageUri(imageUri)
                .setCommands(commands)
                .setMounts(mounts)
                .setFlags(flags.collect { flag -> flag.toString() })
                .setEntrypoint(entrypoint)
                .setEnvironment(getEnvironment())
    }

    Pipeline createPipeline(List<Action> actions, Resources resources) {
        new Pipeline().setActions(actions).setResources(resources)
    }

    Operation submitPipeline(GooglePipelinesSubmitRequest req) {

        final actions = [createStagingAction(req),
                       createMainAction(req),
                       createUnstagingAction(req)]

        final pipeline = createPipeline( actions, createResources(req) )

        runPipeline(pipeline, ["taskName" : req.taskName])
    }

    protected Resources createResources(GooglePipelinesSubmitRequest req) {
        configureResources(
                req.machineType,
                req.project,
                req.zone,
                req.region,
                req.diskName,
                req.diskSizeGb,
                [GooglePipelinesHelper.SCOPE_CLOUD_PLATFORM], req.preemptible)
    }

    protected Action createMainAction(GooglePipelinesSubmitRequest req) {
        createAction(
                "$req.taskName-main",
                req.containerImage,
                ['bash', '-c', req.mainScript],
                [req.sharedMount],
                [ActionFlags.IGNORE_EXIT_STATUS])
    }

    protected Action createStagingAction(GooglePipelinesSubmitRequest req) {
        createAction(
                "$req.taskName-stage",
                req.fileCopyImage,
                ["bash", "-c", req.stagingScript],
                [req.sharedMount],
                [ActionFlags.ALWAYS_RUN, ActionFlags.IGNORE_EXIT_STATUS])
    }

    protected Action createUnstagingAction(GooglePipelinesSubmitRequest req) {
        createAction(
                "$req.taskName-unstage",
                req.fileCopyImage,
                ["bash", "-c", req.unstagingScript],
                [req.sharedMount],
                [ActionFlags.ALWAYS_RUN, ActionFlags.IGNORE_EXIT_STATUS])
    }

    Resources configureResources(String machineType, String projectId, List<String> zone, List<String> region, String diskName, Integer diskSizeGb=null, List<String> scopes = null, boolean preEmptible = false) {

        def disk = new Disk()
        disk.setName(diskName)
        disk.setSizeGb(diskSizeGb)

        def serviceAccount = new ServiceAccount()
        if (scopes)
            serviceAccount.setScopes(scopes)

        def vm = new VirtualMachine()
                .setMachineType(machineType)
                .setDisks([disk])
                .setServiceAccount(serviceAccount)
                .setPreemptible(preEmptible)


        new Resources()
                .setProjectId(projectId)
                .setZones(zone)
                .setRegions(region)
                .setVirtualMachine(vm)
    }

    Mount configureMount(String diskName, String mountPath, boolean readOnly = false) {
        new Mount().setDisk(diskName).setPath(mountPath).setReadOnly(readOnly)
    }

    Operation checkOperationStatus(Operation operation) {
        init()
        try {
            genomicsClient.projects().operations().get(operation.getName()).execute()
        }
        catch( IOException e ) {
            log.warn("Invalid server response fetching operation status: $operation", e)
            return null
        }
    }

    void cancelOperation(Operation operation) {
        init()
        try {
            genomicsClient.projects().operations().cancel(operation.getName(), new CancelOperationRequest()).execute()
        }
        catch( IOException e ) {
            log.warn("Invalid server response cancelling operation: $operation", e)
        }
    }

    Operation runPipeline(Pipeline pipeline, Map<String,String> labels = [:]) {
        init()
        try {
            genomicsClient.pipelines().run(new RunPipelineRequest().setPipeline(pipeline).setLabels(labels)).execute()
        }
        catch( IOException e ) {
            log.warn("Invalid server response running pipeline: $pipeline", e)
            return null
        }
    }
}
