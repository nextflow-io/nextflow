/*
 * Copyright 2019, Google Inc
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

package nextflow.cloud.google.lifesciences

import java.nio.file.Path

import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.http.HttpRequestInitializer
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.services.lifesciences.v2beta.CloudLifeSciences
import com.google.api.services.lifesciences.v2beta.model.Accelerator
import com.google.api.services.lifesciences.v2beta.model.Action
import com.google.api.services.lifesciences.v2beta.model.CancelOperationRequest
import com.google.api.services.lifesciences.v2beta.model.Disk
import com.google.api.services.lifesciences.v2beta.model.Mount
import com.google.api.services.lifesciences.v2beta.model.Network
import com.google.api.services.lifesciences.v2beta.model.Operation
import com.google.api.services.lifesciences.v2beta.model.Pipeline
import com.google.api.services.lifesciences.v2beta.model.Resources
import com.google.api.services.lifesciences.v2beta.model.RunPipelineRequest
import com.google.api.services.lifesciences.v2beta.model.ServiceAccount
import com.google.api.services.lifesciences.v2beta.model.VirtualMachine
import com.google.auth.http.HttpCredentialsAdapter
import com.google.auth.oauth2.GoogleCredentials
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun
import nextflow.util.Escape

/**
 * Helper class for Google Pipelines.
 *
 * @author Ólafur Haukur Flygenring <olafurh@wuxinextcode.com>
 * @author  Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class GoogleLifeSciencesHelper {

    /*
     * Avail location see https://cloud.google.com/life-sciences/docs/concepts/locations
     */
    public static final List<String> DEFAULT_LOCATIONS  = ['us-central1','europe-west2']

    public static final String SSH_DAEMON_NAME = 'ssh-daemon'
    public static final String DEFAULT_APP_NAME = "Nextflow/GLS"
    public static final String SCOPE_CLOUD_PLATFORM = "https://www.googleapis.com/auth/cloud-platform"

    CloudLifeSciences client
    GoogleCredentials credentials
    final String applicationName
    GoogleLifeSciencesConfig config

    /**
     * As defined by Google pipeline API
     * https://cloud.google.com/genomics/reference/rest/Shared.Types/Flag
     */
    static enum ActionFlags {
        FLAG_UNSPECIFIED,
        IGNORE_EXIT_STATUS,
        RUN_IN_BACKGROUND,
        ALWAYS_RUN,
        ENABLE_FUSE,
        PUBLISH_EXPOSED_PORTS,
        DISABLE_IMAGE_PREFETCH,
        DISABLE_STANDARD_ERROR_CAPTURE
    }

    /* only for testing purpose */
    protected GoogleLifeSciencesHelper(GoogleCredentials credential = null, String name = DEFAULT_APP_NAME) {
        this.credentials = credential
        this.applicationName = name
    }

    GoogleLifeSciencesHelper(GoogleLifeSciencesConfig config) {
        this.config = config
        this.applicationName = DEFAULT_APP_NAME
    }

    static String sanitizeName(String name) {
        name.replaceAll(/[^a-zA-Z0-9\-_]+/, '-').take(63)
    }

    GoogleLifeSciencesHelper init() {
        if (!credentials)
            credentials = GoogleCredentials.getApplicationDefault()

        if (credentials.createScopedRequired()) {
            credentials = credentials.createScoped([SCOPE_CLOUD_PLATFORM])
        }

        HttpRequestInitializer requestInitializer = new HttpCredentialsAdapter(credentials);
        client = new CloudLifeSciences.Builder(GoogleNetHttpTransport.newTrustedTransport(), JacksonFactory.defaultInstance, requestInitializer)
                .setApplicationName(applicationName)
                .build()

        return this
    }

    Map<String, String> getEnvironment() {
        final result = new HashMap(1)
        if( config.debugMode )
            result.put('NXF_DEBUG', config.debugMode.toString())
        return result
    }

    Action createAction(String name, String imageUri, List<String> commands, List<Mount> mounts, List<ActionFlags> flags = [], String entrypoint = null) {
        final action = new Action()
                .setContainerName(name)
                .setImageUri(imageUri)
                .setCommands(commands)
                .setMounts(mounts)
                .setEntrypoint(entrypoint)
                .setEnvironment(getEnvironment())

        setFlags(action,flags)
        return action
    }

    protected Action setFlags(Action action, List<ActionFlags> flags) {
        if( ActionFlags.IGNORE_EXIT_STATUS in flags )
            action.setIgnoreExitStatus(true)
        if( ActionFlags.RUN_IN_BACKGROUND in flags )
            action.setRunInBackground(true)
        if( ActionFlags.ALWAYS_RUN in flags)
            action.setAlwaysRun(true)
        if( ActionFlags.ENABLE_FUSE )
            action.setEnableFuse(true)
        if( ActionFlags.PUBLISH_EXPOSED_PORTS in flags )
            action.setPublishExposedPorts(true)
        if( ActionFlags.DISABLE_IMAGE_PREFETCH in flags )
            action.setDisableImagePrefetch(true)
        if( ActionFlags.DISABLE_STANDARD_ERROR_CAPTURE )
            action.setDisableStandardErrorCapture(true)

        return action
    }

    Pipeline createPipeline(List<Action> actions, Resources resources) {
        new Pipeline().setActions(actions).setResources(resources)
    }

    Operation submitPipeline(GoogleLifeSciencesSubmitRequest req) {

        final actions = new ArrayList(5)
        if( config.sshDaemon ) {
            actions.add(createSshDaemonAction(req))
        }
        actions.add(createStagingAction(req))
        actions.add(createMainAction(req))
        actions.add(createUnstagingAction(req))

        final pipeline = createPipeline( actions, createResources(req) )

        runPipeline(req.project, req.location, pipeline, ["taskName" : req.taskName])
    }

    protected Resources createResources(GoogleLifeSciencesSubmitRequest req) {
        def disk = new Disk()
        disk.setName(req.diskName)
        disk.setSizeGb(req.diskSizeGb)

        def serviceAccount = new ServiceAccount().setScopes( [SCOPE_CLOUD_PLATFORM] )

        def vm = new VirtualMachine()
                .setMachineType(req.machineType)
                .setDisks([disk])
                .setServiceAccount(serviceAccount)
                .setPreemptible(req.preemptible)

        if( req.usePrivateAddress ) {
            vm.setNetwork( new Network().setUsePrivateAddress(true) )
        }

        if( req.bootDiskSizeGb ) {
            vm.setBootDiskSizeGb(req.bootDiskSizeGb)
        }
        
        if( req.cpuPlatform ) {
            vm.setCpuPlatform(req.cpuPlatform)
        }

        if( req.accelerator ) {
            final acc = new Accelerator().setType(req.accelerator.type).setCount(req.accelerator.request)
            final list = new ArrayList(1)
            list.add(acc)
            vm.setAccelerators(list)
        }

        final result = new Resources()
                .setZones(req.zone)
                .setRegions(req.region)
                .setVirtualMachine(vm)

        log.trace "[GLS] task=$req.taskName; VM resources=$result"
        return result
    }

    protected Action createMainAction(GoogleLifeSciencesSubmitRequest req) {
        // flag pipefail is required otherwise the command exit status is not returned
        List<String> cmd = ['-o', 'pipefail', '-c', getMainScript(req.workDir)]
        if( !req.entryPoint )
            cmd.add(0, 'bash')

        createAction(
                "$req.taskName-main",
                req.containerImage,
                cmd,
                [req.sharedMount],
                Collections.<ActionFlags>emptyList(),
                req.entryPoint)
    }

    protected Action createStagingAction(GoogleLifeSciencesSubmitRequest req) {
        createAction(
                "$req.taskName-stage",
                config.copyImage,
                ["bash", "-c", getStagingScript(req.workDir)],
                [req.sharedMount] )
    }

    protected Action createUnstagingAction(GoogleLifeSciencesSubmitRequest req) {
        createAction(
                "$req.taskName-unstage",
                config.copyImage,
                ["bash", "-c", getUnstagingScript(req.workDir)],
                [req.sharedMount],
                [ActionFlags.ALWAYS_RUN, ActionFlags.IGNORE_EXIT_STATUS])
    }

    protected Action createSshDaemonAction(GoogleLifeSciencesSubmitRequest req) {
        new Action()
            .setContainerName(SSH_DAEMON_NAME)
            .setImageUri(config.sshImage)
            .setEntrypoint('ssh-server')
            .setMounts([req.sharedMount])
            .setPortMappings(['22':22])
            .setAlwaysRun(true)
            .setRunInBackground(true)
            .setIgnoreExitStatus(true)
    }

    Operation checkOperationStatus(Operation operation) {
        try {
            client.projects().locations().operations().get(operation.getName()).execute()
        }
        catch( IOException e ) {
            log.warn("Invalid server response fetching operation status: ${operation.getName()}", e)
            return null
        }
    }

    void cancelOperation(Operation operation) {
        try {
            client.projects().locations().operations().cancel(operation.getName(), new CancelOperationRequest()).execute()
        }
        catch( IOException e ) {
            log.debug("Invalid server response cancelling operation: ${operation.getName()} | ${e.message}")
        }
    }

    Operation runPipeline(String project, String location, Pipeline pipeline, Map<String,String> labels = [:]) {
        final parent = "projects/$project/locations/$location"
        try {
            client
                    .projects()
                    .locations()
                    .pipelines()
                    .run(parent,
                            new RunPipelineRequest().setPipeline(pipeline) .setLabels(labels) )
                    .execute()
        }
        catch( IOException e ) {
            log.warn("Invalid server response running pipeline: $pipeline", e)
            return null
        }
    }

    static String getLocalTaskDir(Path workDir) { Escape.path(workDir) }

    static String getRemoteTaskDir(Path workDir) { Escape.uriPath(workDir) }

    String getStagingScript(Path workDir) {
        final localTaskDir = getLocalTaskDir(workDir)
        final remoteTaskDir = getRemoteTaskDir(workDir)

        String result = 'set -x; '
        result += '{ '
        result += "cd ${localTaskDir}; "
        result += "gsutil -m -q cp $remoteTaskDir/${TaskRun.CMD_RUN} .; "
        result += "bash ${TaskRun.CMD_RUN} nxf_stage; "
        result += '[[ $NXF_DEBUG -gt 0 ]] && ls -lah $PWD || true; '
        result += "} &> $localTaskDir/${TaskRun.CMD_LOG}"
        return result
    }

    String getMainScript(Path workDir) {
        final localTaskDir = getLocalTaskDir(workDir)
        return "{ cd $localTaskDir; bash ${TaskRun.CMD_RUN}; } >> $localTaskDir/${TaskRun.CMD_LOG} 2>&1"
    }

    String getUnstagingScript(Path workDir) {
        final localTaskDir = getLocalTaskDir(workDir)
        final remoteTaskDir = getRemoteTaskDir(workDir)
        def result = 'set -x; '
        result += "trap 'err=\$?; exec 1>&2; gsutil -m -q cp -R $localTaskDir/${TaskRun.CMD_LOG} ${remoteTaskDir}/${TaskRun.CMD_LOG} || true; [[ \$err -gt 0 || \$GOOGLE_LAST_EXIT_STATUS -gt 0 || \$NXF_DEBUG -gt 0 ]] && { ls -lah $localTaskDir || true; gsutil -m -q cp -R /google/ ${remoteTaskDir}; } || rm -rf $localTaskDir; exit \$err' EXIT; "
        result += "{ cd $localTaskDir; bash ${TaskRun.CMD_RUN} nxf_unstage; } >> $localTaskDir/${TaskRun.CMD_LOG} 2>&1"
        return result
    }

    void checkValidLocation() {
        if( config.location in DEFAULT_LOCATIONS ) {
            // that's fine, just return
            return
        }
        // check the current list of location to make sure the specified one really exists
        final availLocations = client
                .projects()
                .locations()
                .list("projects/$config.project")
                .execute()
                .getLocations()
                *.getLocationId()
        if( config.location in availLocations ) {
            return
        }
        // show a warning message
        log.warn "The specified Google Life Sciences location is not available: \"$config.location\" -- Please choose open of the following ${availLocations.join(',')}"
    }
}
