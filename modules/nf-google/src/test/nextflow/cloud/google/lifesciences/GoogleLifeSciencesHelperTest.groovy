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
package nextflow.cloud.google.lifesciences

import static nextflow.cloud.google.lifesciences.GoogleLifeSciencesHelper.*

import com.google.api.services.lifesciences.v2beta.CloudLifeSciences
import com.google.api.services.lifesciences.v2beta.model.Action
import com.google.api.services.lifesciences.v2beta.model.CancelOperationRequest
import com.google.api.services.lifesciences.v2beta.model.Mount
import com.google.api.services.lifesciences.v2beta.model.Operation
import com.google.api.services.lifesciences.v2beta.model.Pipeline
import com.google.api.services.lifesciences.v2beta.model.Resources
import com.google.api.services.lifesciences.v2beta.model.RunPipelineRequest
import com.google.auth.oauth2.GoogleCredentials
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import nextflow.cloud.google.GoogleSpecification
import nextflow.executor.res.AcceleratorResource
import spock.lang.Unroll

class GoogleLifeSciencesHelperTest extends GoogleSpecification {

    def mockClient = Mock(CloudLifeSciences) {
        projects() >> Mock(CloudLifeSciences.Projects) {
            locations() >> Mock(CloudLifeSciences.Projects.Locations) {
                operations() >> Mock(CloudLifeSciences.Projects.Locations.Operations) {
                    get(_) >> Mock(CloudLifeSciences.Projects.Locations.Operations.Get)
                    cancel(_,_) >> Mock(CloudLifeSciences.Projects.Locations.Operations.Cancel)
                }
                pipelines() >> Mock(CloudLifeSciences.Projects.Locations.Pipelines) {
                    run(_,_) >> Mock(CloudLifeSciences.Projects.Locations.Pipelines.Run)
                }
            }
        }
    }

    def 'should sanitize names' () {
        given:
        def longName = "x" * 65
        def spaces = "this is a name"
        def helper = new GoogleLifeSciencesHelper()

        when:
        def longSanitized = helper.sanitizeName(longName)
        def spacesSanitized = helper.sanitizeName(spaces)

        then:
        longSanitized == "x" * 63
        spacesSanitized == "this-is-a-name"
    }

    def 'should initialize the genomics client'() {
        given:
        def stubCredential = Stub(GoogleCredentials)
        def helper = new GoogleLifeSciencesHelper(stubCredential,"testName")

        when:
        helper.init()

        then:
        helper.applicationName == "testName"
        helper.client
    }

    def 'should construct a pipeline action'() {
        given:
        def actionName = "actionName"
        def imageName = "imageName"
        def commands = ["command1","command2"]
        def mounts = []
        def flags = [GoogleLifeSciencesHelper.ActionFlags.ALWAYS_RUN]
        def entryPoint = "entryPoint"
        def helper = new GoogleLifeSciencesHelper()

        when:
        helper.config = Mock(GoogleLifeSciencesConfig)
        def action1 = helper.createAction(actionName,imageName,commands,mounts,flags,entryPoint)
        then:
        with(action1) {
            getContainerName() == actionName
            getImageUri() == imageName
            getCommands() == commands
            getMounts() == mounts
            getEntrypoint() == entryPoint
            getAlwaysRun()
            getEnvironment() == [:]
        }

        when:
        helper.config = Mock(GoogleLifeSciencesConfig) { getDebugMode() >> 2 }
        def action2 = helper.createAction(actionName,imageName,commands,mounts)
        then:

        with(action2) {
            getContainerName() == actionName
            getImageUri() == imageName
            getCommands() == commands
            getMounts() == mounts
            getEnvironment() == [NXF_DEBUG: '2']
            !getEntrypoint()
            !getAlwaysRun()
        }
    }

    def 'should construct a pipeline'() {
        given:
        def actions = [new Action(),new Action()]
        def resources = new Resources().setRegions(['foo'])
        def helper = new GoogleLifeSciencesHelper()

        when:
        def pipe = helper.createPipeline(actions,resources)

        then:
        pipe.getActions().size() == 2
        pipe.getResources().getRegions() == ['foo']

    }

    def 'should configure resources correctly'() {
        given:
        def SCOPES = ["https://www.googleapis.com/auth/cloud-platform"]
        def type = "testType"
        def zone = ["testZone1","testZone2"]
        def region = ["testRegion1","testRegion2"]
        def diskName = "testDisk"
        def preEmptible = true
        def acc = new AcceleratorResource(request: 4, type: 'nvidia-tesla-k80')
        def helper = new GoogleLifeSciencesHelper()

        when:
        def resources1 = helper.createResources(new GoogleLifeSciencesSubmitRequest(
                machineType:type,
                zone:zone,
                diskName: diskName,
                diskSizeGb: 100,
                preemptible: true))
        then:
        with(resources1) {
            getVirtualMachine().getMachineType() == type
            getZones() == zone
            getRegions() == null
            getVirtualMachine().getDisks().get(0).getName() == diskName
            getVirtualMachine().getDisks().get(0).getSizeGb() == 100
            getVirtualMachine().getServiceAccount().getScopes() == SCOPES
            getVirtualMachine().getPreemptible() == preEmptible
            !getVirtualMachine().getAccelerators()
            !getVirtualMachine().getNetwork()?.getUsePrivateAddress()
        }

        when:
        def resources2 = helper.createResources(new GoogleLifeSciencesSubmitRequest(
                machineType:type,
                region:region,
                diskName:diskName,
                diskSizeGb: 200,
                preemptible: true))
        then:
        with(resources2) {
            getVirtualMachine().getMachineType() == type
            getZones() == null
            getRegions() == region
            getVirtualMachine().getDisks().get(0).getName() == diskName
            getVirtualMachine().getDisks().get(0).getSizeGb() == 200
            getVirtualMachine().getServiceAccount().getScopes() == SCOPES
            getVirtualMachine().getPreemptible() == preEmptible
            !getVirtualMachine().getAccelerators()
            !getVirtualMachine().getNetwork()?.getUsePrivateAddress()
        }

        when:
        def resources3 = helper.createResources( new GoogleLifeSciencesSubmitRequest(
                    machineType:type,
                    zone:zone,
                    diskName:diskName,
                    preemptible: false,
                    accelerator: acc,
                    bootDiskSizeGb: 75,
                    cpuPlatform: 'Intel Skylake',
                    usePrivateAddress: true ))
        then:
        with(resources3) {
            getVirtualMachine().getMachineType() == type
            getZones() == zone
            getVirtualMachine().getDisks().get(0).getName() == diskName
            getVirtualMachine().getServiceAccount().getScopes() == SCOPES
            !getVirtualMachine().getDisks().get(0).getSizeGb()
            !getVirtualMachine().getPreemptible()
            getVirtualMachine().getAccelerators().size()==1
            getVirtualMachine().getAccelerators()[0].getCount()==4
            getVirtualMachine().getAccelerators()[0].getType()=='nvidia-tesla-k80'
            getVirtualMachine().getBootDiskSizeGb() == 75
            getVirtualMachine().getCpuPlatform() == 'Intel Skylake'
            getVirtualMachine().getNetwork().getUsePrivateAddress()
        }
    }


    def 'should check on the status of an operation'() {
        given:
        def op = new Operation().setName("testOperation")
        def helper = new GoogleLifeSciencesHelper()
        helper.client = mockClient

        when:
        def opStatus = helper.checkOperationStatus(op)

        then:
        1 * mockClient.projects().locations().operations().get(op.getName()).execute() >> new Operation().setName("status").setDone(true)
        opStatus.getName() == "status"
        opStatus.getDone()

    }

    def 'should cancel an operation'() {
        given:
        def op = new Operation().setName("testOperation").setDone(false)
        def helper = new GoogleLifeSciencesHelper()
        helper.client = mockClient

        when:
        helper.cancelOperation(op)

        then:
        1 * mockClient.projects().locations().operations().cancel(op.getName(),new CancelOperationRequest()).execute() >> {}
    }

    def 'should run a pipeline'() {
        given:
        def pipe = new Pipeline()
        def request = new RunPipelineRequest().setPipeline(pipe)

        def helper = new GoogleLifeSciencesHelper()
        helper.client = mockClient

        when:
        def runOp = helper.runPipeline('PRJ-X', 'LOC-Y', pipe)

        then:
        1 * mockClient.projects().locations().pipelines().run('projects/PRJ-X/locations/LOC-Y', request).execute() >> new Operation().setName("runOp")

        runOp.getName() == "runOp"
    }

    def 'should create main action with entry' () {
        given:
        def mount = new Mount()
        def workDir = mockGsPath('gs://foo/work')
        def helper = Spy(GoogleLifeSciencesHelper)
        helper.config = Mock(GoogleLifeSciencesConfig)
        and:
        def req = Mock(GoogleLifeSciencesSubmitRequest) {
            getTaskName() >> 'foo'
            getContainerImage() >>'my/image'
            getSharedMount() >> mount
            getWorkDir() >> workDir
            getEntryPoint() >> '/bin/sh'
        }

        when:
        def action = helper.createMainAction(req)

        then:
        1 * helper.getMainScript(workDir) >> 'main.sh'
        and:
        1 * helper.createAction(
                'foo-main',
                'my/image',
                ['-o','pipefail','-c', 'main.sh'],
                [mount],
                [],
                '/bin/sh' )

        and:
        action.getContainerName() == 'foo-main'
        action.getImageUri() == 'my/image'
        action.getCommands() == ['-o', 'pipefail','-c', 'main.sh']
        action.getMounts() == [mount]
        action.getEntrypoint() == '/bin/sh'
        action.getEnvironment() == [:]
    }

    def 'should create action with no entry' () {
        given:
        def mount = new Mount()
        def workDir = mockGsPath('gs://foo/work')
        def helper = Spy(GoogleLifeSciencesHelper)
        helper.config = Mock(GoogleLifeSciencesConfig)
        and:
        def req = Mock(GoogleLifeSciencesSubmitRequest) {
            getTaskName() >> 'foo'
            getContainerImage() >>'my/image'
            getSharedMount() >> mount
            getWorkDir() >> workDir
            getEntryPoint() >> null
        }

        when:
        def action = helper.createMainAction(req)

        then:
        1 * helper.getMainScript(workDir) >> 'main.sh'
        and:
        1 * helper.createAction(
                'foo-main',
                'my/image',
                ['bash','-o','pipefail','-c', 'main.sh'],
                [mount],
                [],
                null )

        and:
        action.getContainerName() == 'foo-main'
        action.getImageUri() == 'my/image'
        action.getCommands() == ['bash', '-o', 'pipefail','-c', 'main.sh']
        action.getMounts() == [mount]
        action.getEnvironment() == [:]
        action.getEntrypoint() == null
    }


    def 'should create staging action' () {
        given:
        def workDir = mockGsPath('gs://foo/work')
        def mount = new Mount()
        and:
        def helper = Spy(GoogleLifeSciencesHelper)
        helper.config = Mock(GoogleLifeSciencesConfig) { getCopyImage() >> 'alpine' }
        and:
        def req = Mock(GoogleLifeSciencesSubmitRequest) {
            getTaskName() >> 'bar'
            getSharedMount() >> mount
            getWorkDir() >> workDir
        }

        when:
        def action = helper.createStagingAction(req)
        then:
        1 * helper.getStagingScript(workDir) >> 'stage.sh'
        and:
        1 * helper.createAction(
                'bar-stage',
                'alpine',
                ['bash', '-c', 'stage.sh'],
                [mount] )

        action.getContainerName() == 'bar-stage'
        action.getImageUri() == 'alpine'
        action.getCommands() == ['bash', '-c', 'stage.sh']
        action.getMounts() == [mount]
        action.getEntrypoint() == null
        action.getEnvironment() == [:]
        !action.getAlwaysRun()
        !action.getIgnoreExitStatus()
    }

    def 'should create unstaging action' () {
        given:
        def workDir = mockGsPath('gs://foo/work')
        def mount = new Mount()
        and:
        def helper = Spy(GoogleLifeSciencesHelper)
        helper.config = Mock(GoogleLifeSciencesConfig) { getCopyImage() >> 'alpine' }
        and:
        def req = Mock(GoogleLifeSciencesSubmitRequest) {
            getTaskName() >> 'bar'
            getSharedMount() >> mount
            getWorkDir() >> workDir
        }

        when:
        def action = helper.createUnstagingAction(req)

        then:
        1 * helper.getUnstagingScript(workDir) >> 'unstage.sh'
        and:
        1 * helper.createAction(
                'bar-unstage',
                'alpine',
                ['bash', '-c', 'unstage.sh'],
                [mount],
                [ActionFlags.ALWAYS_RUN, ActionFlags.IGNORE_EXIT_STATUS]
        )

        action.getContainerName() == 'bar-unstage'
        action.getImageUri() == 'alpine'
        action.getCommands() == ['bash', '-c', 'unstage.sh']
        action.getMounts() == [mount]
        action.getEntrypoint() == null
        action.getEnvironment() == [:]
        action.getAlwaysRun()
        action.getIgnoreExitStatus()
    }

    def 'should submit pipeline request' () {
        given:
        def helper = Spy(GoogleLifeSciencesHelper)
        helper.config = Mock(GoogleLifeSciencesConfig) { getCopyImage() >> 'alpine' }
        def stage = GroovyMock(Action)
        def unstage = GroovyMock(Action)
        def main = GroovyMock(Action)
        def res = GroovyMock(Resources)
        def operation = GroovyMock(Operation)
        and:

        def pipeline = new Pipeline()
        def req = new GoogleLifeSciencesSubmitRequest(location: 'LOC-1', project: 'PRJ-X', taskName: 'foo')

        when:
        def result = helper.submitPipeline(req)

        then:
        req.taskName >> 'foo'
        1 * helper.createStagingAction(req) >> stage
        1 * helper.createUnstagingAction(req) >> unstage
        1 * helper.createMainAction(req) >> main
        1 * helper.createResources(req) >> res
        1 * helper.createPipeline([stage, main, unstage], res) >> pipeline
        1 * helper.runPipeline('PRJ-X','LOC-1', pipeline, [taskName: 'foo']) >> operation
        and:
        result == operation
    }

    def 'should set action flags' () {
        given:
        def helper = new GoogleLifeSciencesHelper()

        expect:
        helper
            .setFlags(new Action(), [ActionFlags.IGNORE_EXIT_STATUS])
            .getIgnoreExitStatus()

        helper
                .setFlags(new Action(), [ActionFlags.RUN_IN_BACKGROUND])
                .getRunInBackground()

        helper
                .setFlags(new Action(), [ActionFlags.ALWAYS_RUN])
                .getAlwaysRun()


        helper
                .setFlags(new Action(), [ActionFlags.ENABLE_FUSE])
                .getEnableFuse()

        helper
                .setFlags(new Action(), [ActionFlags.PUBLISH_EXPOSED_PORTS])
                .getPublishExposedPorts()

        helper
                .setFlags(new Action(), [ActionFlags.DISABLE_IMAGE_PREFETCH])
                .getDisableImagePrefetch()

        helper
                .setFlags(new Action(), [ActionFlags.DISABLE_STANDARD_ERROR_CAPTURE])
                .getDisableStandardErrorCapture()
    }

    def 'should create wrapper scripts' () {
        given:
        def dir = CloudStorageFileSystem.forBucket("my-bucket").getPath('/work/dir')
        def helper = new GoogleLifeSciencesHelper()

        when:
        def stage = helper.getStagingScript(dir)
        then:
        stage ==
                'set -x; { cd /work/dir; gsutil -m -q cp gs://my-bucket/work/dir/.command.run .; bash .command.run nxf_stage; [[ $NXF_DEBUG -gt 0 ]] && ls -lah $PWD || true; } &> /work/dir/.command.log'
        when:
        def main = helper.getMainScript(dir)
        then:
        main ==
                '{ cd /work/dir; bash .command.run; } >> /work/dir/.command.log 2>&1'

        when:
        def unstage = helper.getUnstagingScript(dir)
        then:
        unstage ==
                'set -x; trap \'err=$?; exec 1>&2; gsutil -m -q cp -R /work/dir/.command.log gs://my-bucket/work/dir/.command.log || true; [[ $err -gt 0 || $GOOGLE_LAST_EXIT_STATUS -gt 0 || $NXF_DEBUG -gt 0 ]] && { ls -lah /work/dir || true; gsutil -m -q cp -R /google/ gs://my-bucket/work/dir; } || rm -rf /work/dir; exit $err\' EXIT; { cd /work/dir; bash .command.run nxf_unstage; } >> /work/dir/.command.log 2>&1'
    }

    @Unroll
    def 'should get local dir' () {
        expect:
        getLocalTaskDir(mockGsPath(PATH)) == EXPECTED
        where:
        PATH                        | EXPECTED
        'gs://my-bucket/work/dir'   | '/work/dir'
        'gs://my-bucket/work dir'   | '/work\\ dir'

    }

    @Unroll
    def 'should get remote dir' () {
        expect:
        getRemoteTaskDir(mockGsPath(PATH)) == EXPECTED
        where:
        PATH                        | EXPECTED
        'gs://my-bucket/work/dir'   | 'gs://my-bucket/work/dir'
        'gs://my-bucket/work dir'   | 'gs://my-bucket/work\\ dir'
        'gs://my bucket/work/dir'   | 'gs://my\\ bucket/work/dir'
    }

}

