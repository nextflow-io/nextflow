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

import static GoogleLifeSciencesHelper.*

import com.google.api.services.lifesciences.v2beta.CloudLifeSciences
import com.google.api.services.lifesciences.v2beta.model.Action
import com.google.api.services.lifesciences.v2beta.model.CancelOperationRequest
import com.google.api.services.lifesciences.v2beta.model.Mount
import com.google.api.services.lifesciences.v2beta.model.Operation
import com.google.api.services.lifesciences.v2beta.model.Pipeline
import com.google.api.services.lifesciences.v2beta.model.Resources
import com.google.api.services.lifesciences.v2beta.model.RunPipelineRequest
import com.google.auth.oauth2.GoogleCredentials
import nextflow.executor.res.AcceleratorResource
import spock.lang.Specification

class GoogleLifeSciencesHelperTest extends Specification {

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
        def action1 = helper.createAction(actionName,imageName,commands,mounts,flags,entryPoint)
        def action2 = helper.createAction(actionName,imageName,commands,mounts)

        then:
        with(action1) {
            getContainerName() == actionName
            getImageUri() == imageName
            getCommands() == commands
            getMounts() == mounts
            getEntrypoint() == entryPoint
            getAlwaysRun()
        }

        with(action2) {
            getContainerName() == actionName
            getImageUri() == imageName
            getCommands() == commands
            getMounts() == mounts
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
        }

        when:
        def resources3 = helper.createResources( new GoogleLifeSciencesSubmitRequest(
                    machineType:type,
                    zone:zone,
                    diskName:diskName,
                    preemptible: false,
                    accelerator: acc,
                    bootDiskSizeGb: 75 ))
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

    def 'should create main action' () {
        given:
        def helper = Spy(GoogleLifeSciencesHelper)
        def req = Mock(GoogleLifeSciencesSubmitRequest)
        def mount = new Mount()

        when:
        def action = helper.createMainAction(req)

        then:
        req.taskName >> 'foo'
        req.containerImage >>'my/image'
        req.mainScript >> 'something.sh'
        req.sharedMount >> mount

        1 * helper.createAction(
                'foo-main',
                'my/image',
                ['bash', '-c', 'something.sh'],
                [mount] )

        action.getContainerName() == 'foo-main'
        action.getImageUri() == 'my/image'
        action.getCommands() == ['bash', '-c', 'something.sh']
        action.getMounts() == [mount]
        action.getEntrypoint() == null
        action.getEnvironment() == [:]
    }

    def 'should create staging action' () {
        given:
        def helper = Spy(GoogleLifeSciencesHelper)
        def req = Mock(GoogleLifeSciencesSubmitRequest)
        def mount = new Mount()

        when:
        def action = helper.createStagingAction(req)

        then:
        req.taskName >> 'bar'
        req.fileCopyImage >> 'alpine'
        req.stagingScript >> 'copy this and that'
        req.sharedMount >> mount

        1 * helper.createAction(
                'bar-stage',
                'alpine',
                ['bash', '-c', 'copy this and that'],
                [mount] )

        action.getContainerName() == 'bar-stage'
        action.getImageUri() == 'alpine'
        action.getCommands() == ['bash', '-c', 'copy this and that']
        action.getMounts() == [mount]
        action.getEntrypoint() == null
        action.getEnvironment() == [:]
        !action.getAlwaysRun()
        !action.getIgnoreExitStatus()
    }

    def 'should create unstaging action' () {
        given:
        def helper = Spy(GoogleLifeSciencesHelper)
        def req = Mock(GoogleLifeSciencesSubmitRequest)
        def mount = new Mount()

        when:
        def action = helper.createUnstagingAction(req)

        then:
        req.taskName >> 'bar'
        req.fileCopyImage >> 'alpine'
        req.unstagingScript >> 'upload command here'
        req.sharedMount >> mount

        1 * helper.createAction(
                'bar-unstage',
                'alpine',
                ['bash', '-c', 'upload command here'],
                [mount],
                [ActionFlags.ALWAYS_RUN, ActionFlags.IGNORE_EXIT_STATUS]
        )

        action.getContainerName() == 'bar-unstage'
        action.getImageUri() == 'alpine'
        action.getCommands() == ['bash', '-c', 'upload command here']
        action.getMounts() == [mount]
        action.getEntrypoint() == null
        action.getEnvironment() == [:]
        action.getAlwaysRun()
        action.getIgnoreExitStatus()
    }

    def 'should submit pipeline request' () {
        given:
        def helper = Spy(GoogleLifeSciencesHelper)
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

    def 'should create gs path' () {
        expect:
        gsPath(PATH).toUriString() == PATH

        where:
        _ | PATH
        _ | 'gs://foo'
        _ | 'gs://foo/bar'      
        _ | 'gs://foo/b a r'
        _ | 'gs://f o o/bar'
        _ | 'gs://f_o_o/bar'
    }
}

