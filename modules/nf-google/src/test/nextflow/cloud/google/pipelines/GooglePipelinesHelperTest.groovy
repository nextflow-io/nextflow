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

import spock.lang.Specification

import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.google.api.services.genomics.v2alpha1.Genomics
import com.google.api.services.genomics.v2alpha1.model.Action
import com.google.api.services.genomics.v2alpha1.model.CancelOperationRequest
import com.google.api.services.genomics.v2alpha1.model.Mount
import com.google.api.services.genomics.v2alpha1.model.Operation
import com.google.api.services.genomics.v2alpha1.model.Pipeline
import com.google.api.services.genomics.v2alpha1.model.Resources
import com.google.api.services.genomics.v2alpha1.model.RunPipelineRequest
import static nextflow.cloud.google.pipelines.GooglePipelinesHelper.ActionFlags.ALWAYS_RUN
import static nextflow.cloud.google.pipelines.GooglePipelinesHelper.ActionFlags.IGNORE_EXIT_STATUS

class GooglePipelinesHelperTest extends Specification {

    def mockClient = Mock(Genomics) {
        projects() >> Mock(Genomics.Projects) {
            operations() >> Mock(Genomics.Projects.Operations) {
                get(_) >> Mock(Genomics.Projects.Operations.Get)
                cancel(_,_) >> Mock(Genomics.Projects.Operations.Cancel)
            }
        }
        pipelines() >> Mock(Genomics.Pipelines) {
            run(_) >> Mock(Genomics.Pipelines.Run)
        }
    }

    def 'should sanitize names' () {
        given:
        def longName = "x" * 65
        def spaces = "this is a name"
        def helper = new GooglePipelinesHelper()

        when:
        def longSanitized = helper.sanitizeName(longName)
        def spacesSanitized = helper.sanitizeName(spaces)

        then:
        longSanitized == "x" * 63
        spacesSanitized == "this-is-a-name"
    }

    def 'should initialize the genomics client'() {
        given:
        def stubCredential = Stub(GoogleCredential)
        def helper = new GooglePipelinesHelper(stubCredential,"testName")

        when:
        helper.init()

        then:
        helper.applicationName == "testName"
        helper.genomicsClient
    }

    def 'should construct a pipeline action'() {
        given:
        def actionName = "actionName"
        def imageName = "imageName"
        def commands = ["command1","command2"]
        def mounts = []
        def flags = [GooglePipelinesHelper.ActionFlags.ALWAYS_RUN]
        def entryPoint = "entryPoint"
        def helper = new GooglePipelinesHelper()


        when:
        def action1 = helper.createAction(actionName,imageName,commands,mounts,flags,entryPoint)
        def action2 = helper.createAction(actionName,imageName,commands,mounts)

        then:
        with(action1) {
            getName() == actionName
            getImageUri() == imageName
            getCommands() == commands
            getMounts() == mounts
            getFlags() == flags.collect{it.toString()}
            getEntrypoint() == entryPoint
        }

        with(action2) {
            getName() == actionName
            getImageUri() == imageName
            getCommands() == commands
            getMounts() == mounts
            !getFlags()
            !getEntrypoint()
        }
    }

    def 'should construct a pipeline'() {
        given:
        def actions = [new Action(),new Action()]
        def resources = new Resources().setProjectId("testId")
        def helper = new GooglePipelinesHelper()

        when:
        def pipe = helper.createPipeline(actions,resources)

        then:
        pipe.getActions().size() == 2
        pipe.getResources().getProjectId() == "testId"

    }

    def 'should configure resources correctly'() {
        given:
        def type = "testType"
        def projectId = "testProject"
        def zone = ["testZone1","testZone2"]
        def region = ["testRegion1","testRegion2"]
        def diskName = "testDisk"
        def scopes = ["scope1","scope2"]
        def preEmptible = true
        def helper = new GooglePipelinesHelper()

        when:
        def resources1 = helper.configureResources(type,projectId,zone,null,diskName,100, scopes,preEmptible)
        def resources2 = helper.configureResources(type,projectId,null,region,diskName,200,scopes,preEmptible)
        def resources3 = helper.configureResources(type,projectId,zone,null,diskName)

        then:
        with(resources1) {
            getVirtualMachine().getMachineType() == type
            getProjectId() == projectId
            getZones() == zone
            getRegions() == null
            getVirtualMachine().getDisks().get(0).getName() == diskName
            getVirtualMachine().getDisks().get(0).getSizeGb() == 100
            getVirtualMachine().getServiceAccount().getScopes() == scopes
            getVirtualMachine().getPreemptible() == preEmptible
        }

        with(resources2) {
            getVirtualMachine().getMachineType() == type
            getProjectId() == projectId
            getZones() == null
            getRegions() == region
            getVirtualMachine().getDisks().get(0).getName() == diskName
            getVirtualMachine().getDisks().get(0).getSizeGb() == 200
            getVirtualMachine().getServiceAccount().getScopes() == scopes
            getVirtualMachine().getPreemptible() == preEmptible
        }

        with(resources3) {
            getVirtualMachine().getMachineType() == type
            getProjectId() == projectId
            getZones() == zone
            getVirtualMachine().getDisks().get(0).getName() == diskName
            !getVirtualMachine().getDisks().get(0).getSizeGb()
            !getVirtualMachine().getServiceAccount().getScopes()
            !getVirtualMachine().getPreemptible()
        }
    }

    def 'should set up a mount correctly'() {
        given:
        def diskName = "testDisk"
        def mountPath = "testPath"
        def readOnly = true
        def helper = new GooglePipelinesHelper()

        when:
        def mount1 = helper.configureMount(diskName,mountPath,readOnly)
        def mount2 = helper.configureMount(diskName,mountPath)

        then:
        with(mount1) {
            getDisk() == diskName
            getPath() == mountPath
            getReadOnly() == readOnly
        }

        with(mount2) {
            getDisk() == diskName
            getPath() == mountPath
            !getReadOnly()
        }
    }

    def 'should check on the status of an operation'() {
        given:
        def op = new Operation().setName("testOperation")
        def helper = new GooglePipelinesHelper()
        helper.genomicsClient = mockClient

        when:
        def opStatus = helper.checkOperationStatus(op)

        then:
        1 * mockClient.projects().operations().get(op.getName()).execute() >> new Operation().setName("status").setDone(true)
        opStatus.getName() == "status"
        opStatus.getDone()

    }

    def 'should cancel an operation'() {
        given:
        def op = new Operation().setName("testOperation").setDone(false)
        def helper = new GooglePipelinesHelper()
        helper.genomicsClient = mockClient

        when:
        helper.cancelOperation(op)

        then:
        1 * mockClient.projects().operations().cancel(op.getName(),new CancelOperationRequest()).execute() >> {}
    }

    def 'should run a pipeline'() {
        given:
        def pipe = new Pipeline()
        def pipeRequest = new RunPipelineRequest().setPipeline(pipe)

        def helper = new GooglePipelinesHelper()
        helper.genomicsClient = mockClient

        when:
        def runOp = helper.runPipeline(pipe)

        then:
        1 * mockClient.pipelines().run(pipeRequest).execute() >> new Operation().setName("runOp")

        runOp.getName() == "runOp"
    }

    def 'should create a resource request' () {
        given:
        def helper = Spy(GooglePipelinesHelper)
        def req = Mock(GooglePipelinesSubmitRequest)

        when:
        def res = helper.createResources(req)
        then:
        req.machineType >> 'n1-abc'
        req.project >> 'my-project'
        req.zone >> ['my-zone']
        req.region >> ['my-region']
        req.diskName >> 'my-disk'
        req.diskSizeGb >> 500
        req.preemptible >> true

        1 * helper.configureResources(
                'n1-abc',
                'my-project',
                ['my-zone'],
                ['my-region'],
                'my-disk',
                500,
                [GooglePipelinesHelper.SCOPE_CLOUD_PLATFORM],
                true )

        res.getProjectId() == 'my-project'
        res.getZones() ==['my-zone']
        res.getRegions() == ['my-region']
        res.getVirtualMachine().getPreemptible()
        res.getVirtualMachine().getMachineType() == 'n1-abc'
        res.getVirtualMachine().getDisks()[0].getName() == 'my-disk'
        res.getVirtualMachine().getDisks()[0].getSizeGb() == 500
    }

    def 'should create main action' () {
        given:
        def helper = Spy(GooglePipelinesHelper)
        def req = Mock(GooglePipelinesSubmitRequest)
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
                [mount],
                [IGNORE_EXIT_STATUS]
        )

        action.getName() == 'foo-main'
        action.getImageUri() == 'my/image'
        action.getCommands() == ['bash', '-c', 'something.sh']
        action.getMounts() == [mount]
        action.getFlags() == ['IGNORE_EXIT_STATUS']
        action.getEntrypoint() == null
        action.getEnvironment() == [:]
    }

    def 'should create staging action' () {
        given:
        def helper = Spy(GooglePipelinesHelper)
        def req = Mock(GooglePipelinesSubmitRequest)
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
                [mount],
                [ALWAYS_RUN, IGNORE_EXIT_STATUS]
        )

        action.getName() == 'bar-stage'
        action.getImageUri() == 'alpine'
        action.getCommands() == ['bash', '-c', 'copy this and that']
        action.getMounts() == [mount]
        action.getFlags() == ['ALWAYS_RUN', 'IGNORE_EXIT_STATUS']
        action.getEntrypoint() == null
        action.getEnvironment() == [:]
    }

    def 'should create unstaging action' () {
        given:
        def helper = Spy(GooglePipelinesHelper)
        def req = Mock(GooglePipelinesSubmitRequest)
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
                [ALWAYS_RUN, IGNORE_EXIT_STATUS]
        )

        action.getName() == 'bar-unstage'
        action.getImageUri() == 'alpine'
        action.getCommands() == ['bash', '-c', 'upload command here']
        action.getMounts() == [mount]
        action.getFlags() == ['ALWAYS_RUN', 'IGNORE_EXIT_STATUS']
        action.getEntrypoint() == null
        action.getEnvironment() == [:]
    }

    def 'should submit pipeline request' () {
        given:
        def helper = Spy(GooglePipelinesHelper)
        def req = Mock(GooglePipelinesSubmitRequest)
        def stage = GroovyMock(Action)
        def unstage = GroovyMock(Action)
        def main = GroovyMock(Action)
        def res = GroovyMock(Resources)
        def pipeline = GroovyMock(Pipeline)
        def operation = GroovyMock(Operation)

        when:
        def result = helper.submitPipeline(req)

        then:
        req.taskName >> 'foo'
        1 * helper.createStagingAction(req) >> stage
        1 * helper.createUnstagingAction(req) >> unstage
        1 * helper.createMainAction(req) >> main
        1 * helper.createResources(req) >> res
        1 * helper.createPipeline([stage, main, unstage], res) >> pipeline
        1 * helper.runPipeline(pipeline, [taskName: 'foo']) >> operation
        result == operation
    }
}
