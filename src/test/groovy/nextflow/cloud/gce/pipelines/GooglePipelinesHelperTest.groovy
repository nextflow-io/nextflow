package nextflow.cloud.gce.pipelines

import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.google.api.services.genomics.v2alpha1.Genomics
import com.google.api.services.genomics.v2alpha1.model.Action
import com.google.api.services.genomics.v2alpha1.model.CancelOperationRequest
import com.google.api.services.genomics.v2alpha1.model.Operation
import com.google.api.services.genomics.v2alpha1.model.Pipeline
import com.google.api.services.genomics.v2alpha1.model.Resources
import com.google.api.services.genomics.v2alpha1.model.RunPipelineRequest
import spock.lang.Shared
import spock.lang.Specification

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
        def zone = "testZone"
        def diskName = "testDisk"
        def scopes = ["scope1","scope2"]
        def preEmptible = true
        def helper = new GooglePipelinesHelper()

        when:
        def resources1 = helper.configureResources(type,projectId,zone,diskName,scopes,preEmptible)
        def resources2 = helper.configureResources(type,projectId,zone,diskName)

        then:
        with(resources1) {
            getVirtualMachine().getMachineType() == type
            getProjectId() == projectId
            getZones() == [zone]
            getVirtualMachine().getDisks().get(0).getName() == diskName
            getVirtualMachine().getServiceAccount().getScopes() == scopes
            getVirtualMachine().getPreemptible() == preEmptible
        }

        with(resources2) {
            getVirtualMachine().getMachineType() == type
            getProjectId() == projectId
            getZones() == [zone]
            getVirtualMachine().getDisks().get(0).getName() == diskName
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
}
