package nextflow.cloud.gce.pipelines

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import nextflow.Session
import nextflow.exception.AbortOperationException
import spock.lang.Shared
import spock.lang.Specification

import java.nio.file.Path

class GooglePipelinesExecutorTest extends Specification {

    @Shared
    def validConfig = [
            "gce" : [
            "project" : "testProject",
            "zone" : "testZone"
            ]
    ]

    def 'should abort operation when the wordir is not a CloudStoragePath'() {
        given:
        def session = Stub(Session)
        session.workDir = Stub(Path)
        def executor = new GooglePipelinesExecutor()
        executor.session = session

        when:
        executor.register()

        then:
        def error = thrown(AbortOperationException)
        error.getMessage().contains("GCE bucket must be provided as a working directory")
    }

    def 'should abort operation when required configuration keys are missing'() {
        given:
        def session = Stub(Session)
        def path = CloudStorageFileSystem.forBucket("test").getPath("/")
        session.workDir >> path
        session.config >> [
                "gce" : [
                        (key) : configValue
                ]
        ]


        def executor = new GooglePipelinesExecutor()
        executor.session = session

        when:
        executor.register()

        then:
        def error = thrown(AbortOperationException)
        !error.getMessage().contains(configKey.toString())

        where:
        key         |   configKey     |   configValue
        "project"   | "gce.project"   |   "testProject"
        "zone"      | "gce.zone"      |   "testZone"
    }

    def 'should register successfully'()  {
        given:
        def session = Stub(Session)
        def helper = Mock(GooglePipelinesHelper)
        def path = CloudStorageFileSystem.forBucket("test").getPath("/")
        session.workDir >> path
        session.config >> validConfig
        def executor = new GooglePipelinesExecutor(helper)
        executor.session = session

        when:
        executor.register()

        then:
        executor.pipelineConfig.project == validConfig.gce?.project
    }

    def 'should be containerNative'() {
        when:
        def executor = new GooglePipelinesExecutor()

        then:
        executor.isContainerNative()
    }
}
