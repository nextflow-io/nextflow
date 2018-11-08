package nextflow.cloud.gce.pipelines

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import nextflow.Session
import nextflow.exception.AbortOperationException
import spock.lang.Ignore
import spock.lang.Shared
import spock.lang.Specification

import java.nio.file.Path

class GooglePipelinesExecutorTest extends Specification {

    @Shared
    def validZoneConfig = [
            "gce" : [
            "project" : "testProject",
            "zone" : "testZone1,testZone2"
            ]
    ]

    @Shared
    def validRegionConfig = [
            "gce" : [
                    "project" : "testProject",
                    "region" : "testRegion1,testRegion2"
            ]
    ]

    def 'should abort operation when the workdir is not a CloudStoragePath'() {
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

    def 'should abort operation when project is not specified'() {
        given:
        def session = Stub(Session)
        def path = CloudStorageFileSystem.forBucket("test").getPath("/")
        session.workDir >> path
        session.config >> [
                "gce" : [
                        "zone" : "testZone"
                ]
        ]
        def executor = new GooglePipelinesExecutor()
        executor.session = session

        when:
        executor.register()

        then:
        def error = thrown(AbortOperationException)
        error.getMessage() == "Required config value 'gce.project' for executor null is not defined. Please add it to your process or nextflow configuration file."
    }

    def 'should abort operation when neither zone or region are specified'() {
        given:
        def session = Stub(Session)
        def path = CloudStorageFileSystem.forBucket("test").getPath("/")
        session.workDir >> path
        session.config >> [
                "gce" : [
                        "project" : "testproject"
                ]
        ]
        def executor = new GooglePipelinesExecutor()
        executor.session = session

        when:
        executor.register()

        then:
        def error = thrown(AbortOperationException)
        error.getMessage().contains("Missing configuration value 'gce.zone' or 'gce.region'")
    }


    def 'should abort operation when both zone and region are specified'() {
        given:
        def session = Stub(Session)
        def path = CloudStorageFileSystem.forBucket("test").getPath("/")
        session.workDir >> path
        session.config >> [
                "gce" : [
                        "project" : "testproject",
                        "zone" : "testZone",
                        "region" : "testRegion"
                ]
        ]
        def executor = new GooglePipelinesExecutor()
        executor.session = session

        when:
        executor.register()

        then:
        def error = thrown(AbortOperationException)
        error.getMessage().contains("You can't specify both 'gce.zone' and 'gce.region' configuration parameters. Please remove one of them from your configuration.")
    }



    @Ignore
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

    def 'should register successfully with zone'()  {
        given:
        def session = Stub(Session)
        def helper = Mock(GooglePipelinesHelper)
        def path = CloudStorageFileSystem.forBucket("test").getPath("/")
        session.workDir >> path
        session.config >> validZoneConfig
        def executor = new GooglePipelinesExecutor(helper)
        executor.session = session

        when:
        executor.register()

        then:
        executor.pipelineConfig.project == validZoneConfig.gce?.project
        executor.pipelineConfig.zone == validZoneConfig.gce?.zone?.split(",")?.toList()
    }

    def 'should register successfully with region'()  {
        given:
        def session = Stub(Session)
        def helper = Mock(GooglePipelinesHelper)
        def path = CloudStorageFileSystem.forBucket("test").getPath("/")
        session.workDir >> path
        session.config >> validRegionConfig
        def executor = new GooglePipelinesExecutor(helper)
        executor.session = session

        when:
        executor.register()

        then:
        executor.pipelineConfig.project == validRegionConfig.gce?.project
        executor.pipelineConfig.region == validRegionConfig.gce?.region?.split(",")?.toList()
    }

    def 'should be containerNative'() {
        when:
        def executor = new GooglePipelinesExecutor()

        then:
        executor.isContainerNative()
    }
}
