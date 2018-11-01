package nextflow.cloud.gce.pipelines

import groovy.transform.CompileStatic

import java.nio.file.Path

@CompileStatic
class GooglePipelinesConfiguration {
    String project
    String zone
    String vmInstanceType
    boolean preemptible
    Path remoteBinDir

    GooglePipelinesConfiguration(String project, String zone, String vmInstanceType, Path remoteBinDir = null, boolean preemptible = false) {
        this.project = project
        this.zone = zone
        this.vmInstanceType = vmInstanceType
        this.remoteBinDir = remoteBinDir
        this.preemptible = preemptible
    }
}