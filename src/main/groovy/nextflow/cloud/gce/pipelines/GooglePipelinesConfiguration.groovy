package nextflow.cloud.gce.pipelines

import groovy.transform.CompileStatic

import java.nio.file.Path

@CompileStatic
class GooglePipelinesConfiguration {
    String project
    List<String> zone
    List<String> region
    String vmInstanceType
    boolean preemptible
    Path remoteBinDir

    GooglePipelinesConfiguration(String project, List<String> zone,List<String> region, String vmInstanceType, Path remoteBinDir = null, boolean preemptible = false) {
        this.project = project
        this.zone = zone
        this.region = region
        this.vmInstanceType = vmInstanceType
        this.remoteBinDir = remoteBinDir
        this.preemptible = preemptible
    }


    @Override
    String toString() {
        return "GooglePipelinesConfiguration{" +
                "project='" + project + '\'' +
                ", zone=" + zone +
                ", region=" + region +
                ", vmInstanceType='" + vmInstanceType + '\'' +
                ", preemptible=" + preemptible +
                ", remoteBinDir=" + remoteBinDir +
                '}';
    }
}