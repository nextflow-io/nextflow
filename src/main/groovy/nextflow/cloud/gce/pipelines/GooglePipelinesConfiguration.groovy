package nextflow.cloud.gce.pipelines

import groovy.transform.CompileStatic

@CompileStatic
class GooglePipelinesConfiguration {
    String project
    String zone
    String vmInstanceType
    boolean preemptible

    GooglePipelinesConfiguration(String project, String zone, String vmInstanceType, boolean preemptible = false) {
        this.project = project
        this.zone = zone
        this.vmInstanceType = vmInstanceType
        this.preemptible = preemptible
    }


    @Override
    String toString() {
        "GooglePipelinesConfiguration{" +
                "project='" + project + '\'' +
                ", zone='" + zone + '\'' +
                ", vmInstanceType='" + vmInstanceType + '\'' +
                ", preemptible=" + preemptible +
                '}'
    }
}