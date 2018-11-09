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
                '}'
    }
}