/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.cloud

import groovy.transform.stc.ClosureParams
import groovy.transform.stc.SimpleType
import nextflow.cloud.types.CloudInstanceStatus
import nextflow.cloud.types.CloudInstanceType
import nextflow.util.ServiceName

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ServiceName('fake')
class FakeCloudDriver implements CloudDriver {


    @Override
    List<String> launchInstances(int instanceCount, LaunchConfig config) {
        return null
    }

    @Override
    void waitInstanceStatus(Collection<String> instanceIds, CloudInstanceStatus status) {

    }

    @Override
    void tagInstances(Collection<String> instanceIds, Map<String, String> tags) {

    }

    @Override
    void eachSpotPrice(List<String> instanceTypes, @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudSpotPrice']) Closure callback) {

    }

    @Override
    void eachInstanceWithTags(Map tags, Closure callback) {

    }

    @Override
    void eachInstanceWithIds(List<String> instanceIds, Closure callback) {

    }

    @Override
    void eachInstance(Closure callback) {

    }

    @Override
    List<String> listPrivateIPs(String clusterName) {
        return null
    }

    @Override
    void terminateInstances(Collection<String> instanceIds) {

    }

    @Override
    String getLocalInstanceId() {
        return null
    }

    @Override
    String getLocalTerminationNotice() {
        return null
    }

    @Override
    CloudInstanceType describeInstanceType(String instanceType) {
        return null
    }

    void validate(LaunchConfig config) {

    }
}
