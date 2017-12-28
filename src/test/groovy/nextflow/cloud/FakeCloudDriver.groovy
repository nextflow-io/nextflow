/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
