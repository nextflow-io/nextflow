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
import groovy.transform.CompileStatic
import nextflow.config.CascadingConfig
import nextflow.config.ConfigField
import nextflow.util.MemoryUnit
/**
 * Model a cloud instance launch configuration properties
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class LaunchConfig extends CascadingConfig<String,Object> {

    LaunchConfig() {}

    LaunchConfig(Map config, CascadingConfig<String,Object> fallback) {
        super(config, fallback)
    }

    @ConfigField(_private = true)
    String getRole() {
        getAttribute('role')
    }

    @ConfigField(_private = true)
    String getClusterName() {
        getAttribute('clusterName')
    }

    @ConfigField(_private = true)
    boolean getCreateUser() {
        getAttribute('createUser')
    }

    boolean isSpot() { spotPrice != null }

    @ConfigField
    String getUserName() {
        getAttribute('userName')
    }

    @ConfigField
    String getImageId() {
        getAttribute('imageId')
    }

    @ConfigField
    String getInstanceType() {
        getAttribute('instanceType')
    }

    @ConfigField
    String getKeyName() {
        getAttribute('keyName')
    }

    @ConfigField
    String getKeyHash() {
        getAttribute('keyHash')
    }

    @ConfigField
    String getSubnetId() {
        getAttribute('subnetId')
    }

    @ConfigField
    String getInstanceRole() {
        getAttribute('instanceRole')
    }

    LaunchConfig setInstanceRole(String profile) {
        setAttribute('instanceRole',profile)
        return this
    }

    @ConfigField
    List<String> getSecurityGroup() {
        getSecurityGroups()
    }

    @ConfigField
    List<String> getSecurityGroups() {
        String groups = getAttribute('securityGroup')
        if( !groups )
            groups = getAttribute('securityGroups')

        return groups ? groups.tokenize(',') : Collections.<String>emptyList()
    }

    @ConfigField
    MemoryUnit getBootStorageSize() {
        getAttribute('bootStorageSize') as MemoryUnit
    }


    @ConfigField
    String getInstanceStorageMount() {
        getAttribute('instanceStorageMount')
    }

    @ConfigField
    String getInstanceStorageDevice() {
        getAttribute('instanceStorageDevice')
    }

    @ConfigField
    String getSharedStorageMount() {
        getAttribute('sharedStorageMount') ?: '/mnt/efs'
    }

    @ConfigField
    String getSharedStorageId() {
        getAttribute('sharedStorageId')
    }

    @ConfigField
    String getSpotPrice() {
        getAttribute('spotPrice')
    }

    @ConfigField
    CloudConfig.Nextflow getNextflow() {
        (CloudConfig.Nextflow)getAttribute('nextflow')
    }

    @ConfigField
    List<String> getDockerPull() {
        def val = getAttribute('dockerPull')
        if( !val )
            return Collections.<String>emptyList()

        if( val instanceof List )
            return (List)val

        if( val instanceof CharSequence ) {
            return [val.toString()]
        }

        throw new IllegalArgumentException("Valid a valid `dockerPull` argument: $val [${val.class.name}]")
    }

    @Override
    ConfigObject toConfigObject() {
        def copy = new LinkedHashMap(this.config)
        copy.remove('role')
        toConfigObject0(copy)
    }

    ConfigObject toCloudConfigObject() {
        def current = this
        while( current && current.getClass() != CloudConfig ) {
            current = this.parent
        }

        def obj = current ? current.toConfigObject() : new ConfigObject()
        def result = new ConfigObject()
        result.put('cloud', obj)
        return result
    }

    String renderCloudConfigObject() {
        toCloudConfigObject().prettyPrint()
    }


}
