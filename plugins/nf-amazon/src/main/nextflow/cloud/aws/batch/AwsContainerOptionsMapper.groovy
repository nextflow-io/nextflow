/*
 * Copyright 2013-2024, Seqera Labs
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
package nextflow.cloud.aws.batch

import com.amazonaws.services.batch.model.ContainerProperties
import com.amazonaws.services.batch.model.KeyValuePair
import com.amazonaws.services.batch.model.LinuxParameters
import com.amazonaws.services.batch.model.Tmpfs
import com.amazonaws.services.batch.model.Ulimit
import groovy.transform.CompileStatic
import nextflow.util.CmdLineOptionMap
import nextflow.util.MemoryUnit

/**
 * Maps task container options to AWS container properties
 *
 * @see https://docs.docker.com/engine/reference/commandline/run/
 * @see https://docs.aws.amazon.com/batch/latest/APIReference/API_ContainerProperties.html
 *
 * @author Manuele Simi <manuele.simi@gmail.com>
 */
@CompileStatic
class AwsContainerOptionsMapper {

    @Deprecated
    static ContainerProperties createContainerOpts(CmdLineOptionMap options) {
        createContainerProperties(options)
    }

    static ContainerProperties createContainerProperties(CmdLineOptionMap options) {
        final containerProperties = new ContainerProperties()
        if ( options?.hasOptions() ) {
            checkPrivileged(options, containerProperties)
            checkEnvVars(options, containerProperties)
            checkUser(options, containerProperties)
            checkReadOnly(options, containerProperties)
            checkUlimit(options, containerProperties)
            LinuxParameters params = checkLinuxParameters(options)
            if ( params != null )
                containerProperties.setLinuxParameters(params)
        }
        return containerProperties
    }

    protected static void checkPrivileged(CmdLineOptionMap options, ContainerProperties containerProperties) {
        if ( findOptionWithBooleanValue(options, 'privileged') )
            containerProperties.setPrivileged(true);
    }

    protected static void checkEnvVars(CmdLineOptionMap options, ContainerProperties containerProperties) {
        final keyValuePairs = new ArrayList<KeyValuePair>()
        List<String> values = findOptionWithMultipleValues(options, 'env')
        values.addAll(findOptionWithMultipleValues(options, 'e'))
        for( String it : values ) {
            final tokens = it.tokenize('=')
            keyValuePairs << new KeyValuePair().withName(tokens[0]).withValue(tokens.size() == 2 ? tokens[1] : null)
        }
        if ( keyValuePairs )
            containerProperties.setEnvironment(keyValuePairs)
    }

    protected static void checkUser(CmdLineOptionMap options, ContainerProperties containerProperties) {
        String user = findOptionWithSingleValue(options, 'u')
        if ( !user)
            user = findOptionWithSingleValue(options, 'user')
        if ( user )
            containerProperties.setUser(user)
    }

    protected static void checkReadOnly(CmdLineOptionMap options, ContainerProperties containerProperties) {
        if ( findOptionWithBooleanValue(options, 'read-only') )
            containerProperties.setReadonlyRootFilesystem(true);
    }

    protected static void checkUlimit(CmdLineOptionMap options, ContainerProperties containerProperties) {
        final ulimits = new ArrayList<Ulimit>()
        findOptionWithMultipleValues(options, 'ulimit').each { value ->
            final tokens = value.tokenize('=')
            final limits = tokens[1].tokenize(':')
            if ( limits.size() > 1 )
                ulimits << new Ulimit().withName(tokens[0])
                        .withSoftLimit(limits[0] as Integer).withHardLimit(limits[1] as Integer)
            else
                ulimits << new Ulimit().withName(tokens[0]).withSoftLimit(limits[0] as Integer)
        }
        if ( ulimits.size() )
            containerProperties.setUlimits(ulimits)
    }

    protected static LinuxParameters checkLinuxParameters(CmdLineOptionMap options) {
        final params = new LinuxParameters()
        boolean atLeastOneSet = false

        // shared Memory Size
        def value = findOptionWithSingleValue(options, 'shm-size')
        if ( value ) {
            final sharedMemorySize = MemoryUnit.of(value)
            params.setSharedMemorySize(sharedMemorySize.mega as Integer)
            atLeastOneSet = true
        }

        // tmpfs mounts, e.g --tmpfs /run:rw,noexec,nosuid,size=64
        final tmpfs = new ArrayList<Tmpfs>()
        findOptionWithMultipleValues(options, 'tmpfs').each { ovalue ->
            def matcher = ovalue =~ /^(?<path>.*):(?<options>.*?),size=(?<sizeMiB>.*)$/
            if (matcher.matches()) {
                tmpfs << new Tmpfs().withContainerPath(matcher.group('path'))
                        .withSize(matcher.group('sizeMiB') as Integer)
                        .withMountOptions(matcher.group('options').tokenize(','))
            } else {
                throw new IllegalArgumentException("Found a malformed value '${ovalue}' for --tmpfs option")
            }
        }
        if ( tmpfs ) {
            params.setTmpfs(tmpfs)
            atLeastOneSet = true
        }

        // swap limit equal to memory plus swap
        value = findOptionWithSingleValue(options, 'memory-swap')
        if ( value ) {
            params.setMaxSwap(value as Integer)
            atLeastOneSet = true
        }

        // run an init inside the container
        if ( findOptionWithBooleanValue(options, 'init') ) {
            params.setInitProcessEnabled(true)
            atLeastOneSet = true
        }

        // tune container memory swappiness
        value = findOptionWithSingleValue(options, 'memory-swappiness')
        if ( value ) {
            params.setSwappiness(value as Integer)
            atLeastOneSet = true
        }

        return atLeastOneSet ? params : null
    }

    /**
     * Finds the value of an option
     * @param name the name of the option
     * @return the value, if any, or empty
     */
    protected static String findOptionWithSingleValue(CmdLineOptionMap options, String name) {
        options.getFirstValueOrDefault(name,null) as String
    }

    /**
     * Finds the values of an option that can be repeated
     * @param name the name of the option
     * @return the list of values
     */
    protected static List<String> findOptionWithMultipleValues(CmdLineOptionMap options, String name) {
        options.getValues(name)
    }

    /**
     * Checks if a boolean flag exists
     * @param name the name of the flag
     * @return true if it exists, false otherwise
     */
    protected static boolean findOptionWithBooleanValue(CmdLineOptionMap options, String name) {
        options.exists(name) ? options.getFirstValue(name) as Boolean : false
    }
}
