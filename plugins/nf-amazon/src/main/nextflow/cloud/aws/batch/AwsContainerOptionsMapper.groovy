/*
 * Copyright 2013-2026, Seqera Labs
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

import nextflow.cloud.aws.batch.model.ContainerPropertiesModel
import software.amazon.awssdk.services.batch.model.KeyValuePair
import software.amazon.awssdk.services.batch.model.LinuxParameters
import software.amazon.awssdk.services.batch.model.Tmpfs
import software.amazon.awssdk.services.batch.model.Ulimit
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
    static ContainerPropertiesModel createContainerOpts(CmdLineOptionMap options) {
        return createContainerProperties(options)
    }

    static ContainerPropertiesModel createContainerProperties(CmdLineOptionMap options) {
        final containerProperties = new ContainerPropertiesModel()
        if ( options?.hasOptions() ) {
            checkPrivileged(options, containerProperties)
            checkEnvVars(options, containerProperties)
            checkUser(options, containerProperties)
            checkReadOnly(options, containerProperties)
            checkUlimit(options, containerProperties)
            LinuxParameters params = checkLinuxParameters(options)
            if ( params != null )
                containerProperties.linuxParameters(params)
        }
        return containerProperties
    }

    protected static void checkPrivileged(CmdLineOptionMap options, ContainerPropertiesModel containerProperties) {
        if ( findOptionWithBooleanValue(options, 'privileged') )
            containerProperties.privileged(true)
    }

    protected static void checkEnvVars(CmdLineOptionMap options, ContainerPropertiesModel containerProperties) {
        final keyValuePairs = new ArrayList<KeyValuePair>()
        List<String> values = findOptionWithMultipleValues(options, 'env')
        values.addAll(findOptionWithMultipleValues(options, 'e'))
        for( String it : values ) {
            final tokens = it.tokenize('=')
            keyValuePairs << KeyValuePair.builder().name(tokens[0]).value(tokens.size() == 2 ? tokens[1] : null).build()
        }
        if ( keyValuePairs )
            containerProperties.environment(keyValuePairs)
    }

    protected static void checkUser(CmdLineOptionMap options, ContainerPropertiesModel containerProperties) {
        String user = findOptionWithSingleValue(options, 'u')
        if ( !user)
            user = findOptionWithSingleValue(options, 'user')
        if ( user )
            containerProperties.user(user)
    }

    protected static void checkReadOnly(CmdLineOptionMap options, ContainerPropertiesModel containerProperties) {
        if ( findOptionWithBooleanValue(options, 'read-only') )
            containerProperties.readonlyRootFilesystem(true);
    }

    protected static void checkUlimit(CmdLineOptionMap options, ContainerPropertiesModel containerProperties) {
        final ulimits = new ArrayList<Ulimit>()
        findOptionWithMultipleValues(options, 'ulimit').each { value ->
            final tokens = value.tokenize('=')
            final limits = tokens[1].tokenize(':')
            if ( limits.size() > 1 )
                ulimits << Ulimit.builder().name(tokens[0]).softLimit(limits[0] as Integer).hardLimit(limits[1] as Integer).build()
            else
                ulimits << Ulimit.builder().name(tokens[0]).softLimit(limits[0] as Integer).build()
        }
        if ( ulimits.size() )
            containerProperties.ulimits(ulimits)
    }

    protected static LinuxParameters checkLinuxParameters(CmdLineOptionMap options) {
        final params = LinuxParameters.builder()
        boolean atLeastOneSet = false

        // shared Memory Size
        def value = findOptionWithSingleValue(options, 'shm-size')
        if ( value ) {
            final sharedMemorySize = MemoryUnit.of(value)
            params.sharedMemorySize(sharedMemorySize.mega as Integer)
            atLeastOneSet = true
        }

        // tmpfs mounts, e.g --tmpfs /run:rw,noexec,nosuid,size=64
        final tmpfs = new ArrayList<Tmpfs>()
        findOptionWithMultipleValues(options, 'tmpfs').each { ovalue ->
            def matcher = ovalue =~ /^(?<path>.*):(?<options>.*?),size=(?<sizeMiB>.*)$/
            if (matcher.matches()) {
                tmpfs << Tmpfs.builder().containerPath(matcher.group('path'))
                        .size(matcher.group('sizeMiB') as Integer)
                        .mountOptions(matcher.group('options').tokenize(','))
                        .build()
            } else {
                throw new IllegalArgumentException("Found a malformed value '${ovalue}' for --tmpfs option")
            }
        }
        if ( tmpfs ) {
            params.tmpfs(tmpfs)
            atLeastOneSet = true
        }

        // swap limit equal to memory plus swap
        value = findOptionWithSingleValue(options, 'memory-swap')
        if ( value ) {
            params.maxSwap(value as Integer)
            atLeastOneSet = true
        }

        // run an init inside the container
        if ( findOptionWithBooleanValue(options, 'init') ) {
            params.initProcessEnabled(true)
            atLeastOneSet = true
        }

        // tune container memory swappiness
        value = findOptionWithSingleValue(options, 'memory-swappiness')
        if ( value ) {
            params.swappiness(value as Integer)
            atLeastOneSet = true
        }

        return atLeastOneSet ? params.build() : null
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
