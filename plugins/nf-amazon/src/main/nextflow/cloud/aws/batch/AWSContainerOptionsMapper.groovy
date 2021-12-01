package nextflow.cloud.aws.batch

import com.amazonaws.services.batch.model.ContainerProperties
import com.amazonaws.services.batch.model.KeyValuePair
import com.amazonaws.services.batch.model.LinuxParameters
import com.amazonaws.services.batch.model.Tmpfs
import com.amazonaws.services.batch.model.Ulimit

/**
 * Maps task container options to AWS container properties
 *
 * @see https://docs.docker.com/engine/reference/commandline/run/
 * @see https://docs.aws.amazon.com/batch/latest/APIReference/API_ContainerProperties.html
 *
 * @author Manuele Simi <manuele.simi@gmail.com>
 */
class AWSContainerOptionsMapper {

    def options

    protected AWSContainerOptionsMapper(String containerOptions) {
        options = containerOptions.trim().split("\\s+")
    }

    protected ContainerProperties addProperties(ContainerProperties containerProperties) {
        if (options.size() > 0) {
            checkPrivileged(containerProperties)
            checkEnvVars(containerProperties)
            checkUser(containerProperties)
            checkReadOnly(containerProperties)
            checkUlimit(containerProperties)
            checkLinuxParameters(containerProperties)
        }
        return containerProperties
    }

    protected void checkPrivileged(ContainerProperties containerProperties) {
        if ( findOptionWithBooleanValue('--privileged') )
            containerProperties.setPrivileged(true);
    }

    protected void checkEnvVars(ContainerProperties containerProperties) {
        final keyValuePairs = new ArrayList<KeyValuePair>()
        def values = findOptionWithMultipleValues('--env')
        values.addAll(findOptionWithMultipleValues('-e'))
        values.each { value ->
            final tokens = value.split('=')
            keyValuePairs << new KeyValuePair().withName(tokens[0]).withValue(tokens.size() == 2 ? tokens[1] : null)
        }
        if ( keyValuePairs.size() > 0 )
            containerProperties.setEnvironment(keyValuePairs)
    }

    protected void checkUser(ContainerProperties containerProperties) {
        def user = findOptionWithSingleValue('-u')
        if ( !user )
            user = findOptionWithSingleValue('--user')
        if ( user )
            containerProperties.setUser(user)
    }

    protected void checkReadOnly(ContainerProperties containerProperties) {
        if ( findOptionWithBooleanValue('--read-only') )
            containerProperties.setReadonlyRootFilesystem(true);
    }

    protected void checkUlimit(ContainerProperties containerProperties) {
        final ulimits = new ArrayList<Ulimit>()
        findOptionWithMultipleValues('--ulimit').each { value ->
            final tokens = value.tokenize('=')
            final limits = tokens[1].tokenize(':')
            if (limits.size() > 1)
                ulimits << new Ulimit().withName(tokens[0])
                        .withSoftLimit(limits[0] as Integer).withHardLimit(limits[1] as Integer)
            else
                ulimits << new Ulimit().withName(tokens[0]).withSoftLimit(limits[0] as Integer)
        }
        if (ulimits.size() > 0)
            containerProperties.setUlimits(ulimits)
    }

    protected void checkLinuxParameters(ContainerProperties containerProperties) {
        final params = new LinuxParameters()

        // shared Memory Size
        def value = findOptionWithSingleValue('--shm-size')
        if (value)
            params.setSharedMemorySize(value as Integer)

        // tmpfs mounts
        final tmpfs = new ArrayList<Tmpfs>()
        findOptionWithMultipleValues('--tmpfs').each { path ->
            tmpfs << new Tmpfs().withContainerPath(path)
        }
        if ( tmpfs.size() > 0 )
            params.setTmpfs(tmpfs)

        // swap limit equal to memory plus swap
        value = findOptionWithSingleValue('--memory-swap')
        if ( value )
            params.setMaxSwap(value as Integer)

        // run an init inside the container
        value = findOptionWithBooleanValue('--init')
        if ( value )
            params.setInitProcessEnabled(value)

        // tune container memory swappiness
        value = findOptionWithSingleValue('--memory-swappiness')
        if ( value )
            params.setSwappiness(value as Integer)

        containerProperties.setLinuxParameters(params)
    }

    /**
     * Finds the value of an option
     * @param name the name of the option
     * @return the value, if any, or empty
     */
    protected def findOptionWithSingleValue(def name) {
        def index = options.findIndexOf({ it == name })
        if (index != -1) {
            if ( !isValidValue(options[index + 1] as String) )
                throw new IllegalArgumentException("Found a malformed option '${name}' for the job container")
            return options[index + 1] as String
        }
        return ''
    }

    /**
     * Finds the values of an option that can be repeated
     * @param name the name of the option
     * @return the list of values
     */
    protected def findOptionWithMultipleValues(String name) {
        final values = new ArrayList<String>()
        options.findIndexValues{ it == name}.collect { it as Integer }
                .each { index ->
                    if ( !isValidValue(options[index + 1]) )
                        throw new IllegalArgumentException("Found a malformed option ${name} for the job container")
                    values << options[index + 1]
                }
        return  values
    }

    /**
     * Checks if a boolean flag exists
     * @param name the name of the flag
     * @return true if it exists, false otherwise
     */
    protected def findOptionWithBooleanValue(def name) {
        options.find { it == name }
    }

    /**
     * Checks if the value of an option is valid
     * @param value the value to check
     * @return true if the value is valid, false otherwise
     */
    protected boolean isValidValue(String value) {
        !value.startsWith('--') && !value.startsWith('-')
    }
}