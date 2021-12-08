package nextflow.cloud.aws.batch

import com.amazonaws.services.batch.model.ContainerProperties
import com.amazonaws.services.batch.model.KeyValuePair
import com.amazonaws.services.batch.model.LinuxParameters
import com.amazonaws.services.batch.model.Tmpfs
import com.amazonaws.services.batch.model.Ulimit
import groovy.transform.CompileStatic
import nextflow.util.CmdLineOptionMap

/**
 * Maps task container options to AWS container properties
 *
 * @see <a href="https://docs.docker.com/engine/reference/commandline/run/">Docker run</a>
 * @see <a href="https://docs.aws.amazon.com/batch/latest/APIReference/API_ContainerProperties.html">API Container Properties</a>
 * @author Manuele Simi <manuele.simi@gmail.com>
 */
@CompileStatic
class AwsContainerOptionsMapper {

    final CmdLineOptionMap options

    protected AwsContainerOptionsMapper(CmdLineOptionMap containerOptions) {
        options = containerOptions
    }

    protected ContainerProperties addProperties(ContainerProperties containerProperties) {
        if ( options.hasOptions() ) {
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
        if ( findOptionWithBooleanValue('privileged') )
            containerProperties.setPrivileged(true);
    }

    protected void checkEnvVars(ContainerProperties containerProperties) {
        final keyValuePairs = new ArrayList<KeyValuePair>()
        List<String> values = findOptionWithMultipleValues('env')
        values.addAll(findOptionWithMultipleValues('e'))
        values.each { String value ->
            final tokens = value.tokenize('=')
            keyValuePairs << new KeyValuePair().withName(tokens[0]).withValue(tokens.size() == 2 ? tokens[1] : null)
        }
        if ( keyValuePairs.size() > 0 )
            containerProperties.setEnvironment(keyValuePairs)
    }

    protected void checkUser(ContainerProperties containerProperties) {
        String user = findOptionWithSingleValue('u')
        if ( !user || user.length() == 0 )
            user = findOptionWithSingleValue('user')
        if ( user && user.length() > 0 )
            containerProperties.setUser(user)
    }

    protected void checkReadOnly(ContainerProperties containerProperties) {
        if ( findOptionWithBooleanValue('read-only') )
            containerProperties.setReadonlyRootFilesystem(true);
    }

    protected void checkUlimit(ContainerProperties containerProperties) {
        final ulimits = new ArrayList<Ulimit>()
        findOptionWithMultipleValues('ulimit').each { value ->
            final tokens = value.tokenize('=')
            final limits = tokens[1].tokenize(':')
            if ( limits.size() > 1 )
                ulimits << new Ulimit().withName(tokens[0])
                        .withSoftLimit(limits[0] as Integer).withHardLimit(limits[1] as Integer)
            else
                ulimits << new Ulimit().withName(tokens[0]).withSoftLimit(limits[0] as Integer)
        }
        if ( ulimits.size() > 0 )
            containerProperties.setUlimits(ulimits)
    }

    protected void checkLinuxParameters(ContainerProperties containerProperties) {
        final params = new LinuxParameters()
        boolean atLeastOneSet = false

        // shared Memory Size
        def value = findOptionWithSingleValue('shm-size')
        if ( value && value.length() > 0 ) {
            params.setSharedMemorySize(value as Integer)
            atLeastOneSet = true
        }

        // tmpfs mounts, e.g --tmpfs /run:rw,noexec,nosuid,size=64
        final tmpfs = new ArrayList<Tmpfs>()
        findOptionWithMultipleValues('tmpfs').each { ovalue ->
            def matcher = ovalue =~ /^(?<path>.*):(?<options>.*?),size=(?<sizeMiB>.*)$/
            if (matcher.matches()) {
                tmpfs << new Tmpfs().withContainerPath(matcher.group('path'))
                        .withSize(matcher.group('sizeMiB') as Integer)
                        .withMountOptions(matcher.group('options').tokenize(','))
            } else {
                throw new IllegalArgumentException("Found a malformed value '${ovalue}' for --tmpfs option")
            }
        }
        if ( tmpfs.size() > 0 ) {
            params.setTmpfs(tmpfs)
            atLeastOneSet = true
        }

        // swap limit equal to memory plus swap
        value = findOptionWithSingleValue('memory-swap')
        if ( value && value.length() > 0) {
            params.setMaxSwap(value as Integer)
            atLeastOneSet = true
        }

        // run an init inside the container
        if ( findOptionWithBooleanValue('init') ) {
            params.setInitProcessEnabled(true)
            atLeastOneSet = true
        }

        // tune container memory swappiness
        value = findOptionWithSingleValue('memory-swappiness')
        if ( value && value.length() > 0 ) {
            params.setSwappiness(value as Integer)
            atLeastOneSet = true
        }

        if ( atLeastOneSet )
            containerProperties.setLinuxParameters(params)
    }

    /**
     * Finds the value of an option
     * @param name the name of the option
     * @return the value, if any, or empty
     */
    protected String findOptionWithSingleValue(String name) {
        options.getFirstValueOrDefault(name,'') as String
    }

    /**
     * Finds the values of an option that can be repeated
     * @param name the name of the option
     * @return the list of values
     */
    protected List<String> findOptionWithMultipleValues(String name) {
        options.getValues(name)
    }

    /**
     * Checks if a boolean flag exists
     * @param name the name of the flag
     * @return true if it exists, false otherwise
     */
    protected boolean findOptionWithBooleanValue(String name) {
        options.exists(name) ? options.getFirstValue(name) as Boolean : false
    }
}
