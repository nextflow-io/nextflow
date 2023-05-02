/*
 * Copyright 2013-2023, Seqera Labs
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

import static nextflow.cloud.aws.batch.AwsContainerOptionsMapper.*

import java.nio.file.Paths
import java.time.Instant

import com.amazonaws.services.batch.AWSBatch
import com.amazonaws.services.batch.model.AWSBatchException
import com.amazonaws.services.batch.model.ClientException
import com.amazonaws.services.batch.model.ContainerOverrides
import com.amazonaws.services.batch.model.ContainerProperties
import com.amazonaws.services.batch.model.DescribeJobDefinitionsRequest
import com.amazonaws.services.batch.model.DescribeJobDefinitionsResult
import com.amazonaws.services.batch.model.EvaluateOnExit
import com.amazonaws.services.batch.model.Host
import com.amazonaws.services.batch.model.JobDefinition
import com.amazonaws.services.batch.model.JobDefinitionType
import com.amazonaws.services.batch.model.JobTimeout
import com.amazonaws.services.batch.model.KeyValuePair
import com.amazonaws.services.batch.model.LogConfiguration
import com.amazonaws.services.batch.model.MountPoint
import com.amazonaws.services.batch.model.RegisterJobDefinitionRequest
import com.amazonaws.services.batch.model.RegisterJobDefinitionResult
import com.amazonaws.services.batch.model.ResourceRequirement
import com.amazonaws.services.batch.model.ResourceType
import com.amazonaws.services.batch.model.RetryStrategy
import com.amazonaws.services.batch.model.SubmitJobRequest
import com.amazonaws.services.batch.model.SubmitJobResult
import com.amazonaws.services.batch.model.Volume
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import nextflow.Const
import nextflow.container.ContainerNameValidator
import nextflow.exception.ProcessSubmitException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.fusion.FusionAwareTask
import nextflow.processor.TaskRun
import nextflow.util.CacheHelper
import org.slf4j.Logger
import org.slf4j.LoggerFactory
/**
 * Implementation of submit job requests for tasks.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
trait SubmitJobAware extends FusionAwareTask {

    static private Logger log = LoggerFactory.getLogger(SubmitJobAware)

    static private Map<String,String> jobDefinitions = [:]

    Map<String,String> environment = System.getenv()

    AWSBatch bypassProxy(AWSBatch client) {
        client instanceof AwsBatchProxy ? client.getClient() : client
    }

    /**
     * Get the AWS Batch client.
     */
    abstract AWSBatch getClient()

    /**
     * Get the Nextflow configuration options for AWS Batch.
     */
    abstract AwsOptions getAwsOptions()

    /**
     * Get the submit command for a job request.
     */
    abstract List<String> getSubmitCommand()

    /**
     * Create a new submit job request for a task.
     *
     * @param task
     */
    SubmitJobRequest newSubmitRequest(TaskRun task) {

        // create the request object
        final request = new SubmitJobRequest()
        request.setJobName(normalizeJobName(task.name))
        request.setJobQueue(getJobQueue(task))
        request.setJobDefinition(getJobDefinition(task))

        // set task timeout
        final time = task.config.getTime()
        if( time ) {
            def secs = time.toSeconds() as Integer
            if( secs < 60 ) {
                secs = 60   // Batch minimal allowed timeout is 60 seconds
            }
            request.setTimeout(new JobTimeout().withAttemptDurationSeconds(secs))
        }

        // set the container overrides
        final container = new ContainerOverrides()

        // set the submit command
        container.setCommand(getSubmitCommand())

        // set the environment vars
        final env = getEnvironmentVars()
        if( env )
            container.setEnvironment(env)

        // set the task resources
        final resources = new ArrayList<ResourceRequirement>(5)

        if( task.config.getMemory() ) {
            final mega = (int)task.config.getMemory().toMega()
            if( mega >= 4 )
                resources << new ResourceRequirement().withType(ResourceType.MEMORY).withValue(mega.toString())
            else
                log.warn "Ignoring task ${task.lazyName()} memory directive: ${task.config.getMemory()} -- AWS Batch job memory request cannot be lower than 4 MB"
        }

        if( task.config.getCpus() > 1 )
            resources << new ResourceRequirement().withType(ResourceType.VCPU).withValue(task.config.getCpus().toString())

        final accelerator = task.config.getAccelerator()
        if( accelerator ) {
            final type = accelerator.type?.toUpperCase() ?: 'GPU'
            final count = accelerator.request?.toString() ?: '1'

            resources << new ResourceRequirement().withType(type).withValue(count)
        }

        if( resources )
            container.withResourceRequirements(resources)

        request.setContainerOverrides(container)

        // set the resource labels
        final labels = task.config.getResourceLabels()
        if( labels ) {
            request.setTags(labels)
            request.setPropagateTags(true)
        }

        // set the share identifier
        if( awsOptions.shareIdentifier ) {
            request.setShareIdentifier(awsOptions.shareIdentifier)
            request.setSchedulingPriorityOverride(awsOptions.schedulingPriority)
        }

        // retry on spot reclamation
        // https://aws.amazon.com/blogs/compute/introducing-retry-strategies-for-aws-batch/
        if( awsOptions.maxSpotAttempts > 0 ) {
            final retryStrategy = new RetryStrategy()
                .withAttempts(awsOptions.maxSpotAttempts)
                .withEvaluateOnExit(
                    // retry the job when an EC2 instance is terminated
                    new EvaluateOnExit().withAction('RETRY').withOnStatusReason('Host EC2*'),
                    // delegate all other exit conditions to Nextflow
                    new EvaluateOnExit().withAction('EXIT').withOnReason('*')
                )

            request.setRetryStrategy(retryStrategy)
        }

        return request
    }

    /**
     * Remove invalid characters from a job name.
     *
     * @param name
     */
    String normalizeJobName(String name) {
        final result = name.replaceAll(' ','_').replaceAll(/[^a-zA-Z0-9_]/,'')
        result.size() > 128 ? result.substring(0,128) : result
    }

    /**
     * Get the Batch queue name for a task.
     *
     * @param task
     */
    String getJobQueue(TaskRun task) {
        final queue = task.config.queue?.toString()
        if( !queue )
            throw new ProcessUnrecoverableException("Missing AWS Batch job queue -- provide it by using the process `queue` directive")

        return queue
    }

    /**
     * Get the list of environment variables for a job request.
     */
    List<KeyValuePair> getEnvironmentVars() {
        def vars = [] as List<KeyValuePair>

        if( environment.containsKey('NXF_DEBUG') )
            vars << new KeyValuePair().withName('NXF_DEBUG').withValue(environment['NXF_DEBUG'])

        if( awsOptions.retryMode && awsOptions.retryMode in AwsOptions.VALID_RETRY_MODES )
            vars << new KeyValuePair().withName('AWS_RETRY_MODE').withValue(awsOptions.retryMode)

        if( awsOptions.maxTransferAttempts ) {
            vars << new KeyValuePair().withName('AWS_MAX_ATTEMPTS').withValue(awsOptions.maxTransferAttempts as String)
            vars << new KeyValuePair().withName('AWS_METADATA_SERVICE_NUM_ATTEMPTS').withValue(awsOptions.maxTransferAttempts as String)
        }

        if( fusionEnabled() ) {
            for( Map.Entry<String,String> it : fusionLauncher().fusionEnv() ) {
                vars << new KeyValuePair().withName(it.key).withValue(it.value)
            }
        }

        return vars
    }

    /**
     * Get the Batch job definition name for a task.
     *
     * @param task
     */
    String getJobDefinition(TaskRun task) {
        final container = task.getContainer()
        if( !container )
            throw new ProcessUnrecoverableException("Invalid AWS Batch job definition -- provide a Docker image name or a Batch job definition name")

        if( container.startsWith('job-definition://'))
            return container.substring(17)

        resolveJobDefinition(container)
    }

    /**
     * Get the Batch job definition name associated with a Docker container image.
     *
     * @param container
     */
    String resolveJobDefinition(String container) {
        final int DEFAULT_BACK_OFF_BASE = 3
        final int DEFAULT_BACK_OFF_DELAY = 250
        final int MAX_ATTEMPTS = 5
        int attempt = 0
        while( true ) {
            try {
                return resolveJobDefinition0(container)
            }
            catch (ClientException e) {
                if( e.statusCode != 404 || attempt++ > MAX_ATTEMPTS)
                    throw e

                final delay = (Math.pow(DEFAULT_BACK_OFF_BASE, attempt) as long) * DEFAULT_BACK_OFF_DELAY
                log.debug "Got AWS Client exception on Batch resolve job definition - message=$e.message; waiting for ${delay}ms (attempt=$attempt)"
                Thread.sleep(delay)
            }
        }
    }

    private String resolveJobDefinition0(String container) {
        final request = makeJobDefRequest(container)
        final token = request.getParameters().get('nf-token')
        final jobKey = "$container:$token".toString()
        if( jobDefinitions.containsKey(jobKey) )
            return jobDefinitions[jobKey]

        synchronized (jobDefinitions) {
            if( jobDefinitions.containsKey(jobKey) )
                return jobDefinitions[jobKey]

            def msg
            def name = findJobDefinition(request.jobDefinitionName, token)
            if( name ) {
                msg = "[AWS BATCH] Found job definition name=$name; container=$container"
            }
            else {
                name = registerJobDefinition(request)
                msg = "[AWS BATCH] Created job definition name=$name; container=$container"
            }
            // log the request
            if( log.isTraceEnabled() )
                log.debug "[AWS BATCH] $msg; request=${request.toString().indent()}"
            else
                log.debug "[AWS BATCH] $msg"

            jobDefinitions[jobKey] = name
            return name
        }
    }

    /**
     * Create a Batch job definition request for a container image.
     *
     * @param container
     */
    RegisterJobDefinitionRequest makeJobDefRequest(String container) {
        // create the job definition request
        final hashingTokens = new ArrayList()
        final request = makeJobDefRequest0(container, hashingTokens)

        // create a unique id for the job definition
        final hash = computeUniqueToken(hashingTokens)
        request.setParameters(['nf-token': hash])

        return request
    }

    /**
     * Create a unique id from a list of hashable tokens that represent
     * a unique object.
     *
     * @param uniq
     */
    private String computeUniqueToken(List uniq) {
        return CacheHelper.hasher(uniq).hash().toString()
    }

    /**
     * Create a job definition request for a container image.
     *
     * The hashing tokens collect values that should be used to create a
     * unique job definition id for the job request.
     *
     * @param image
     * @param hashingTokens
     */
    private RegisterJobDefinitionRequest makeJobDefRequest0(String image, List hashingTokens) {
        final name = normalizeJobDefinitionName(image)
        final opts = getAwsOptions()

        final request = new RegisterJobDefinitionRequest()
        request.setJobDefinitionName(name)
        request.setType(JobDefinitionType.Container)

        // apply the container options from the task configuration
        final containerOpts = task.getConfig().getContainerOptionsMap()
        final container = createContainerProperties(containerOpts)

        // set the container configuration
        // the actual command, cpus, and memory are overridden when the job is executed
        container
            .withImage(image)
            .withCommand('true')
            .withResourceRequirements(
                new ResourceRequirement().withType(ResourceType.VCPU).withValue('1'),
                new ResourceRequirement().withType(ResourceType.MEMORY).withValue('1024')
            )

        // set the job role
        final jobRole = opts.getJobRole()
        if( jobRole )
            container.setJobRoleArn(jobRole)

        // set the logs group
        final logsGroup = opts.getLogsGroup()
        if( logsGroup )
            container.setLogConfiguration(getLogConfiguration(logsGroup, opts.getRegion()))

        // set privilged mode if fusion is enabled
        if( fusionEnabled() )
            container.setPrivileged(true)

        // add the volume mounts
        final mountsMap = new LinkedHashMap( 10)
        if( opts.cliPath ) {
            def path = Paths.get(opts.cliPath).parent.parent.toString()
            mountsMap.put('aws-cli', "$path:$path:ro")
        }

        int c = 0
        final volumes = opts.getVolumes()
        for( String vol : volumes ) {
            mountsMap.put("vol-${++c}", vol)
        }

        if( mountsMap )
            addVolumeMountsToContainer(mountsMap, container)

        // set the container options
        request.setContainerProperties(container)

        // add to the hashing tokens all values that contribute to the
        // uniqueness of the job definition
        hashingTokens.add(name)
        hashingTokens.add(container.toString())
        if( containerOpts )
            hashingTokens.add(containerOpts)

        return request
    }

    /**
     * Normalize a name to be compliant with the Batch job definition format.
     *
     * @param name
     */
    String normalizeJobDefinitionName(String name) {
        if( !name )
            return null
        if( !ContainerNameValidator.isValidImageName(name) )
            throw new IllegalArgumentException("Invalid container image name: $name")

        def result = name.replaceAll(/[^a-zA-Z0-9\-_]+/,'-')
        // Batch job definition length cannot exceed 128 characters
        // take first 40 chars + a unique MD5 hash (32 chars)
        if( result.length() > 125 ) {
            result = result.substring(0,40) + '-' + name.md5()
        }

        return "nf-" + result
    }

    @Memoized 
    LogConfiguration getLogConfiguration(String name, String region) {
        new LogConfiguration()
            .withLogDriver('awslogs')
            .withOptions([
                'awslogs-region': region,
                'awslogs-group': name
            ])
    }

    void addVolumeMountsToContainer(Map<String,String> mountsMap, ContainerProperties container) {
        final mounts = new ArrayList<MountPoint>(mountsMap.size())
        final volumes = new ArrayList<Volume>(mountsMap.size())
        for( Map.Entry<String,String> entry : mountsMap.entrySet() ) {
            final mountName = entry.key
            final parts = entry.value.tokenize(':')
            final containerPath = parts[0]
            final hostPath = parts.size() > 1 ? parts[1] : containerPath
            final readOnly = parts.size() > 2 ? parts[2]=='ro' : false
            if( parts.size() > 3 )
                throw new IllegalArgumentException("Not a valid volume mount syntax: $entry.value")

            def mount = new MountPoint()
                    .withSourceVolume(mountName)
                    .withContainerPath(hostPath)
                    .withReadOnly(readOnly)
            mounts << mount

            def vol = new Volume()
                    .withName(mountName)
                    .withHost(new Host()
                    .withSourcePath(containerPath))
            volumes << vol
        }

        if( mountsMap ) {
            container.setMountPoints(mounts)
            container.setVolumes(volumes)
        }
    }

    /**
     * Search for a Batch job definition in ACTIVE status for the given name and job definition id.
     *
     * @param name
     * @param jobId
     * @return The fully qualified Batch job definition name, e.g. {@code my-job-definition:3}
     */
    String findJobDefinition(String name, String jobId) {
        log.trace "[AWS BATCH] checking job definition with name=$name; jobid=$jobId"
        final request = new DescribeJobDefinitionsRequest().withJobDefinitionName(name)
        // bypass the proxy because this method is invoked during a
        // job submit request that's already in a separate thread pool request
        // therefore it's private by a TooManyRequestsException
        final response = describeJobDefinitions0(bypassProxy(client), request)
        final jobs = response.getJobDefinitions()
        if( jobs.size() == 0 )
            return null

        def job = jobs.find { JobDefinition it -> it.status == 'ACTIVE' && it.parameters?.'nf-token' == jobId }
        return job ? "${name}:${job.revision}" : null
    }

    static private DescribeJobDefinitionsResult describeJobDefinitions0(AWSBatch client, DescribeJobDefinitionsRequest req) {
        try {
            client.describeJobDefinitions(req)
        }
        catch( AWSBatchException e ) {
            if( e.statusCode >= 500 )
                // raise a process exception so that nextflow can try to recover it
                throw new ProcessSubmitException("Failed to describe job definitions: ${req.jobDefinitions} - Reason: ${e.errorCode}", e)
            else
                // status code < 500 are not expected to be recoverable, just throw it again
                throw e
        }
    }

    /**
     * Register a new Batch job definition.
     *
     * @param request
     * @return The fully qualified Batch job definition name, e.g. {@code my-job-definition:3}
     */
    String registerJobDefinition(RegisterJobDefinitionRequest request) {
        // add nextflow tags
        request.addTagsEntry('nextflow.io/createdAt', Instant.now().toString())
        request.addTagsEntry('nextflow.io/version', Const.APP_VER)
        // create the job def
        final res = registerJobDefinition0(bypassProxy(client), request) // bypass the client proxy! see #1024
        return "${res.jobDefinitionName}:$res.revision"
    }

    static private RegisterJobDefinitionResult registerJobDefinition0(AWSBatch client, RegisterJobDefinitionRequest req) {
        try {
            return client.registerJobDefinition(req)
        }
        catch( AWSBatchException e ) {
            if( e.statusCode >= 500 )
                // raise a process exception so that nextflow can try to recover it
                throw new ProcessSubmitException("Failed to register job definition: ${req.jobDefinitionName} - Reason: ${e.errorCode}", e)
            else
                // status code < 500 are not expected to be recoverable, just throw it again
                throw e
        }
    }

    static SubmitJobResult submitJobRequest(AWSBatch client, SubmitJobRequest request) {
        try {
            return client.submitJob(request)
        }
        catch( AWSBatchException e ) {
            if( e.statusCode >= 500 )
                // raise a process exception so that nextflow can try to recover it
                throw new ProcessSubmitException("Failed to submit job: ${request.jobName} - Reason: ${e.errorCode}", e)
            else
                // status code < 500 are not expected to be recoverable, just throw it again
                throw e
        }
    }

}
