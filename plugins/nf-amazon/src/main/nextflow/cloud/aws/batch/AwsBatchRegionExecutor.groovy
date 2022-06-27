/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import com.amazonaws.services.batch.AWSBatch
import com.amazonaws.services.ecs.model.AccessDeniedException
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.cloud.aws.AmazonClientFactory
import nextflow.cloud.types.CloudMachineInfo
import nextflow.processor.TaskConfig
import nextflow.util.ThrottlingExecutor

import java.nio.file.Path

/**
 * AWS Batch executor
 * https://aws.amazon.com/batch/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsBatchRegionExecutor {

    AwsBatchRegionExecutor(ThrottlingExecutor submitter, ThrottlingExecutor reaper, Path remoteBinDir){
        this.submitter = submitter
        this.reaper = reaper
        this.remoteBinDir = remoteBinDir
    }

    /**
     * Reference to current running session
     */
    @PackageScope
    private Session session

    Session getSession(){
        this.session
    }

    /**
     * Proxy to throttle AWS batch client requests
     */
    @PackageScope
    private AwsBatchProxy client

    /** Helper class to resolve Batch related metadata */
    private AwsBatchHelper helper

    /**
     * executor service to throttle service requests
     */
    private ThrottlingExecutor submitter

    /**
     * Executor service to throttle cancel requests
     */
    private ThrottlingExecutor reaper

    /**
     * A S3 path where executable scripts need to be uploaded
     */
    private Path remoteBinDir = null

    private AwsOptions awsOptions

    AwsOptions getAwsOptions() {  awsOptions  }


    protected void createAwsClient(Session session, TaskConfig taskConfig) {
        this.session = session

        Map config = session.config.deepClone()
        // AmazonClientFactory is looking into top `region` key but AwsOptions into `aws.region`
        if( taskConfig?.region ) {
            config.region = taskConfig.region
            if( config.aws && config.aws instanceof Map){
                (config.aws as Map).region = taskConfig.region
            }
        }

        /*
         * retrieve config and credentials and create AWS client
         */
        final driver = new AmazonClientFactory(config)

        /*
         * create a proxy for the aws batch client that manages the request throttling
         */
        client = new AwsBatchProxy(driver.batchClient, submitter)
        helper = new AwsBatchHelper(client, driver)
        // create the options object
        awsOptions = new AwsOptions(this, config, remoteBinDir)
        log.debug "[AWS BATCH] Executor options=$awsOptions"
    }

    Path getRemoteBinDir() {
        this.remoteBinDir
    }

    AWSBatch getClient() {
        client
    }

    ThrottlingExecutor getReaper() { reaper }


    CloudMachineInfo getMachineInfoByQueueAndTaskArn(String queue, String taskArn) {
        try {
            return helper?.getCloudInfoByQueueAndTaskArn(queue, taskArn)
        }
        catch ( AccessDeniedException e ) {
            log.warn "Unable to retrieve AWS Batch instance type | ${e.message}"
            // disable it since user has not permission to access this info
            awsOptions.fetchInstanceType = false
            return null
        }
        catch( Exception e ) {
            log.warn "Unable to retrieve AWS batch instance type for queue=$queue; task=$taskArn | ${e.message}", e
            return null
        }
    }

    String getJobOutputStream(String jobId) {
        try {
            return helper.getTaskLogStream(jobId)
        }
        catch (Exception e) {
            log.debug "Unable to retrieve AWS Cloudwatch logs for Batch Job id=$jobId | ${e.message}", e
            return null
        }
    }

}









