(compute-storage-page)=

# Compute and storage

This page describes common compute and storage errors and strategies to resolve them.

(aws-compute-storage)=

## Amazon Web Services

### Job queue not found

**`JobQueue <QUEUE> not found`**

This error occurs when Nextflow cannot locate the specified AWS Batch job queue. It usually happens when job queues do not exist, are not enabled, or there is a region mismatch between the configuration and the AWS Batch environment.

To resolve this error, ensure you have defined an AWS region in your `nextflow.config` file and that it matches your Batch environment region.

### Process terminated for an unknown reason

**`Process terminated for an unknown reason -- Likely it has been terminated by the external system`**

This error typically occurs when AWS Batch is unable to execute the process script. The most common reason is that the specified Docker container image has a non-standard entrypoint that prevents the execution of the Bash launcher script required by Nextflow to run the job. Another possible cause is an issue with the AWS CLI failing to run correctly within the job environment.

To resolve this error, ensure the Docker container image used for the job does not have a custom entrypoint overriding or preventing Bash from launching and that the AWS CLI is properly installed.

Check the following logs for more detailed error information:

- `.nextflow.log` file
- Job execution log in the AWS Batch dashboard
- CloudWatch logs found in the `/aws/batch/job` log group

### Process stalled in RUNNABLE status

If a process execution is stalled in the RUNNABLE status you may see an output similar to the following:

```
executor >  awsbatch (1)
process > <PROCESS> (1) [  0%] 0 of ....
```

This error occurs when a job remains stuck in the RUNNABLE state in AWS Batch and never progresses to execution. In the AWS Console, the job will be listed as RUNNABLE indefinitely, indicating that it’s waiting to be scheduled but cannot proceed. The root cause is often related to issues with the Compute Environment, Docker configuration, or network settings.

See [Why is my AWS Batch job stuck in RUNNABLE status?](https://repost.aws/knowledge-center/batch-job-stuck-runnable-status) for several resolutions and tips to investigate this error.
