(spot-retry-page)=

# Spot instance failures and retries

This page describes changes in how Nextflow handles spot instance failures on AWS and Google Cloud, the impact of those changes, and how to configure spot retry behavior for your pipelines. It applies to Nextflow v24.10 and later.

## Retry behavior

Up to version 24.10, Nextflow silently retried spot instance failures up to `5` times when using AWS Batch or Google Batch. These retries were controlled by cloud-specific configuration parameters (e.g., `aws.batch.maxSpotAttempts`) and happened in cloud infrastructure without explicit visibility to Nextflow.

<h3>Before Nextflow v24.10</h3>

By default, Nextflow would instruct AWS and Google to automatically retry jobs lost to spot reclamation up to `5` times. Retries were handled by the cloud provider _within_ a Nextflow task. It was often unclear that tasks were restarted as there was no explicit message. Task runtimes and associated cloud costs were increased because they included the runtime of the reclaimed and retried tasks. Due to the high likelihood of reclamation before completion, long-running tasks running on spot instances frequently required retries, leading to inefficient allocation of resources and higher costs.

<h3>After Nextflow v24.10</h3>

The default spot reclamation retry setting changed to `0` on AWS and Google. By default, no _internal_ retries are attempted on these platforms. Spot reclamations leads to an immediate failure that is exposed to Nextflow in the same way as other generic failures (for example, returning `exit code 1` on AWS). Nextflow treats these failures like any other job failure unless a retry strategy is configured.

## Impact on Existing Workflows

If you have been relying on silent spot retries (the previous default behavior), you may now see more tasks fail with the following characteristics:

- **AWS**: Generic failure with `exit code 1`. You may see messages indicating the host machine was terminated.
- **Google**: Spot reclamation typically produces a specific code, but is now surfaced as a recognizable task failure in Nextflow logs.

Since the default for spot retries is now `0`, you must actively enable a retry strategy if you want Nextflow to handle reclaimed spot instances automatically.

## Possible actions

### Do nothing

If you do not configure anything, you will observe more pipeline failures when spot instances are reclaimed. This approach provides clearer visibility into failures. Failed tasks can be re-run with the `-resume` option. However, frequent task reclamation may lead to a higher failure rate and each retry requires manual intervention.

:::{note}
If you resume the pipeline using the resume option, it will pick up at the point the pipeline was interrupted and start with a retry of that task.
:::

### Re-enable spot retries

You can re-enable spot retries at the provider level in your Nextflow configuration:

```
// nextflow.config
aws {
    batch {
        maxSpotAttempts = 5
    }
}

google {
    batch {
        maxSpotAttempts = 5
    }
}
```

The above example sets the maximum number of spot retries to `5` for both AWS and Google.

### Make spot failures visible but retry them

You can set `maxRetries` to enable Nextflow-level retries for any failure:

```
// nextflow.config
process {
    maxRetries = 5
}
```

The above example sets retries to `5` for any failures, not just failures at the provider level.

### Use Fusion snapshots

If you have long-running tasks where losing progress due to spot reclamation is costly, consider [Fusion snapshots](https://docs.seqera.io/fusion/guide/snapshots). Fusion snapshots allows you to resume a partially completed task on a new machine if a spot instance is reclaimed and reduce wasted compute time.

Key features of Fusion snapshots:

- Jobs do not need to start from scratch after reclamation.
- Especially useful for tasks that take many hours or days to complete.
- May significantly reduce costs and run times in high-reclamation environments.

See [Fusion snapshots](https://docs.seqera.io/fusion/guide/snapshots) for more information.

## Best practices

Best practices for spot instance failures and retries:

- **Evaluate job duration**: If your tasks are very long (for example, multi-hour, multi-day), spot instances can cause repeated interruptions. Consider using on-demand instances or Fusion snapshots.
- **Set sensible retry limits**: If you enable spot retries, choose a retry count that balances the cost savings of spot usage against the overhead of restarting tasks.
- **Monitor logs and exit codes**: Failures due to spot reclamation will now appear in Nextflow logs. Monitoring failures and fine-tune your strategy.
- **Consider partial usage of spot**: Some workflows may mix on-demand instances for critical or long tasks and spot instances for shorter less critical tasks. This can optimize cost while minimizing wasted compute time.
