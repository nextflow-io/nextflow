# Dataset Upload Test Workflow

This directory contains a test workflow for validating the automatic dataset upload feature to Seqera Platform.

## Prerequisites

- Nextflow 25.10.0 or later (required for workflow output feature)
- `TOWER_ACCESS_TOKEN` environment variable set
- Access to Seqera Platform (https://cloud.seqera.io)

**Note**: This workflow uses the `output {}` block feature which requires Nextflow 25.10.0+. Enable the feature flag in your config:

```groovy
nextflow.preview.output = true
```

## Running the Test

```bash
# Make sure your token is set
export TOWER_ACCESS_TOKEN="your-token-here"

# Run the workflow
nextflow run test-workflow.nf

# Or with a custom run name
nextflow run test-workflow.nf -name my-test-run
```

## What to Expect

1. **Workflow Execution**:

    - Creates a CSV file with sample processing results
    - Publishes the CSV as a workflow output named `analysis_results`

2. **Dataset Upload**:

    - On workflow completion, the plugin will:
        - Create a new dataset named `{runName}-outputs`
        - Upload the `results.csv` index file to the dataset

3. **Log Messages**:
   Look for these messages in the console output:

    ```
    [INFO] Creating new dataset: {runName}-outputs
    [INFO] Created dataset '{name}' with ID: {id}
    [INFO] Uploading index file for output 'analysis_results' to dataset {id}: ...
    [INFO] Successfully uploaded index file for output 'analysis_results' to dataset {id}
    ```

4. **Verify in Platform**:
    - Go to https://cloud.seqera.io
    - Navigate to the "Datasets" section
    - Find the dataset named `{runName}-outputs`
    - Verify the CSV file is uploaded with header row

## Configuration Options

Edit `nextflow.config` to test different configurations:

### Disable for Specific Outputs

```groovy
tower.datasets {
    perOutput {
        'analysis_results' {
            enabled = false
        }
    }
}
```

### Use Existing Dataset

```groovy
tower.datasets {
    perOutput {
        'analysis_results' {
            datasetId = 'your-dataset-id'
        }
    }
}
```

### Disable Auto-Create

```groovy
tower.datasets {
    createMode = 'existing'  // Only use existing datasets
}
```

## Troubleshooting

### No Dataset Created

- Check that `TOWER_ACCESS_TOKEN` is set correctly
- Verify `tower.datasets.enabled = true` in config
- Look for error messages in the workflow log

### Upload Failed

- Check network connectivity to api.cloud.seqera.io
- Verify workspace permissions
- Check that the index file exists in `results/` directory

### Dataset Not Found in UI

- Refresh the Datasets page in the Platform UI
- Check you're in the correct workspace
- Search by the workflow run name

## Cleanup

After testing, you can delete the test dataset from the Platform UI or using tower-cli:

```bash
tw datasets delete -n "{runName}-outputs"
```
