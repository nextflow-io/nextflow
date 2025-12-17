# AWS CodeCommit plugin for Nextflow

The AWS CodeCommit plugin provides integration with AWS CodeCommit. It enables Nextflow to pull pipeline scripts directly from CodeCommit repositories.

## Get Started

To use this plugin, add it to your `nextflow.config`:

```groovy
plugins {
    id 'nf-codecommit'
}
```

The plugin enables Nextflow to recognize CodeCommit repository URLs and authenticate using your AWS credentials.

Run a pipeline directly from CodeCommit:

```bash
nextflow run codecommit://my-repo/main.nf
```

## Examples

### Running a Pipeline from CodeCommit

```bash
nextflow run codecommit://my-pipeline-repo/main.nf
```

### Specifying a Branch or Tag

```bash
nextflow run codecommit://my-pipeline-repo/main.nf -r develop
```

### Using a Specific AWS Region

```bash
nextflow run codecommit://my-pipeline-repo/main.nf -hub codecommit -hub-opts region=eu-west-1
```

### Configuration with AWS Region

```groovy
plugins {
    id 'nf-codecommit'
}

aws {
    region = 'us-east-1'
}
```

### Using AWS Profiles

Set the AWS profile via environment variable:

```bash
export AWS_PROFILE=my-profile
nextflow run codecommit://my-repo/main.nf
```

## Resources

- [AWS CodeCommit Documentation](https://docs.aws.amazon.com/codecommit/)
- [Nextflow Pipeline Sharing](https://nextflow.io/docs/latest/sharing.html)

## License

[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)
