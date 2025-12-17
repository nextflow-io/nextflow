# Amazon Web Services plugin for Nextflow

The Amazon Web Services (AWS) plugin provides support for AWS, including AWS Batch as a compute executor, S3 as a file system, and Fusion file system for high-performance data operations.

## Get started

To use this plugin, add it to your `nextflow.config`:

```groovy
plugins {
    id 'nf-amazon'
}
```

Configure your AWS credentials using environment variables, AWS CLI profiles, or IAM roles. Then set up the executor and work directory:

```groovy
process.executor = 'awsbatch'
process.queue = '<YOUR BATCH QUEUE>'
workDir = 's3://<YOUR BUCKET>/work'

aws {
    region = 'us-east-1'
    batch {
        cliPath = '/home/ec2-user/miniconda/bin/aws'
    }
}
```

## Examples

### Basic AWS Batch configuration

```groovy
plugins {
    id 'nf-amazon'
}

process.executor = 'awsbatch'
process.queue = 'my-batch-queue'
workDir = 's3://my-bucket/work'

aws {
    region = 'eu-west-1'
    batch {
        cliPath = '/home/ec2-user/miniconda/bin/aws'
        jobRole = 'arn:aws:iam::123456789:role/MyBatchJobRole'
    }
}
```

### Using Fusion file system

```groovy
fusion {
    enabled = true
}

wave {
    enabled = true
}

process.executor = 'awsbatch'
workDir = 's3://my-bucket/work'
```

### S3 Storage Options

```groovy
aws {
    client {
        maxConnections = 20
        connectionTimeout = 10000
        storageEncryption = 'AES256'
    }
    region = 'us-east-1'
}
```

## Resources

- [AWS Batch Executor Documentation](https://nextflow.io/docs/latest/aws.html)
- [Amazon S3 Storage Documentation](https://nextflow.io/docs/latest/aws.html#s3-storage)

## License

[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)
