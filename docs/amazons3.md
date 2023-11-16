(amazons3-page)=

# Amazon S3 storage

Nextflow includes support for Amazon S3 storage. Files stored in an S3 bucket can be accessed transparently in your pipeline script like any other file in the local file system.

## S3 path

In order to access an S3 file, you only need to prefix the file path with the `s3` schema and the `bucket` name where it is stored.

For example, if you need to access the file `/data/sequences.fa` stored in a bucket named `my-bucket`, that file can be accessed using the following fully qualified path:

```
s3://my-bucket/data/sequences.fa
```

The usual file operations can be applied to a path handle with the above notation. For example, the content of an S3 file can be printed as follows:

```groovy
println file('s3://my-bucket/data/sequences.fa').text
```

See the {ref}`script-file-io` section to learn more about available file operations.

## Security credentials

Amazon access credentials can be provided in two ways:

1. Using AWS access and secret keys in your pipeline configuration.
2. Using IAM roles to grant access to S3 storage on Amazon EC2 instances.

### AWS access and secret keys

The AWS access and secret keys can be specified by using the `aws` section in the `nextflow.config` configuration file as shown below:

```groovy
aws {
    accessKey = '<Your AWS access key>'
    secretKey = '<Your AWS secret key>'
    region = '<AWS region identifier>'
}
```

If the access credentials are not found in the above file, Nextflow looks for AWS credentials in the following order:

1. The `nextflow.config` file in the pipeline execution directory
2. The environment variables `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY`
3. The environment variables `AWS_ACCESS_KEY` and `AWS_SECRET_KEY`
4. The `default` profile in the AWS credentials file located at `~/.aws/credentials`
5. The `default` profile in the AWS client configuration file located at `~/.aws/config`
6. The temporary AWS credentials provided by an IAM instance role. See [IAM Roles](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/iam-roles-for-amazon-ec2.html) documentation for details.

More information regarding [AWS Security Credentials](http://docs.aws.amazon.com/general/latest/gr/aws-security-credentials.html) are available in the AWS documentation.

### IAM roles with Amazon EC2 instances

When running your pipeline in an EC2 instance, IAM roles can be used to grant access to AWS resources.

In this scenario, you only need to launch the EC2 instance with an IAM role which includes the `AmazonS3FullAccess` policy. Nextflow will detect and automatically acquire the permission to access S3 storage, without any further configuration.

Learn more about [Using IAM Roles to Delegate Permissions to Applications that Run on Amazon EC2](http://docs.aws.amazon.com/IAM/latest/UserGuide/roles-usingrole-ec2instance.html) in the Amazon documentation.

## China regions

To use an AWS China region, make sure to specify the corresponding AWS API S3 endpoint in the Nextflow configuration file as shown below:

```groovy
aws { 
    client {
        endpoint = "https://s3.cn-north-1.amazonaws.com.cn"        
    }
}
```

Read more about AWS API endpoints in the [AWS documentation](https://docs.aws.amazon.com/general/latest/gr/s3.html)

## S3-compatible storage

To use S3-compatible object storage such as [Ceph](https://ceph.io) or [Minio](https://min.io) specify the endpoint of 
your storage provider and enable the [S3 path style access](https://docs.aws.amazon.com/AmazonS3/latest/userguide/VirtualHosting.html#path-style-access) 
in your Nextflow configuration as shown below:


```groovy
aws {
    accessKey = '<Your access key>'
    secretKey = '<Your secret key>'
    client {
        endpoint = '<Your storage endpoint URL>'
        s3PathStyleAccess = true
    }
}
```

## Advanced configuration

Read {ref}`AWS configuration<config-aws>` section to learn more about advanced S3 client configuration options.
