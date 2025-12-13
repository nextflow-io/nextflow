# Kubernetes plugin for Nextflow

This plugin provides native Kubernetes execution capability for Nextflow pipelines, with support for pod management, volume mounting, and resource allocation.

## Get Started

To use this plugin, add it to your `nextflow.config`:

```groovy
plugins {
    id 'nf-k8s'
}
```

Configure the Kubernetes executor:

```groovy
process.executor = 'k8s'

k8s {
    namespace = 'default'
    serviceAccount = 'nextflow'
    storageClaimName = 'nextflow-pvc'
}
```

The plugin automatically detects the Kubernetes configuration from:
- In-cluster configuration (when running inside a pod)
- `~/.kube/config` file
- `KUBECONFIG` environment variable

## Examples

### Basic Kubernetes Configuration

```groovy
plugins {
    id 'nf-k8s'
}

process.executor = 'k8s'

k8s {
    namespace = 'nextflow'
    serviceAccount = 'nextflow-sa'
    storageClaimName = 'nf-workdir-pvc'
    storageMountPath = '/workspace'
}

workDir = '/workspace/work'
```

### Pod Configuration

```groovy
k8s {
    namespace = 'nextflow'
    pod = [
        [volumeClaim: 'data-pvc', mountPath: '/data'],
        [secret: 'aws-credentials', mountPath: '/root/.aws']
    ]
}
```

### Resource Requests

```groovy
process {
    executor = 'k8s'
    cpus = 2
    memory = '4 GB'

    pod = [[label: 'app', value: 'nextflow']]
}
```

## Resources

- [Kubernetes Executor Documentation](https://nextflow.io/docs/latest/kubernetes.html)

## License

[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)
