# Seqera Platform plugin for Nextflow

## Summary

The Seqera Platform plugin provides integration with Seqera Platform (formerly Tower) for workflow monitoring, resource tracking, and centralized pipeline management.

## Get Started

To use this plugin, add it to your `nextflow.config`:

```groovy
plugins {
    id 'nf-tower'
}
```

Configure your Seqera Platform access token:

```groovy
tower {
    enabled = true
    accessToken = '<YOUR ACCESS TOKEN>'
}
```

Alternatively, set the access token via environment variable:

```bash
export TOWER_ACCESS_TOKEN='<YOUR ACCESS TOKEN>'
```

Then enable Tower in your config:

```groovy
tower {
    enabled = true
}
```

## Examples

### Basic Configuration

```groovy
plugins {
    id 'nf-tower'
}

tower {
    enabled = true
    accessToken = '<SEQERA ACCESS TOKEN>'
}
```

### Using a Specific Workspace

```groovy
tower {
    enabled = true
    accessToken = '<SEQERA ACCESS TOKEN>'
    workspaceId = '1234567890'
}
```

### Custom Seqera Platform Endpoint

```groovy
tower {
    enabled = true
    accessToken = '<SEQERA ACCESS TOKEN>'
    endpoint = 'https://tower.mycompany.com/api'
}
```

### Reports and cache management

```groovy
tower {
    enabled = true
    reports {
        enabled = true
    }
}
```

## Resources

- [Seqera Platform Documentation](https://docs.seqera.io/)
- [Nextflow Tower Configuration](https://nextflow.io/docs/latest/tracing.html#seqera-platform)

## License

[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)
