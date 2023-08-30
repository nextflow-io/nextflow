
# `nextflow.cache`

The `nextflow.cache` package implements the cache database of previously executed tasks.

## Class Diagram

```{mermaid} diagrams/nextflow.cache.mmd
```

```{note}
Some classes may be excluded from the above diagram for brevity.
```

## Notes

The cache database uses [Kryo](https://github.com/EsotericSoftware/kryo) to serialize and deserialize task data. Each key-value pair in the cache database corresponds to a task. The key is the task hash, and the value consists of (1) the task `TraceRecord`, (2) the `TaskContext`, and (3) the task reference count.

The default cache store is backed by [LevelDB](https://mvnrepository.com/artifact/org.iq80.leveldb/leveldb) and is stored in `.nextflow/cache/<session-id>` relative to the launch directory.

The cloud cache store is backed by remote object storage such as Amazon S3, Azure Blob Storage, and Google Cloud Storage. It stores each task entry as a separate object.
