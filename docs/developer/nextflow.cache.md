
# `nextflow.cache`

The `nextflow.cache` package implements the cache database of previously executed tasks.

## Class Diagram

```{mermaid} diagrams/nextflow.cache.mmd
```

```{note}
Some classes may be excluded from the above diagram for brevity.
```

## Notes

The default cache database uses [Kryo](https://github.com/EsotericSoftware/kryo) to serialize and deserialize task data and [LevelDB](https://mvnrepository.com/artifact/org.iq80.leveldb/leveldb) as the storage backend. The cache database is stored in `.nextflow/cache/<session-id>` relative to the launch directory.

Each key-value pair in the cache database corresponds to a task. The key is the task hash, and the value consists of (1) the task `TraceRecord`, (2) the `TaskContext`, and (3) the task reference count.
