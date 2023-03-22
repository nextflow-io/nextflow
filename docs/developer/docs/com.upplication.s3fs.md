
# Package `com.upplication.s3fs`

The `com.upplication.s3fs` package implements the S3 filesystem.

## Class Diagram

```mermaid
--8<-- "class-diagrams/com.upplication.s3fs.mmd"
```

!!! note
    Some classes may be excluded from the above diagrams for brevity.

## Notes

The S3 filesystem translates Java Path API calls into S3 API calls, which allows Nextflow to interact with both local files and S3 objects through the same interface.
