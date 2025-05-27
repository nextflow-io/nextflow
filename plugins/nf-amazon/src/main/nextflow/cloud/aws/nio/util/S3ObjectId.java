package nextflow.cloud.aws.nio.util;

public class S3ObjectId {
    private final String bucket;
    private final String key;
    private final String versionId;

    public S3ObjectId(String bucket, String key, String versionId) {
        this.bucket = bucket;
        this.key = key;
        this.versionId = versionId;
    }

    public S3ObjectId(String bucket, String key) {
        this(bucket, key, null);
    }

    public String bucket() {
        return bucket;
    }

    public String key() {
        return key;
    }

    public String versionId() {
        return versionId;
    }
}
