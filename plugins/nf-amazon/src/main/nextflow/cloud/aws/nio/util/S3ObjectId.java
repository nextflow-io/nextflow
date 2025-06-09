/*
 * Copyright 2020-2025, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

package nextflow.cloud.aws.nio.util;

/**
 * Class to mimic Old V1 S3ObjectId
 */
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
