/*
 * Copyright 2013-2023, Seqera Labs
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
 */

package nextflow.cloud.google.nio

import java.nio.file.FileSystem

import com.google.cloud.storage.StorageOptions
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import groovy.transform.CompileStatic
import groovy.transform.Memoized
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class GsFileSystem extends FileSystem {

    @Delegate
    CloudStorageFileSystem target

    private GsFileSystemProvider provider

    GsFileSystem(CloudStorageFileSystem target) {
        this.target = target
        this.provider = new GsFileSystemProvider(target.provider())
    }

    @Override
    GsPath getPath(String path, String... more) {
        return new GsPath(this, target.getPath(path, more))
    }

    @Override
    GsFileSystemProvider provider() {
        return provider
    }

    static GsFileSystem forBucket(String bucket) {
        final target = CloudStorageFileSystem.forBucket(bucket)
        return new GsFileSystem(target)
    }

    static GsFileSystem forBucket(String bucket, CloudStorageConfiguration config, StorageOptions storageOptions) {
        final target = CloudStorageFileSystem.forBucket(bucket, config, storageOptions)
        return new GsFileSystem(target)
    }

}
