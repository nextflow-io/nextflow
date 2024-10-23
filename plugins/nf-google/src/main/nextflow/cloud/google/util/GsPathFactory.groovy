/*
 * Copyright 2019, Google Inc
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

package nextflow.cloud.google.util

import java.nio.file.Path

import com.google.api.gax.retrying.RetrySettings
import com.google.cloud.storage.StorageOptions
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.CompileStatic
import nextflow.Global
import nextflow.Session
import nextflow.cloud.google.GoogleOpts
import nextflow.cloud.google.lifesciences.GoogleLifeSciencesFileCopyStrategy
import nextflow.file.FileSystemPathFactory
/**
 * Implements FileSystemPathFactory interface for Google storage
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class GsPathFactory extends FileSystemPathFactory {

    private GoogleOpts googleOpts;

    private CloudStorageConfiguration storageConfig

    private StorageOptions storageOptions

    static private GoogleOpts getGoogleOpts() {
        final session = (Session) Global.getSession()
        if (!session)
            throw new IllegalStateException("Cannot initialize GsPathFactory: missing session")
        return GoogleOpts.fromSession(session)
    }

    static protected CloudStorageConfiguration getCloudStorageConfig(GoogleOpts googleOpts) {
        final builder = CloudStorageConfiguration.builder()
        if (googleOpts.enableRequesterPaysBuckets) {
            builder.userProject(googleOpts.getProjectId())
        }
        return builder.build()
    }

    static protected StorageOptions getCloudStorageOptions(GoogleOpts opts) {
        final transportOptions = StorageOptions.getDefaultHttpTransportOptions().toBuilder()
        if( opts.httpConnectTimeout )
            transportOptions.setConnectTimeout( (int)opts.httpConnectTimeout.toMillis() )
        if( opts.httpReadTimeout )
            transportOptions.setReadTimeout( (int)opts.httpReadTimeout.toMillis() )

        RetrySettings retrySettings =
            StorageOptions.getDefaultRetrySettings()
                .toBuilder()
                .setMaxAttempts(opts.storageOpts.retryPolicy.maxAttempts)
                .setRetryDelayMultiplier(opts.storageOpts.retryPolicy.multiplier)
                .setTotalTimeout(org.threeten.bp.Duration.ofSeconds(opts.storageOpts.retryPolicy.maxDelaySecs()))
                .build()

        return StorageOptions.getDefaultInstance()
            .toBuilder()
            .setTransportOptions(transportOptions.build())
            .setRetrySettings(retrySettings)
            .build()
    }

    private void init() {
        synchronized (this) {
            if( googleOpts!=null )
                return
            this.googleOpts = getGoogleOpts()
            this.storageConfig = getCloudStorageConfig(googleOpts)
            this.storageOptions = getCloudStorageOptions(googleOpts)
        }
    }

    @Override
    protected Path parseUri(String uri) {
        if( !uri.startsWith('gs://') )
            return null
        init()
        final str = uri.substring(5)
        final p = str.indexOf('/')
        final ret = p == -1
            ? CloudStorageFileSystem.forBucket(str, storageConfig, storageOptions).getPath('')
            : CloudStorageFileSystem.forBucket(str.substring(0,p), storageConfig, storageOptions).getPath(str.substring(p))
        return ret.normalize()
    }

    @Override
    protected String toUriString(Path path) {
        if( path instanceof CloudStoragePath ) {
            return "gs://${path.bucket()}$path".toString()
        }
        return null
    }

    @Override
    protected String getBashLib(Path path) {
        if( path instanceof CloudStoragePath ) {
            return GsBashLib.fromSession( Global.session as Session )
        }
        return null
    }

    @Override
    protected String getUploadCmd(String source, Path target) {
        if( target instanceof CloudStoragePath ) {
            GoogleLifeSciencesFileCopyStrategy.uploadCmd(source,target)
        }
        return null
    }
}
