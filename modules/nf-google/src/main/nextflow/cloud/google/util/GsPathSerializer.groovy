/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import com.esotericsoftware.kryo.Kryo
import com.esotericsoftware.kryo.Serializer
import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.file.FileHelper
import nextflow.util.SerializerRegistrant

/**
 * Serializer for a {@link CloudStoragePath}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class GsPathSerializer extends Serializer<CloudStoragePath> implements SerializerRegistrant {

    @Override
    void write(Kryo kryo, Output output, CloudStoragePath target) {
        def path = target.toString()
        if( !path.startsWith('/') )  // <-- it looks a bug in the google nio library, in some case the path return is not absolute
            path = '/' + path
        path = target.bucket() + path
        log.trace "Google CloudStoragePath serialisation > path=$path"
        output.writeString(path)
    }

    @Override
    CloudStoragePath read(Kryo kryo, Input input, Class<CloudStoragePath> type) {
        final path = input.readString()
        log.trace "Google CloudStoragePath de-serialization > path=$path"
        def uri = URI.create( CloudStorageFileSystem.URI_SCHEME + '://' + path )
        def fs = FileHelper.getOrCreateFileSystemFor(uri)
        (CloudStoragePath) fs.provider().getPath(uri)
    }

    @Override
    void register(Map<Class, Object> serializers) {
        serializers.put(CloudStoragePath, GsPathSerializer)
    }
}