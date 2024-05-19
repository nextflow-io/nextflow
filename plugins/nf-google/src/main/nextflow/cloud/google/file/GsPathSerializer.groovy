/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.cloud.google.file

import com.esotericsoftware.kryo.Kryo
import com.esotericsoftware.kryo.Serializer
import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.google.nio.GsPath
import nextflow.util.SerializerRegistrant
import org.pf4j.Extension
/**
 * Serializer for a {@link GsPath}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Extension
@CompileStatic
class GsPathSerializer extends Serializer<GsPath> implements SerializerRegistrant {

    @Override
    void write(Kryo kryo, Output output, GsPath target) {
        def path = target.toString()
        if( !path.startsWith('/') )  // <-- it looks a bug in the google nio library, in some case the path returned is not absolute
            path = '/' + path
        path = target.bucket() + path
        log.trace "Google GsPath serialisation > path=$path"
        output.writeString(path)
    }

    @Override
    GsPath read(Kryo kryo, Input input, Class<GsPath> type) {
        final path = input.readString()
        log.trace "Google GsPath de-serialization > path=$path"
        def uri = CloudStorageFileSystem.URI_SCHEME + '://' + path
        (GsPath) GsPathFactory.parse(uri)
    }

    @Override
    void register(Map<Class, Object> serializers) {
        serializers.put(GsPath, GsPathSerializer)
    }
}
