/*
 * Copyright 2021, Microsoft Corp
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
package nextflow.cloud.azure.file

import com.esotericsoftware.kryo.Kryo
import com.esotericsoftware.kryo.Serializer
import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.azure.nio.AzPath
import nextflow.file.FileHelper
import nextflow.util.SerializerRegistrant

/**
 * Implements Serializer for {@link AzPath} objects
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AzPathSerializer extends Serializer<AzPath> implements SerializerRegistrant {

    @Override
    void write(Kryo kryo, Output output, AzPath path) {
        log.trace "Azure Blob storage path serialisation > path=$path"
        output.writeString(path.toUriString())
    }

    @Override
    AzPath read(Kryo kryo, Input input, Class<AzPath> type) {
        final path = input.readString()
        log.trace "Azure Blob storage path > path=$path"
        return (AzPath)FileHelper.asPath(path)
    }

    @Override
    void register(Map<Class, Object> serializers) {
        serializers.put(AzPath, AzPathSerializer)
    }
}
