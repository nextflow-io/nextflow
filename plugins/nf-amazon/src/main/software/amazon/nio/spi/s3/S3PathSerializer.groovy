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

package software.amazon.nio.spi.s3


import com.esotericsoftware.kryo.Kryo
import com.esotericsoftware.kryo.Serializer
import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.util.SerializerRegistrant
import org.pf4j.Extension

/**
 * Register the S3Path serializer
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Extension
@CompileStatic
class S3PathSerializer extends Serializer<S3Path> implements SerializerRegistrant  {

    @Override
    void register(Map<Class, Object> serializers) {
        serializers.put(S3Path, S3PathSerializer)
    }

    @Override
    void write(Kryo kryo, Output output, S3Path target) {
        final uri = target.toUriString()
        log.trace "S3Path serialization > uri:$uri"
        output.writeString(uri)
    }

    @Override
    S3Path read(Kryo kryo, Input input, Class<S3Path> type) {
        final uri = input.readString()
        if( !uri.startsWith('s3') ) throw new IllegalStateException("Unexpected scheme for S3 path -- offending value '$uri'")
        log.trace "S3Path de-serialization > uri: $uri"
        return (S3Path) S3PathFactory.create(uri)
    }
    
}
