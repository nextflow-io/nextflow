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

package nextflow.cloud.aws.util


import com.esotericsoftware.kryo.Kryo
import com.esotericsoftware.kryo.Serializer
import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import nextflow.cloud.aws.nio.S3Path
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
        final scheme = target.getFileSystem().provider().getScheme()
        final path = target.toString()
        log.trace "S3Path serialization > scheme: $scheme; path: $path"
        output.writeString(scheme)
        output.writeString(path)
    }

    @Override
    S3Path read(Kryo kryo, Input input, Class<S3Path> type) {
        final scheme = input.readString()
        final path = input.readString()
        if( scheme != 's3' ) throw new IllegalStateException("Unexpected scheme for S3 path -- offending value '$scheme'")
        log.trace "S3Path de-serialization > scheme: $scheme; path: $path"
        return (S3Path) S3PathFactory.create("s3://${path}")
    }
    
}
