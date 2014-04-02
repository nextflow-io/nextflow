/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.executor

import java.nio.file.Path

import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import com.hazelcast.config.GlobalSerializerConfig
import com.hazelcast.config.SerializationConfig
import com.hazelcast.config.SerializerConfig
import com.hazelcast.nio.ObjectDataInput
import com.hazelcast.nio.ObjectDataOutput
import com.hazelcast.nio.serialization.StreamSerializer
import groovy.util.logging.Slf4j
import nextflow.util.KryoHelper

/**
 * Serialize a {@code Path} by using the Kryo serializer
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class HzSerializerConfig {

    static GlobalSerializerConfig getGlobalConfig() {
        log.debug "Registering Hazelcast global serializer"
        def result = new GlobalSerializerConfig()
        result.setImplementation( new GlobalSerializer() )
        result
    }

    static SerializerConfig getMapSerializerConfig() {
        log.debug "Registering Hazelcast Map serializer"
        new SerializerConfig()
            .setTypeClass( Map )
            .setImplementation( new MapSerializer() )
    }

    static SerializerConfig getPathConfig() {
        log.debug "Registering Hazelcast Path serializer"
        new SerializerConfig()
                .setTypeClass( Path )
                .setImplementation( new PathSerializer() )

    }

    static SerializerConfig getHzCmdCallConfig() {
        log.debug "Registering Hazelcast HzCmdCall serializer"
        new SerializerConfig()
                .setTypeClass( HzCmdCall )
                .setImplementation( new HzCmdCallSerializer() )

    }

    static SerializerConfig getHzCmdResultConfig() {
        log.debug "Registering Hazelcast HzCmdNotify serializer"
        new SerializerConfig()
                .setTypeClass( HzCmdNotify )
                .setImplementation( new HzCmdNotifySerializer() )

    }

    static void registerAll( SerializationConfig config ) {

        config.addSerializerConfig( getMapSerializerConfig() )
        config.addSerializerConfig( getPathConfig() )
        config.addSerializerConfig( getHzCmdCallConfig() )
        config.addSerializerConfig( getHzCmdResultConfig() )
    }


    static class GlobalSerializer implements StreamSerializer<Object> {

        @Override
        void write(ObjectDataOutput stream, Object obj) throws IOException {
            Output output = new Output((OutputStream) stream);
            KryoHelper.kryo().writeClassAndObject(output, obj)
            output.flush();
        }

        @Override
        Object read(ObjectDataInput stream) throws IOException {
            Input input = new Input((InputStream)stream);
            return KryoHelper.kryo().readClassAndObject(input)
        }

        @Override
        int getTypeId() {
            return -1
        }

        @Override
        void destroy() {

        }
    }

    static class MapSerializer extends GlobalSerializer {
        int getTypeId() { return 101 }
    }

    static class PathSerializer extends GlobalSerializer {
        int getTypeId() { 102 }
    }

    static class HzCmdCallSerializer extends GlobalSerializer {
        int getTypeId() { 103 }
    }

    static class HzCmdNotifySerializer extends GlobalSerializer {
        int getTypeId() { 104 }
    }

}


