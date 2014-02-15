/*
 * Copyright (c) 2012, the authors.
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

package nextflow.fs.dx
import com.esotericsoftware.kryo.Kryo
import com.esotericsoftware.kryo.Serializer
import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import groovy.util.logging.Slf4j
/**
 * Handles serialization / de-serialization a {@code DxPath} object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class DxPathSerializer extends Serializer<DxPath> {

    /**
     * Serialize a {@code DxPath}
     *
     * @param kryo The {@code Kryo} instance
     * @param stream The {@code Output} stream object
     * @param obj The {@code Path} to serialize
     */
    @Override
    void write(Kryo kryo, Output stream, DxPath obj) {
        log.trace "path serialization > $obj"
        kryo.writeObject(stream, obj.fileSystem)
        kryo.writeObjectOrNull( stream, obj.@attributes, Map );
        kryo.writeObject( stream, obj.@type );
        stream.writeInt( obj.asByteArray()?.length ?: 0 )
        stream.write( obj.asByteArray() )
        stream.writeString( obj.@fileId );
    }

    /**
     * De-serialize a {@code DxPath}
     *
     * @param kryo The {@code Kryo} instance
     * @param stream The {@code Input} stream object
     * @param clazz The Class {@code Path} to de-serialize
     * @return a {@code DxPath} object
     */
    @Override
    DxPath read(Kryo kryo, Input stream, Class<DxPath> clazz) {

        DxFileSystem fs = kryo.readObject(stream,DxFileSystem)
        Map attributes = kryo.readObjectOrNull(stream, Map)
        DxPath.PathType type = kryo.readObject(stream, DxPath.PathType)
        int len = stream.readInt()
        byte[] path = stream.readBytes(len)
        String fileId = stream.readString()

        new DxPath( fs, path, type, fileId, attributes )
    }

}

/**
 * Serialize / de-serialize a {@code DxFileSystem} object
 */
class DxFileSystemSerializer extends Serializer<DxFileSystem> {

    /**
     * Serialize a {@code DxFileSystem}
     *
     * @param kryo The {@code Kryo} instance
     * @param stream The {@code Output} stream object
     * @param obj The {@code DxFileSystem} to serialize
     */
    @Override
    void write(Kryo kryo, Output output, DxFileSystem fs) {

        output.writeString(fs.contextId);
        output.writeString(fs.contextName);

    }

    /**
     * De-serialize a {@code DxFileSystem} object
     *
     * @param kryo The {@code Kryo} instance
     * @param stream The {@code Input} stream object
     * @param clazz The Class {@code DxFileSystem} to de-serialize
     * @return a {@code DxFileSystem} object
     */
    @Override
    DxFileSystem read(Kryo kryo, Input input, Class<DxFileSystem> type) {

        DxFileSystemProvider provider = DxFileSystemProvider.defaultInstance();
        String id = input.readString()
        String name = input.readString()
        provider.getOrCreateFileSystem(id,name)

    }
}
