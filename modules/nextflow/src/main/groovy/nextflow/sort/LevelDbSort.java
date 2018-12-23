/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.sort;

import static org.iq80.leveldb.impl.Iq80DBFactory.factory;

import java.io.IOException;
import java.nio.ByteBuffer;

import nextflow.extension.FilesEx;
import nextflow.util.KryoHelper;
import org.iq80.leveldb.DB;
import org.iq80.leveldb.Options;

/**
 * Level-DB based big sort
 *
 * https://github.com/dain/leveldb
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class LevelDbSort<V> extends BigSort<V> {

    private DB db;

    /**
     * Creates a DB instance
     * @return The {code LevelDbSort} object itself
     * @throws IOException
     */
    public LevelDbSort<V> create() throws IOException {
        super.create();
        Options options = new Options().createIfMissing(true);
        db = factory.open(getTempDir().resolve("data").toFile(), options);
        return this;
    }

    private static byte[] bytes(long value) {
        return ByteBuffer.allocate(8).putLong(value).array();
    }

    /**
     * Put a value into the DB
     *
     * @param key   The entry key
     * @param value The entry value
     * @return The length in bytes of the object inserted (including the key)
     */
    @Override
    protected int put(long key, V value) {
        byte[] dbKey = bytes(key);
        byte[] dbValue= KryoHelper.serialize(value);
        db.put(dbKey, dbValue);
        return 8 + dbValue.length;
    }

    /**
     * Retrieve an object from the DB
     * @param key The entry key
     * @return The object instance
     */
    @Override
    @SuppressWarnings("unchecked")
    protected V get(long key) {
        byte[] dbKey = bytes(key);
        byte[] dbValue = db.get(dbKey);
        return (V)KryoHelper.deserialize(dbValue);
    }

    /**
     * Shutdown the DB and remove temporary files
     */
    public void close() {
        FilesEx.closeQuietly(db);
        super.close();
    }
}
