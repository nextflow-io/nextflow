/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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
