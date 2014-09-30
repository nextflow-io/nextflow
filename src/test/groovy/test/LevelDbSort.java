package test;

import static org.iq80.leveldb.impl.Iq80DBFactory.factory;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;

import nextflow.sort.BigSort;
import nextflow.util.KryoHelper;
import org.iq80.leveldb.DB;
import org.iq80.leveldb.Options;

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class LevelDbSort<V> extends BigSort<V> {

    private DB data;

    public LevelDbSort<V> create() throws IOException {
        super.create();
        Options options = new Options().createIfMissing(true);
        data = factory.open(new File(getTempDir(), "data"), options);
        return this;
    }

    private static byte[] bytes(long value) {
        return ByteBuffer.allocate(8).putLong(value).array();
    }

    @Override
    protected int put(long key, V value) {
        byte[] dbKey = bytes(key);
        byte[] dbValue= KryoHelper.serialize(value);
        data.put(dbKey, dbValue);
        return 8 + dbValue.length;
    }

    @Override
    protected V get(long key) {
        byte[] dbKey = bytes(key);
        byte[] dbValue = data.get(dbKey);
        return (V)KryoHelper.deserialize(dbValue);
    }
}
