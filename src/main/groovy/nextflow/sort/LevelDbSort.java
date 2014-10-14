package nextflow.sort;

import static org.iq80.leveldb.impl.Iq80DBFactory.factory;

import java.io.IOException;
import java.nio.ByteBuffer;

import nextflow.extension.FileEx;
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
        FileEx.closeQuietly(db);
        super.close();
    }
}
