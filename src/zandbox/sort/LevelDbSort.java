package test;

import nextflow.sort.BigSort;

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
abstract class LevelDbSort<V> extends BigSort<V> {

//    private DB db;
//
//    public LevelDbSort<V> create() throws IOException {
//        super.create();
//        Options options = new Options().createIfMissing(true);
//        db = factory.open(new File(getTempDir(), "data"), options);
//        return this;
//    }
//
//    private static byte[] bytes(long value) {
//        return ByteBuffer.allocate(8).putLong(value).array();
//    }
//
//    @Override
//    protected int put(long key, V value) {
//        byte[] dbKey = bytes(key);
//        byte[] dbValue= KryoHelper.serialize(value);
//        db.put(dbKey, dbValue);
//        return 8 + dbValue.length;
//    }
//
//    @Override
//    protected V get(long key) {
//        byte[] dbKey = bytes(key);
//        byte[] dbValue = db.get(dbKey);
//        return (V)KryoHelper.deserialize(dbValue);
//    }
//
//    public void close() {
//        FileEx.closeQuietly(db);
//        super.close();
//    }
}
