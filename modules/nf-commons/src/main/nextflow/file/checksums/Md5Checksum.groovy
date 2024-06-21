package nextflow.file.checksums

import java.security.MessageDigest

class Md5Checksum implements BytesChecksum {

    private MessageDigest digest

    private MessageDigest digestLastMarked

    Md5Checksum() {
        this.digest = getDigest()
    }

    @Override
    void update(int b) {
        digest.update((byte) b)
    }

    @Override
    void update(byte[] b, int off, int len) {
        digest.update(b, off, len)
    }

    @Override
    long getValue() {
        throw new UnsupportedOperationException("Use getChecksumBytes() instead.")
    }

    @Override
    byte[] getChecksumBytes() {
        return digest.digest()
    }

    @Override
    void reset() {
        digest = (digestLastMarked == null)
            // This is necessary so that should there be a reset without a
            // preceding mark, the MD5 would still be computed correctly.
            ? getDigest()
            : cloneFrom(digestLastMarked)
    }

    private static MessageDigest getDigest() {
        try {
            return MessageDigest.getInstance("MD5")
        } catch (Exception e) {
            throw new IllegalStateException("Unexpected error creating MD5 checksum", e)
        }
    }

    private static MessageDigest cloneFrom(MessageDigest from) {
        try {
            return (MessageDigest) from.clone()
        } catch (CloneNotSupportedException e) { // should never occur
            throw new IllegalStateException("unexpected", e)
        }
    }
}
