package nextflow.file.checksums

import java.security.MessageDigest

class Sha1Checksum implements BytesChecksum {

    private MessageDigest digest

    private MessageDigest digestLastMarked

    Sha1Checksum() {
        this.digest = getDigest()
    }

    @Override
    public void update(int b) {
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
    void reset() {
        digest = (digestLastMarked == null)
            // This is necessary so that should there be a reset without a
            // preceding mark, the Sha-1 would still be computed correctly.
            ? getDigest()
            : cloneFrom(digestLastMarked)
    }

    private MessageDigest getDigest() {
        try {
            return MessageDigest.getInstance("SHA-1")
        } catch (Exception e) {
            throw new IllegalStateException("Unexpected error creating SHA-1 checksum", e)
        }
    }

    @Override
    byte[] getChecksumBytes() {
        return digest.digest()
    }

    private MessageDigest cloneFrom(MessageDigest from) {
        try {
            return (MessageDigest) from.clone()
        } catch (CloneNotSupportedException e) { // should never occur
            throw new IllegalStateException("unexpected", e)
        }
    }
}