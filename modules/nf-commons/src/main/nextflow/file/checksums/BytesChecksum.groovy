package nextflow.file.checksums

import java.util.zip.Checksum

interface BytesChecksum extends Checksum {

    /**
     * Returns the computed checksum in a byte array rather than the long provided by {@link #getValue()}.
     *
     * @return byte[] containing the checksum
     */
    byte[] getChecksumBytes()

}