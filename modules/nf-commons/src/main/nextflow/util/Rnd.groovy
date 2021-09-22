package nextflow.util

import java.security.SecureRandom

/**
 * Rnd key generator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class Rnd {

    private static final BigInteger B = BigInteger.ONE.shiftLeft(80); // 2^64

    private static final SecureRandom random = new SecureRandom();

    private Rnd() {}

    /**
     * Generate random key base62 encoded ie. containing only [0-9a-zA-Z] characters
     * guaranteed to be unique.
     *
     * Tested with 100 mln iteration -> 0 collision
     *
     * @return A random generated alphanumeric between 9-14 characters
     */
    static String key() {
        byte[] buffer = new byte[10]
        random.nextBytes(buffer)
        def big = new BigInteger(buffer)
        if (big.signum() < 0) {
            big = big.add(B)
        }

        return Base62.encode(big)
    }

    static String hex() {
        byte[] buffer = new byte[10]
        random.nextBytes(buffer)
        return buffer.encodeHex()
    }

}
