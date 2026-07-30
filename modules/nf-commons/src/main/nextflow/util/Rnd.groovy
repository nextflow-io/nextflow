/*
 * Copyright 2013-2026, Seqera Labs
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
