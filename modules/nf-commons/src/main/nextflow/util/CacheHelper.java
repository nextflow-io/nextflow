/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.util;

import com.google.common.hash.HashFunction;
import com.google.common.hash.Hasher;
import org.slf4j.LoggerFactory;

/**
 * Provide helper method to handle caching
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class CacheHelper {

    public enum HashMode {

        STANDARD, DEEP, LENIENT, SHA256;

        private static HashMode defaultValue;

        static {
            if( System.getenv().containsKey("NXF_CACHE_MODE") )
                defaultValue = valueOf(System.getenv().get("NXF_CACHE_MODE"));
        }

        public static HashMode DEFAULT() {
            return defaultValue != null ? defaultValue : STANDARD;
        }

        public static HashMode of( Object obj ) {
            if( obj==null || obj instanceof Boolean )
                return null;
            if( obj instanceof CharSequence ) {
                if( "true".equals(obj) || "false".equals(obj) )
                    return null;
                if( "standard".equals(obj) )
                    return STANDARD;
                if( "lenient".equals(obj) )
                    return LENIENT;
                if( "deep".equals(obj) )
                    return DEEP;
                if( "sha256".equals(obj) )
                    return SHA256;
            }
            LoggerFactory.getLogger(HashMode.class).warn("Unknown cache mode: {}", obj.toString());
            return null;
        }
    }

    public static Hasher hasher( Object value ) {
        return hasher(value, HashMode.STANDARD);
    }

    public static Hasher hasher( Object value, HashMode mode ) {
        return hasher( HashBuilder.defaultHasher(), value, mode );
    }

    public static Hasher hasher( HashFunction function, Object value, HashMode mode ) {
        return hasher( function.newHasher(), value, mode );
    }

    public static Hasher hasher( Hasher hasher, Object value, HashMode mode ) {
        return HashBuilder.hasher(hasher, value, mode);
    }

}
