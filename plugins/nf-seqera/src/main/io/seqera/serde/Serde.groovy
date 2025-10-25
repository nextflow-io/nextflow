/*
 * Copyright 2024, Seqera Labs
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
 *
 */

package io.seqera.serde


import java.time.Instant

import com.squareup.moshi.Moshi
/**
 * Implement model serialization-deserialization helpers
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class Serde {

    static Boolean asBoolean(Object value) {
        if( value==null )
            return null;
        if( value instanceof Boolean )
            return (Boolean) value;
        return Boolean.valueOf(value.toString());
    }

    static Integer asInteger(Object value) {
        if( value==null )
            return null;
        if( value instanceof Number )
            return ((Number) value).intValue();
        if( value instanceof CharSequence )
            return Integer.parseInt(value.toString());
        throw new IllegalArgumentException("Illegal Integer value: " + value);
    }

    static Long asLong(Object value) {
        if( value==null )
            return null;
        if( value instanceof Number )
            return ((Number) value).longValue();
        if( value instanceof CharSequence )
            return Long.parseLong(value.toString());
        throw new IllegalArgumentException("Illegal Long value: " + value);
    }

    static Instant asInstant(Object value) {
        if( value==null )
            return null;
        if( value instanceof Instant )
            return (Instant) value;
        if( value instanceof CharSequence )
            return Instant.parse(value.toString());
        throw new IllegalArgumentException("Illegal Instant value: " + value);
    }


    static <T> T fromJson(String json, Class<T> type) {
        Moshi moshi = new Moshi.Builder().build();
        try {
            return moshi.adapter(type).fromJson(json);
        } catch (IOException e) {
            final String msg = String.format("Unable to parse JSON to %s - offending value: %s", type, json);
            throw new RuntimeException(msg, e);
        }
    }

}
