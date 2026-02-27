/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.serde

/**
 * An interface for encoding and decoding objects between two types.
 *
 * @param <T> the type of the original object to be encoded.
 * @param <S> the type of the encoded representation.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface Encoder<T,S> {

    /**
     * Encodes an object of type {@code T} into its corresponding encoded representation of type {@code S}.
     *
     * @param object the object to encode
     * @return the encoded representation of the object
     */
    S encode(T object)

    /**
     * Decodes an encoded representation of type {@code S} back into its original form of type {@code T}.
     *
     * @param encoded the encoded representation to decode
     * @return the decoded object
     */
    T decode(S encoded)

}
