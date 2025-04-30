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

package nextflow.serde.gson

import java.time.Instant

import com.google.gson.TypeAdapter
import com.google.gson.stream.JsonReader
import com.google.gson.stream.JsonToken
import com.google.gson.stream.JsonWriter
import groovy.transform.CompileStatic

/**
 * Implements a Gson adapter for {@link Instant}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class InstantAdapter extends TypeAdapter<Instant> {
    @Override
    void write(JsonWriter writer, Instant value) throws IOException {
        writer.value(value?.toString())
    }

    @Override
    Instant read(JsonReader reader) throws IOException {
        if( reader.peek() == JsonToken.NULL ) {
            reader.nextNull()
            return null
        }
        return Instant.parse(reader.nextString())
    }
}
