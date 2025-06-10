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

import java.lang.reflect.Type
import java.time.Instant
import java.time.OffsetDateTime

import com.google.gson.Gson
import com.google.gson.GsonBuilder
import com.google.gson.TypeAdapterFactory
import groovy.transform.CompileStatic
import nextflow.serde.Encoder
import nextflow.util.TypeHelper
import org.codehaus.groovy.runtime.GStringImpl

/**
 * Implement a JSON encoder based on Google Gson
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
abstract class GsonEncoder<T> implements Encoder<T, String> {

    private Type type

    private TypeAdapterFactory factory

    private boolean prettyPrint

    private boolean serializeNulls

    private volatile Gson gson

    protected GsonEncoder() {
        this.type = TypeHelper.getGenericType(this, 0)
    }

    GsonEncoder<T> withTypeAdapterFactory(TypeAdapterFactory factory) {
        this.factory = factory
        return this
    }

    GsonEncoder<T> withPrettyPrint(boolean value) {
        this.prettyPrint = value
        return this
    }

    GsonEncoder<T> withSerializeNulls(boolean value) {
        this.serializeNulls = value
        return this
    }

    private Gson gson0() {
        if( gson )
            return gson
        synchronized (this) {
            if( gson )
                return gson
            return gson = create0()
        }
    }

    private Gson create0() {
        final builder = new GsonBuilder()
        builder.registerTypeAdapter(Instant.class, new InstantAdapter())
        builder.registerTypeAdapter(OffsetDateTime.class, new OffsetDateTimeAdapter())
        builder.registerTypeAdapter(GStringImpl.class, new GStringSerializer())
        if( factory )
            builder.registerTypeAdapterFactory(factory)
        if( prettyPrint )
            builder.setPrettyPrinting()
        if( serializeNulls )
            builder.serializeNulls()
        return builder.create()
    }

    @Override
    String encode(T object) {
        return gson0().toJson(object, type)
    }

    @Override
    T decode(String json) {
        gson0().fromJson(json, type)
    }

}
