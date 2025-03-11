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
 *
 */

package nextflow.util

import java.time.Instant

import com.google.gson.Gson
import com.google.gson.GsonBuilder
import groovy.transform.CompileStatic
import groovy.transform.Memoized

/**
 * Implements helper for Gson ser-deserialization
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class GsonHelper {

    @Memoized
    static protected Gson gson() {
        new GsonBuilder()
            .registerTypeAdapter(Instant, new GsonInstantAdapter())
            .create()
    }

    static String toJson(Object obj) {
        return gson().toJson(obj)
    }

    static <T> T fromJson(String json, Class<T> clazz) {
        return gson().fromJson(json,clazz)
    }
}
