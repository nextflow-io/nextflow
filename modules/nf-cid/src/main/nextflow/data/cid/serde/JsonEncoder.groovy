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

package nextflow.data.cid.serde

import java.lang.reflect.Type

import groovy.json.JsonOutput
import groovy.json.JsonSlurper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
abstract class JsonEncoder<V> implements Encoder<String,V> {

    private Type type

    JsonEncoder() {
        this.type = TypeHelper.getGenericType(this, 0)
    }

    @Override
    String encode(V object) {
        return JsonOutput.toJson(object)
    }

    @Override
    V decode(String encoded) {
        return (V) new JsonSlurper().parseText(encoded)
    }
}
