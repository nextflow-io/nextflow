/*
 * Copyright 2022, Google Inc.
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

package nextflow.cloud.google.batch.json

import groovy.json.JsonGenerator
import groovy.json.JsonOutput
import groovy.transform.CompileStatic

/**
 * Simple helper class to handle JSON operation
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class JsonHelper {

    private static final JsonGenerator generator

    static {
        generator = new JsonGenerator.Options()
                .excludeNulls()
                .build()
    }

    static String toJson(Object o, boolean pretty=false) {
        final json = generator.toJson(o)
        return pretty ? JsonOutput.prettyPrint(json) : json
    }

}
