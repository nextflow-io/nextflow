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

import nextflow.data.cid.CidSerializable
import nextflow.data.cid.model.DataType
import nextflow.data.cid.model.Output
import nextflow.data.cid.model.WorkflowResults
import nextflow.data.cid.model.WorkflowRun

import groovy.json.JsonOutput
import groovy.json.JsonSlurper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class JsonEncoder implements Encoder<String> {

    JsonEncoder() {}

    @Override
    String encode(CidSerializable object) {
        return JsonOutput.toJson(object)
    }

    @Override
    CidSerializable decode(String encoded) {
        final object = new JsonSlurper().parseText(encoded) as Map
        final dataType = DataType.valueOf(object.type as String)
        return dataType.clazz.newInstance(object)
    }
}
