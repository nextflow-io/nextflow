/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.extension

import java.nio.file.Path

import groovy.json.JsonBuilder
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import nextflow.extension.FilesEx
import org.codehaus.groovy.runtime.InvokerHelper
/**
 * Helper converters for {@link DumpOp} operator
 *
 * @author : jorge <jorge.aguilera@seqera.io>
 */
@CompileStatic
class DumpHelper {

    static def deepReplaceToString(root, replaceNullWith = "") {
        if (root instanceof List) {
            root.collect {
                if (it instanceof Map) {
                    deepReplaceToString(it, replaceNullWith)
                } else if (it instanceof List) {
                    deepReplaceToString(it, replaceNullWith)
                } else if (it == null) {
                    replaceNullWith
                } else if (it instanceof Path) {
                    FilesEx.toUriString(it)
                } else {
                    it.toString()
                }
            }
        }
        else if (root instanceof Map) {
            root.each {
                if (it.value instanceof Map) {
                    deepReplaceToString(it.value, replaceNullWith)
                } else if (it.value instanceof List) {
                    it.value = deepReplaceToString(it.value, replaceNullWith)
                } else if (it.value == null) {
                    it.value = replaceNullWith
                } else if (it.value instanceof Path) {
                    it.value = FilesEx.toUriString((Path)it.value)
                } else {
                    it.value = it.value.toString()
                }
            }
        }
    }

    static String prettyPrint(input){
        if (input instanceof Map) {
            def converted =  deepReplaceToString(input)
            return JsonOutput.prettyPrint(JsonOutput.toJson(converted))
        } else if (input instanceof List) {
            def converted = deepReplaceToString(input)
            return new JsonBuilder(converted).toPrettyString()
        } else if (input instanceof Path) {
            return FilesEx.toUriString(input)
        } else {
            return InvokerHelper.inspect(input)
        }
    }

}
