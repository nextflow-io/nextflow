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
package nextflow.lineage.serde

import com.google.gson.Gson
import com.google.gson.JsonElement
import com.google.gson.JsonObject
import com.google.gson.JsonParseException
import com.google.gson.JsonParser
import com.google.gson.TypeAdapter
import com.google.gson.reflect.TypeToken
import com.google.gson.stream.JsonReader
import com.google.gson.stream.JsonWriter

import groovy.transform.CompileStatic
import nextflow.lineage.model.v1beta1.FileOutput
import nextflow.lineage.model.v1beta1.LinModel
import nextflow.lineage.model.v1beta1.TaskOutput
import nextflow.lineage.model.v1beta1.TaskRun
import nextflow.lineage.model.v1beta1.Workflow
import nextflow.lineage.model.v1beta1.WorkflowOutput
import nextflow.lineage.model.v1beta1.WorkflowRun
import nextflow.serde.gson.RuntimeTypeAdapterFactory

/**
 * Class to serialize LiSerializable objects including the Lineage model version.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
class LinTypeAdapterFactory<T> extends RuntimeTypeAdapterFactory<T> {
    public static final String VERSION_FIELD = 'version'
    public static final String CURRENT_VERSION = LinModel.VERSION

    LinTypeAdapterFactory() {
        super(LinSerializable.class, "kind", false)
        this.registerSubtype(WorkflowRun, WorkflowRun.simpleName)
            .registerSubtype(WorkflowOutput, WorkflowOutput.simpleName)
            .registerSubtype(Workflow, Workflow.simpleName)
            .registerSubtype(TaskRun, TaskRun.simpleName)
            .registerSubtype(TaskOutput, TaskOutput.simpleName)
            .registerSubtype(FileOutput, FileOutput.simpleName)

    }

    @Override
    <R> TypeAdapter<R> create(Gson gson, TypeToken<R> type) {
        if (!LinSerializable.class.isAssignableFrom(type.rawType)) {
            return null
        }

        def delegate = super.create(gson, type as TypeToken<T>)
        if (delegate == null) {
            return null
        }

        return new TypeAdapter<R>() {
            @Override
            void write(JsonWriter out, R value) throws IOException {
                def json = delegate.toJsonTree(value)
                if (json instanceof JsonObject) {
                    json = addVersion(json)
                }
                gson.toJson(json, out)
            }

            @Override
            R read(JsonReader reader) throws IOException {
                def json = JsonParser.parseReader(reader)
                if (json instanceof JsonObject) {
                    def obj = (JsonObject) json
                    def versionEl = obj.get(VERSION_FIELD)
                    if (versionEl == null || versionEl.asString != CURRENT_VERSION) {
                        throw new JsonParseException("Invalid or missing version")
                    }
                    obj.remove(VERSION_FIELD)
                }
                return delegate.fromJsonTree(json)
            }
        }
    }

    private static JsonObject addVersion(JsonObject json){
        if( json.has(VERSION_FIELD) )
            throw new JsonParseException("object already defines a field named ${VERSION_FIELD}")

        JsonObject clone = new JsonObject();
        clone.addProperty(VERSION_FIELD, CURRENT_VERSION)
        for (Map.Entry<String, JsonElement> e : json.entrySet()) {
            clone.add(e.getKey(), e.getValue());
        }
        return clone
    }

}
