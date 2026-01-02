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
    public static final String SPEC_FIELD = 'spec'
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
                final object = new JsonObject()
                object.addProperty(VERSION_FIELD, CURRENT_VERSION)
                String label = getLabelFromSubtype(value.class)
                if (!label)
                     throw new JsonParseException("Not registered class ${value.class}")
                object.addProperty(getTypeFieldName(), label)
                def json = gson.toJsonTree(value)
                object.add(SPEC_FIELD, json)
                gson.toJson(object, out)
            }

            @Override
            R read(JsonReader reader) throws IOException {
                final obj = JsonParser.parseReader(reader)?.getAsJsonObject()
                if( obj==null )
                    throw new JsonParseException("Parsed JSON object is null")
                final versionEl = obj.get(VERSION_FIELD)
                if (versionEl == null || versionEl.asString != CURRENT_VERSION) {
                    throw new JsonParseException("Invalid or missing '${VERSION_FIELD}' JSON property")
                }
                final typeEl = obj.get(getTypeFieldName())
                if( typeEl==null )
                    throw new JsonParseException("JSON property '${getTypeFieldName()}' not found")
                
                // Check if this is the new format (has 'spec' field) or old format (data at root level)
                final specEl = obj.get(SPEC_FIELD)?.asJsonObject
                if ( specEl != null ) {
                    // New format: data is wrapped in 'spec' field
                    specEl.add(getTypeFieldName(), typeEl)
                    return (R) delegate.fromJsonTree(specEl)
                } else {
                    // Old format: data is at root level, just remove version field
                    obj.remove(VERSION_FIELD)
                    return (R) delegate.fromJsonTree(obj)
                }
            }
        }
    }

}
