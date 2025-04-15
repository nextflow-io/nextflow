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

package nextflow.data.cid

import nextflow.data.cid.model.Annotation
import nextflow.data.cid.model.Checksum
import nextflow.data.cid.model.DataPath
import nextflow.data.cid.model.Parameter
import nextflow.data.cid.model.TaskOutputs
import nextflow.data.cid.model.TaskRun
import nextflow.data.cid.model.Workflow
import nextflow.data.cid.model.WorkflowOutputs
import nextflow.data.cid.model.WorkflowRun

import java.lang.reflect.Field

/**
 * Class to validate if the string refers to a property in the classes of te CID Metadata model.
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class CidPropertyValidator {

    private static List<Class> CID_MODEL_CLASSES = [Workflow, WorkflowRun, WorkflowOutputs, TaskRun, TaskOutputs, DataOutput, DataPath, Parameter, Checksum, Annotation]
    private Set<String> validProperties

    CidPropertyValidator(){
        this.validProperties = new HashSet<String>()
        for( Class clazz: CID_MODEL_CLASSES) {
            for( Field field: clazz.declaredFields) {
                validProperties.add( field.name)
            }
        }
    }

    void validate(Collection<String> properties) {
        for(String property: properties) {
            if (!(property in this.validProperties)) {
                throw new IllegalArgumentException("Property '$property' doesn't exist in the CID model")
            }
        }
    }

    void validateQueryParams (Map<String, String> params){
        for(String key: params.keySet()) {
           validate(key.tokenize('.'))
        }
    }




}
