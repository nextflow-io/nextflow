/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.k8s.client

import groovy.transform.CompileStatic
import nextflow.exception.ShowOnlyExceptionMessage

/**
 * Exception raised when a pod cannot be scheduled because
 * e.g. the container image cannot be pulled, required resources
 * cannot be fulfilled, etc.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class PodUnschedulableException extends Exception implements ShowOnlyExceptionMessage {

    PodUnschedulableException(String message, Throwable cause) {
        super(message,cause)
    }

}
