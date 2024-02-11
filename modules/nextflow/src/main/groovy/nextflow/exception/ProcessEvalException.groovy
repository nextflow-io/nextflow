/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.exception

import groovy.transform.CompileStatic

/**
 * Exception thrown when a command output returns a non-zero exit status
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ProcessEvalException extends RuntimeException implements ShowOnlyExceptionMessage {

    String command
    String output
    int status

    ProcessEvalException(String message, String command, String output, int status) {
        super(message)
        this.command = command
        this.output = output
        this.status = status
    }
}
