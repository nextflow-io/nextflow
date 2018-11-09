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

package nextflow.exception

import groovy.transform.CompileStatic

/**
 *  Exception thrown when a task guard condition e.g. {@code onlyIf} or {@code ignoreIf} raise an error
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class FailedGuardException extends Exception {

    private String source

    FailedGuardException(String message, String source, Throwable cause )  {
        super(message, cause)
        this.source = source
    }

    String getSource() { source }

}
