/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cli

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.exception.AbortOperationException

/**
 * Stub for the removed `console` command, kept to report a clear error.
 * Omitting the command description hides it from the command listing.
 */
@CompileStatic
@Parameters
class CmdConsole extends CmdBase {

    @Parameter
    List<String> args

    String getName() { 'console' }

    void run() {
        throw new AbortOperationException("The `console` command has been removed")
    }
}
