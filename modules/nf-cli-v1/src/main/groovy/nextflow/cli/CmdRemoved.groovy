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
import groovy.transform.CompileStatic
import nextflow.exception.AbortOperationException

/**
 * Base class for the stub of a removed command, which accepts any
 * arguments and reports the removal. Subclasses should omit the command
 * description, so that the command is hidden from the command listing.
 */
@CompileStatic
abstract class CmdRemoved extends CmdBase implements UsageAware {

    @Parameter
    List<String> args

    abstract String getMessage()

    @Override
    void run() {
        abort()
    }

    @Override
    void usage() {
        abort()
    }

    @Override
    void usage(List<String> args) {
        abort()
    }

    private void abort() {
        throw new AbortOperationException(getMessage())
    }
}
