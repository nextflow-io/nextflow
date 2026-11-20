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

import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cli.project.CmdProjectRemove

/**
 * Deprecated alias of `project remove`, hidden from the command listing
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Delete the local copy of a project (deprecated -- use `nextflow project remove`)")
class CmdDrop extends CmdProjectRemove {

    static final public String NAME = 'drop'

    @Override
    String getName() { NAME }

    @Override
    boolean isHidden() { true }

    @Override
    void run() {
        log.warn "The `drop` command is deprecated -- use `nextflow project remove` instead"
        super.run()
    }
}
