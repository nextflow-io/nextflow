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
 */

package nextflow.cli.v2

import picocli.CommandLine.Command

/**
 * CLI `self-update` sub-command (v2)
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Command(
    name = 'self-update',
    description = 'Update nextflow runtime to the latest available version'
)
class SelfUpdateCmd extends AbstractCmd {
    @Override
    void run() {
        // actually it's doing nothing, the update process is managed by the external launcher script
        // this class is only necessary to print the usage help
    }
}
