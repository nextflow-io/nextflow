/*
 * Copyright 2013-2024, Seqera Labs
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

package io.seqera.wave.plugin.cli

import groovy.transform.CompileStatic
import nextflow.cli.CmdBase
import nextflow.cli.CommandExtensionPoint

/**
 * Extension point registration for Wave command
 *
 * @author Edmund Miller <edmund@seqera.io>
 */
@CompileStatic
class WaveCmdExtension implements CommandExtensionPoint {

    @Override
    String getCommandName() {
        return 'wave'
    }

    @Override
    String getCommandDescription() {
        return 'Execute Wave container operations'
    }

    @Override
    int getPriority() {
        return 200
    }

    @Override
    CmdBase createCommand() {
        return new WaveCmd()
    }

}