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

package nextflow.cli.pipeline

import java.nio.file.Path

import com.beust.jcommander.IParameterValidator
import com.beust.jcommander.ParameterException
import groovy.transform.CompileStatic
import nextflow.cli.CmdBase
import nextflow.config.ConfigCmdAdapter
import nextflow.config.RegistryConfig

/**
 * Helpers shared by the pipeline subcommands
 */
@CompileStatic
class PipelineCmdSupport {

    static final String OUTPUT_MODE_HELP = 'Output mode: text, simple (search), json'

    /** Accepts the output modes of list/view (text) and search (simple), plus json */
    static class OutputModeValidator implements IParameterValidator {

        private static final List<String> MODES = List.of('text', 'simple', 'json')

        @Override
        void validate(String name, String value) {
            if( !MODES.contains(value) )
                throw new ParameterException("Output mode must be one of $MODES (found: $value)")
        }
    }

    static RegistryConfig registryConfig(CmdBase cmd, Path baseDir) {
        final config = new ConfigCmdAdapter()
            .setOptions(cmd.launcher.options)
            .setBaseDir(baseDir)
            .build()
        return config.navigate('registry') as RegistryConfig ?: new RegistryConfig()
    }

}
