/*
 * Copyright 2023, Seqera Labs
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

import java.util.concurrent.Callable

import groovy.transform.CompileStatic
import picocli.CommandLine
import picocli.CommandLine.Command

/**
 * Base class for CLI v2 commands
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@Command(
    headerHeading = '%n',
    mixinStandardHelpOptions = true,
    abbreviateSynopsis = true,
    descriptionHeading = '%n',
    commandListHeading = '%nCommands:%n',
    requiredOptionMarker = ((char)'*'),
    usageHelpWidth = 160,
    parameterListHeading = '%nParameters:%n',
    optionListHeading = '%nOptions:%n'
)
abstract class AbstractCmd implements Callable<Integer> {

    @CommandLine.Spec
    CommandLine.Model.CommandSpec spec

    String getCliString() {
        spec.commandLine().getParseResult().originalArgs().join(' ')
    }

}
