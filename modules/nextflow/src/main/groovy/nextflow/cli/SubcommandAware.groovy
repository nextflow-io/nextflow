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

import groovy.transform.CompileStatic

/**
 * Implemented by commands that dispatch to nested sub-commands (e.g. {@code fs},
 * {@code secrets}, {@code module}). This lets {@link CliSchema} describe those
 * sub-commands recursively under {@code -help-json}, since the global flag is
 * resolved against the top-level command and can't drill into a sub-command on
 * its own.
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@CompileStatic
interface SubcommandAware {

    /**
     * @return the nested sub-commands of this command, in a normalised form
     *         that {@link CliSchema} can render.
     */
    List<Subcommand> getSubcommands()

    /**
     * A normalised description of a single nested sub-command.
     *
     * <p>When {@link #command} is set, the sub-command is a full
     * {@link CmdBase} and {@link CliSchema} introspects it for complete
     * option/argument detail (recursing further if it too is
     * {@code SubcommandAware}). Otherwise only {@link #name}, {@link #help}
     * and {@link #aliases} are surfaced — the right granularity for
     * sub-commands whose arguments are parsed manually rather than via
     * JCommander annotations.
     */
    static class Subcommand {
        String name
        String help
        List<String> aliases
        CmdBase command
    }

}
