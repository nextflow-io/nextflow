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

package nextflow.cli

import org.pf4j.ExtensionPoint

/**
 * Extension point interface for plugins that want to register top-level CLI commands.
 * 
 * Plugins implementing this interface can provide commands that appear as first-class
 * citizens in the Nextflow CLI, accessible directly as 'nextflow <command>' rather than
 * through the legacy 'nextflow plugin <pluginId>:<command>' syntax.
 * 
 * Example:
 * - Instead of: nextflow plugin nf-wave:get-container
 * - Enable: nextflow wave get-container
 *
 * @author Edmund Miller <edmund@seqera.io>
 */
interface CommandExtensionPoint extends ExtensionPoint {

    /**
     * Get the primary command name that this plugin registers.
     * This will be the top-level command name (e.g., "wave", "launch").
     *
     * @return The command name that will be available at the top level
     */
    String getCommandName()

    /**
     * Get the command description for help output.
     * This description will appear in the main help listing.
     *
     * @return A brief description of what this command does
     */
    String getCommandDescription()

    /**
     * Get the priority for this command extension.
     * Higher priority commands will take precedence in case of name conflicts.
     * Built-in commands have priority 1000 by default.
     *
     * @return Priority value (higher = higher priority), defaults to 100
     */
    default int getPriority() { return 100 }

    /**
     * Create the actual command instance that will handle execution.
     * This method is called when the command needs to be instantiated.
     *
     * @return A CmdBase instance that will handle the command execution
     */
    CmdBase createCommand()

}