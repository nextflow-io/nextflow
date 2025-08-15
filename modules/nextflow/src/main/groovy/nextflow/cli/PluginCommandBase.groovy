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

import java.nio.file.Paths

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.config.ConfigBuilder
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Abstract base class for plugin commands that want to register as top-level CLI commands.
 * 
 * This class provides the foundation for plugin commands to integrate seamlessly with
 * Nextflow's CLI system, including session management and configuration handling.
 * 
 * Plugin commands extending this class will be discoverable and executable as first-class
 * CLI commands (e.g., 'nextflow wave' instead of 'nextflow plugin nf-wave:command').
 *
 * @author Edmund Miller <edmund@seqera.io>
 */
@CompileStatic
abstract class PluginCommandBase extends CmdBase implements CommandExtensionPoint {

    private static Logger log = LoggerFactory.getLogger(PluginCommandBase)

    protected Session session
    protected String pluginId

    /**
     * Constructor that requires the plugin ID for context
     * 
     * @param pluginId The ID of the plugin providing this command
     */
    PluginCommandBase(String pluginId) {
        this.pluginId = pluginId
    }

    @Override
    String getName() {
        return getCommandName()
    }

    @Override
    CmdBase createCommand() {
        return this
    }

    @Override
    final void run() {
        log.debug("Executing plugin command: ${getCommandName()} from plugin: ${pluginId}")
        
        try {
            // Create session for plugin command
            initializeSession()
            
            // Execute the actual command logic
            execute()
        }
        catch (Exception e) {
            log.error("Error executing plugin command ${getCommandName()}: ${e.message}", e)
            throw e
        }
        finally {
            // Clean up session
            destroySession()
        }
    }

    /**
     * Initialize the Nextflow session for the plugin command.
     * This mirrors the session creation pattern from PluginAbstractExec.
     */
    protected void initializeSession() {
        if (!session) {
            final config = new ConfigBuilder()
                    .setOptions(getLauncher().options)
                    .setBaseDir(Paths.get('.'))
                    .build()
                    
            this.session = new Session(config)
            log.debug("Session initialized for plugin command: ${getCommandName()}")
        }
    }

    /**
     * Clean up the session when command execution is complete.
     */
    protected void destroySession() {
        if (session) {
            session.destroy()
            session = null
            log.debug("Session destroyed for plugin command: ${getCommandName()}")
        }
    }

    /**
     * Get the current session. Session will be initialized if not already created.
     * 
     * @return The current Nextflow session
     */
    protected Session getSession() {
        if (!session) {
            initializeSession()
        }
        return session
    }

    /**
     * Get the plugin ID that provides this command
     * 
     * @return The plugin identifier
     */
    String getPluginId() {
        return pluginId
    }

    /**
     * Abstract method that subclasses must implement to define their command logic.
     * This is called within the session lifecycle management.
     */
    abstract protected void execute()

    /**
     * Default priority for plugin commands. Can be overridden by subclasses.
     */
    @Override
    int getPriority() {
        return 100
    }

}