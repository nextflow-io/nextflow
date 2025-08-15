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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.plugin.Plugins
import org.pf4j.PluginWrapper

/**
 * Utility class for discovering and managing plugin commands that register as top-level CLI commands.
 * 
 * This class works with the plugin system to find CommandExtensionPoint implementations
 * and make them available as first-class CLI commands.
 *
 * @author Edmund Miller <edmund@seqera.io>
 */
@Slf4j
@CompileStatic
class PluginCommandDiscovery {

    /**
     * Discover all plugin commands from loaded plugins
     * 
     * @return List of CommandExtensionPoint instances from all loaded plugins
     */
    static List<CommandExtensionPoint> discoverPluginCommands() {
        try {
            // Get all CommandExtensionPoint extensions from loaded plugins
            return Plugins.getExtensions(CommandExtensionPoint)
        }
        catch (Exception e) {
            log.debug("Error discovering plugin commands: ${e.message}")
            return []
        }
    }

    /**
     * Create a map of command names to their extension points, handling conflicts by priority
     * 
     * @return Map from command name to CommandExtensionPoint, with conflicts resolved by priority
     */
    static Map<String, CommandExtensionPoint> getCommandMap() {
        final commands = discoverPluginCommands()
        final commandMap = [:] as Map<String, CommandExtensionPoint>
        
        // Group commands by name and resolve conflicts by priority
        final grouped = commands.groupBy { it.commandName }
        
        grouped.each { commandName, extensionPoints ->
            if (extensionPoints.size() == 1) {
                commandMap[commandName] = extensionPoints[0]
            }
            else {
                // Multiple plugins providing the same command name - choose by priority
                final winner = extensionPoints.max { it.priority }
                commandMap[commandName] = winner
                
                // Log warning about conflict resolution
                log.warn("Command name conflict for '${commandName}': " +
                        "Choosing plugin with priority ${winner.priority}. " +
                        "Other candidates: ${extensionPoints.findAll { it != winner }.collect { "${it.class.name} (priority: ${it.priority})" }.join(', ')}")
            }
        }
        
        return commandMap
    }

    /**
     * Get all available plugin command names
     * 
     * @return Set of command names provided by plugins
     */
    static Set<String> getAvailableCommandNames() {
        return getCommandMap().keySet()
    }

    /**
     * Create a CmdBase instance for the given command name
     * 
     * @param commandName The command name to create
     * @return CmdBase instance or null if command not found
     */
    static CmdBase createCommand(String commandName) {
        final commandMap = getCommandMap()
        final extensionPoint = commandMap[commandName]
        
        if (!extensionPoint) {
            log.debug("No plugin command found for: ${commandName}")
            return null
        }
        
        try {
            final command = extensionPoint.createCommand()
            log.debug("Created plugin command: ${commandName} from ${extensionPoint.class.name}")
            return command
        }
        catch (Exception e) {
            log.error("Failed to create plugin command ${commandName}: ${e.message}", e)
            return null
        }
    }

    /**
     * Check if a command name is provided by a plugin
     * 
     * @param commandName The command name to check
     * @return true if the command is provided by a plugin
     */
    static boolean isPluginCommand(String commandName) {
        return getAvailableCommandNames().contains(commandName)
    }

    /**
     * Get command descriptions for help output
     * 
     * @return Map from command name to description for all plugin commands
     */
    static Map<String, String> getCommandDescriptions() {
        final commandMap = getCommandMap()
        final descriptions = [:] as Map<String, String>
        
        commandMap.each { commandName, extensionPoint ->
            descriptions[commandName] = extensionPoint.commandDescription
        }
        
        return descriptions
    }

    /**
     * Initialize plugin system if not already initialized.
     * This ensures plugins are loaded before command discovery.
     */
    static void ensurePluginsInitialized() {
        try {
            // Initialize plugins if not already done
            if (!Plugins.manager) {
                log.debug("Initializing plugin system for command discovery")
                Plugins.init()
            }
        }
        catch (Exception e) {
            log.debug("Plugin system initialization failed: ${e.message}")
        }
    }

}