/*
 * Copyright 2013-2025, Seqera Labs
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
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.SysEnv
import nextflow.config.ConfigCmdAdapter
import nextflow.exception.AbortOperationException
import nextflow.platform.PlatformHelper
import nextflow.plugin.Plugins

/**
 * Command-line interface for managing Seqera Platform authentication and configuration.
 *
 * <p>This class provides the {@code nextflow auth} command with four sub-commands for managing
 * authentication credentials and platform settings:
 *
 * <ul>
 *   <li>{@code login} - Authenticate with Seqera Platform using OAuth2 (Cloud) or PAT (Enterprise)
 *   <li>{@code logout} - Revoke access token and remove authentication from local config
 *   <li>{@code config} - Configure platform settings (workspace, monitoring, compute environment)
 *   <li>{@code status} - Display current authentication status and configuration sources
 * </ul>
 *
 * <h2>Authentication Flow</h2>
 * <p>The authentication mechanism varies based on the platform type:
 * <ul>
 *   <li><b>Seqera Platform Cloud</b>: Uses OAuth2 device authorization flow with Auth0
 *   <li><b>Seqera Platform Enterprise</b>: Prompts for Personal Access Token (PAT)
 * </ul>
 *
 * <h2>Configuration Management</h2>
 * <p>Authentication credentials and settings are stored in the user's home directory:
 * <pre>
 *   ~/.nextflow/config                 (includes seqera-auth.config)
 *   ~/.nextflow/seqera-auth.config     (contains tower.* settings)
 * </pre>
 *
 * <p>Configuration values can be sourced from:
 * <ol>
 *   <li>Nextflow config files (highest priority)
 *   <li>Environment variables (TOWER_ACCESS_TOKEN, TOWER_API_ENDPOINT, etc.)
 *   <li>Default values (lowest priority)
 * </ol>
 *
 * <h2>Plugin Architecture</h2>
 * <p>The actual authentication implementation is provided by the {@code nf-tower} plugin through
 * the {@link AuthCommand} extension point. This allows the authentication logic to be maintained
 * separately from the CLI interface.
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 * @see AuthCommand
 * @see nextflow.platform.PlatformHelper
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Manage Seqera Platform authentication")
class CmdAuth extends CmdBase implements UsageAware {

    /**
     * Interface for auth sub-commands that defines the contract for command execution and help text.
     *
     * <p>Each sub-command (login, logout, config, status) implements this interface to provide:
     * <ul>
     *   <li>Command name identification
     *   <li>Argument processing and execution logic
     *   <li>Usage/help text generation
     * </ul>
     */
    interface SubCmd {
        /**
         * @return the name of this sub-command (e.g., "login", "logout")
         */
        String getName()

        /**
         * Executes the sub-command with the provided arguments.
         *
         * @param result the command-line arguments passed to this sub-command
         */
        void apply(List<String> result)

        /**
         * Generates usage/help text for this sub-command.
         *
         * @param result the list to which usage text lines should be added
         */
        void usage(List<String> result)
    }

    static public final String NAME = 'auth'

    private List<SubCmd> commands = []

    private AuthCommand operation

    String getName() {
        return NAME
    }

    @Parameter(hidden = true)
    List<String> args

    @Parameter(names = ['-u', '-url'], description = 'Seqera Platform API endpoint')
    String apiUrl

    CmdAuth() {
        commands.add(new LoginCmd())
        commands.add(new LogoutCmd())
        commands.add(new ConfigCmd())
        commands.add(new StatusCmd())
    }

    void usage() {
        usage(args)
    }

    void usage(List<String> args) {
        List<String> result = []
        if (!args) {
            result << this.getClass().getAnnotation(Parameters).commandDescription()
            result << 'Usage: nextflow auth <sub-command> [options]'
            result << ''
            result << 'Commands:'
            result << '  login    Authenticate with Seqera Platform'
            result << '  logout   Remove authentication and revoke access token'
            result << '  status   Show current authentication status and configuration'
            result << '  config   Configure Seqera Platform settings'
            result << ''
        } else {
            def sub = commands.find { it.name == args[0] }
            if (sub)
                sub.usage(result)
            else {
                throw new AbortOperationException("Unknown auth sub-command: ${args[0]}")
            }
        }
        println result.join('\n').toString()
    }

    @Override
    void run() {
        if (!args) {
            usage()
            return
        }
        // load the Auth command implementation
        this.operation = loadOperation()
        if( !operation )
            throw new IllegalStateException("Unable to load auth extensions.")
        // consume the first argument
        getCmd(args).apply(args.drop(1))
    }

    protected AuthCommand loadOperation(){
        // setup the plugins system and load the secrets provider
        Plugins.init()
        // load the config
        Plugins.start('nf-tower')
        // get Auth command operations implementation from plugins
        return Plugins.getExtension(AuthCommand)
    }

    protected SubCmd getCmd(List<String> args) {
        def cmd = commands.find { it.name == args[0] }
        if (cmd) {
            return cmd
        }

        def matches = commands.collect { it.name }.closest(args[0])
        def msg = "Unknown auth sub-command: ${args[0]}"
        if (matches)
            msg += " -- Did you mean one of these?\n" + matches.collect { "  $it" }.join('\n')
        throw new AbortOperationException(msg)
    }

    /**
     * Implements the {@code nextflow auth login} sub-command for authenticating with Seqera Platform.
     *
     * <p>This command initiates the authentication process, which varies based on the platform type:
     * <ul>
     *   <li><b>Seqera Platform Cloud</b>: Launches OAuth2 device authorization flow
     *       <ol>
     *         <li>Requests device code from Auth0
     *         <li>Displays verification URL and code to user
     *         <li>Opens browser for authentication
     *         <li>Polls for token completion
     *         <li>Generates Personal Access Token from OAuth token
     *       </ol>
     *   </li>
     *   <li><b>Seqera Platform Enterprise</b>: Prompts for manual PAT entry
     *       <ol>
     *         <li>Displays URL for token creation
     *         <li>Prompts user to paste PAT
     *         <li>Validates token with API
     *       </ol>
     *   </li>
     * </ul>
     *
     * <p>Upon successful authentication, the command:
     * <ul>
     *   <li>Saves {@code tower.accessToken} to {@code ~/.nextflow/seqera-auth.config}
     *   <li>Configures {@code tower.endpoint} with the API URL
     *   <li>Enables {@code tower.enabled} for automatic monitoring
     *   <li>Adds {@code includeConfig 'seqera-auth.config'} to main config
     * </ul>
     *
     * <p><b>Usage:</b> {@code nextflow auth login [-u <endpoint>]}
     */
    class LoginCmd implements SubCmd {

        @Override
        String getName() { 'login' }

        @Override
        void apply(List<String> args) {
            if (args.size() > 0) {
                throw new AbortOperationException("Too many arguments for ${name} command")
            }
            operation.login(apiUrl)
        }

        @Override
        void usage(List<String> result) {
            // Read config to get the actual resolved endpoint value
            final builder = new ConfigCmdAdapter().setHomeDir(Const.APP_HOME_DIR).setCurrentDir(Const.APP_HOME_DIR)
            final config = builder.buildConfigObject().flatten()
            final towerConfig = config.findAll { it.key.toString().startsWith('tower.') }
                .collectEntries { k, v -> [(k.toString().substring(6)): v] }
            def defaultEndpoint = PlatformHelper.getEndpoint(towerConfig, SysEnv.get())

            result << 'Authenticate with Seqera Platform'
            result << "Usage: nextflow auth $name [-u <endpoint>]".toString()
            result << ''
            result << 'Options:'
            result << "  -u, -url <endpoint>    Seqera Platform API endpoint (default: ${defaultEndpoint})".toString()
            result << ''
            result << 'This command will:'
            result << '  1. Display a URL and device code for OAuth2 authentication (Cloud) or prompt for PAT (Enterprise)'
            result << '  2. Wait for user to complete authentication in web browser'
            result << '  3. Generate and save access token to home-directory Nextflow config'
            result << '  4. Configure tower.accessToken, tower.endpoint, and tower.enabled settings'
            result << ''
        }
    }

    /**
     * Implements the {@code nextflow auth logout} sub-command for removing authentication credentials.
     *
     * <p>This command performs a safe logout by:
     * <ol>
     *   <li>Checking if {@code tower.accessToken} is configured in {@code ~/.nextflow/seqera-auth.config}
     *   <li>Validating the token with the Seqera Platform API
     *   <li>Displaying current configuration details to the user
     *   <li>Prompting for confirmation before proceeding
     *   <li>Deleting the Personal Access Token from Platform (Cloud only)
     *   <li>Removing {@code seqera-auth.config} file
     *   <li>Removing {@code includeConfig} line from main config
     * </ol>
     *
     * <p><b>Important behavior differences:</b>
     * <ul>
     *   <li><b>Seqera Platform Cloud</b>: Deletes the PAT from the platform via API before removing local config
     *   <li><b>Seqera Platform Enterprise</b>: Only removes local config; does not delete PAT from platform
     * </ul>
     *
     * <p><b>Note:</b> This command only removes credentials from Nextflow config files. If
     * {@code TOWER_ACCESS_TOKEN} environment variable is set, it will remain unchanged and
     * the user will be warned about this.
     *
     * <p><b>Usage:</b> {@code nextflow auth logout}
     */
    class LogoutCmd implements SubCmd {

        @Override
        String getName() { 'logout' }

        @Override
        void apply(List<String> args) {
            if (args.size() > 0) {
                throw new AbortOperationException("Too many arguments for ${name} command")
            }
            operation.logout()
        }

        @Override
        void usage(List<String> result) {
            result << 'Log out and remove Seqera Platform authentication'
            result << "Usage: nextflow auth $name".toString()
            result << ''
            result << 'This command will:'
            result << '  1. Check if tower.accessToken is configured'
            result << '  2. Validate the token with Seqera Platform'
            result << '  3. Delete the PAT from Platform (only if Seqera Platform Cloud)'
            result << '  4. Remove the authentication from Nextflow config'
            result << ''
        }
    }

    /**
     * Implements the {@code nextflow auth config} sub-command for configuring Seqera Platform settings.
     *
     * <p>This interactive command guides users through configuring their Seqera Platform integration.
     * It requires an existing authentication (see {@link LoginCmd}) and allows configuration of:
     *
     * <h3>Workflow Monitoring</h3>
     * <ul>
     *   <li>Configures {@code tower.enabled} setting
     *   <li>When enabled: all workflow runs are automatically monitored
     *   <li>When disabled: monitoring requires {@code -with-tower} flag per run
     * </ul>
     *
     * <h3>Default Workspace</h3>
     * <ul>
     *   <li>Configures {@code tower.workspaceId} setting
     *   <li>Lists all accessible workspaces organized by organization
     *   <li>For many workspaces: uses two-stage selection (org → workspace)
     *   <li>For few workspaces: displays all at once
     *   <li>Option to use Personal workspace (no workspace ID)
     * </ul>
     *
     * <h3>Primary Compute Environment</h3>
     * <ul>
     *   <li>Displays compute environments available in selected workspace
     *   <li>Allows setting a primary compute environment
     *   <li>Primary compute environment is used by default for pipeline execution
     * </ul>
     *
     * <p>All configuration changes are saved to {@code ~/.nextflow/seqera-auth.config}.
     *
     * <p><b>Usage:</b> {@code nextflow auth config}
     */
    class ConfigCmd implements SubCmd {

        @Override
        String getName() { 'config' }

        @Override
        void apply(List<String> args) {
            if (args.size() > 0) {
                throw new AbortOperationException("Too many arguments for ${name} command")
            }

            operation.config()
        }

        @Override
        void usage(List<String> result) {
            result << 'Configure Seqera Platform settings'
            result << "Usage: nextflow auth $name".toString()
            result << ''
            result << 'This command will:'
            result << '  1. Check authentication status'
            result << '  2. Configure tower.enabled setting for workflow monitoring'
            result << '  3. Configure default workspace (tower.workspaceId)'
            result << ''
        }
    }

    /**
     * Implements the {@code nextflow auth status} sub-command for displaying authentication status.
     *
     * <p>This command provides a comprehensive overview of the current Seqera Platform configuration,
     * displaying all settings in a formatted table with their values and sources.
     *
     * <h3>Information Displayed</h3>
     * <table border="1">
     *   <tr><th>Setting</th><th>Description</th><th>Possible Sources</th></tr>
     *   <tr>
     *     <td>API endpoint</td>
     *     <td>Seqera Platform API URL</td>
     *     <td>nextflow config, env var $TOWER_API_ENDPOINT, default</td>
     *   </tr>
     *   <tr>
     *     <td>API connection</td>
     *     <td>Whether the API is reachable</td>
     *     <td>Health check result (✔ OK or ERROR)</td>
     *   </tr>
     *   <tr>
     *     <td>Authentication</td>
     *     <td>Token status and username</td>
     *     <td>nextflow config, env var $TOWER_ACCESS_TOKEN, not set</td>
     *   </tr>
     *   <tr>
     *     <td>Workflow monitoring</td>
     *     <td>Whether monitoring is enabled</td>
     *     <td>nextflow config, default (No)</td>
     *   </tr>
     *   <tr>
     *     <td>Default workspace</td>
     *     <td>Configured workspace ID and name</td>
     *     <td>nextflow config, env var $TOWER_WORKSPACE_ID, default (Personal)</td>
     *   </tr>
     *   <tr>
     *     <td>Primary compute env</td>
     *     <td>Primary compute environment name</td>
     *     <td>Workspace setting</td>
     *   </tr>
     *   <tr>
     *     <td>Default work dir</td>
     *     <td>Working directory from compute env</td>
     *     <td>Compute environment setting</td>
     *   </tr>
     * </table>
     *
     * <p>The source column shows where each setting's value comes from, making it easy to
     * understand configuration precedence and troubleshoot issues.
     *
     * <p><b>Usage:</b> {@code nextflow auth status}
     */
    class StatusCmd implements SubCmd {

        @Override
        String getName() { 'status' }

        @Override
        void apply(List<String> args) {
            if (args.size() > 0) {
                throw new AbortOperationException("Too many arguments for ${name} command")
            }
            operation.status()
        }

        @Override
        void usage(List<String> result) {
            result << 'Show authentication status and configuration'
            result << "Usage: nextflow auth $name".toString()
            result << ''
            result << 'This command shows:'
            result << '  - Authentication status (yes/no) and source'
            result << '  - API endpoint and source'
            result << '  - Monitoring enabled status and source'
            result << '  - Default workspace and source'
            result << '  - System health status (API connection and authentication)'
            result << ''
        }
    }
}
