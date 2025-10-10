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
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException
import nextflow.platform.PlatformHelper
import nextflow.plugin.Plugins
import org.fusesource.jansi.Ansi
import org.pf4j.ExtensionPoint

import java.nio.file.Paths

import static org.fusesource.jansi.Ansi.*

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardOpenOption
import java.util.concurrent.CompletableFuture
import java.util.concurrent.TimeUnit
import java.util.regex.Pattern

/**
 * Implements the 'nextflow auth' commands
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Manage Seqera Platform authentication")
class CmdAuth extends CmdBase implements UsageAware {

    interface SubCmd {
        String getName()
        void apply(List<String> result)
        void usage(List<String> result)
    }

    interface AuthCommand extends ExtensionPoint {
        void login(String url)
        void logout()
        void config()
        void status()
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

    //
    // nextflow auth login
    //
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
            final builder = new ConfigBuilder().setHomeDir(Const.APP_HOME_DIR).setCurrentDir(Const.APP_HOME_DIR)
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

    //
    // nextflow auth config
    //
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
