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
import nextflow.exception.AbortOperationException
import org.fusesource.jansi.Ansi
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

    static public final String NAME = 'auth'

    private List<SubCmd> commands = []

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
            commands.collect { it.name }.sort().each { result << "  $it".toString() }
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

        try {
            def cmd = getCmd(args)
            if (cmd instanceof LoginCmd && apiUrl) {
                cmd.apiUrl = apiUrl
            }
            cmd.apply(args.drop(1))
        } catch (Exception e) {
            throw new AbortOperationException(e.message)
        }
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


    private String promptForApiUrl() {
        System.out.print("Seqera Platform API endpoint [Default https://api.cloud.seqera.io]: ")
        System.out.flush()

        def reader = new BufferedReader(new InputStreamReader(System.in))
        def input = reader.readLine()?.trim()

        return input?.isEmpty() || input == null ? 'https://api.cloud.seqera.io' : input
    }

    private Map getCloudEndpointInfo(String apiUrl) {
        // Check production endpoints
        if (apiUrl == 'https://api.cloud.seqera.io' || apiUrl == 'https://cloud.seqera.io/api') {
            return [isCloud: true, environment: 'prod']
        }
        // Check staging endpoints
        if (apiUrl == 'https://api.cloud.stage-seqera.io' || apiUrl == 'https://cloud.stage-seqera.io/api') {
            return [isCloud: true, environment: 'stage']
        }
        // Check development endpoints
        if (apiUrl == 'https://api.cloud.dev-seqera.io' || apiUrl == 'https://cloud.dev-seqera.io/api') {
            return [isCloud: true, environment: 'dev']
        }
        // Enterprise/other endpoints
        return [isCloud: false, environment: null]
    }

    private boolean isCloudEndpoint(String apiUrl) {
        return getCloudEndpointInfo(apiUrl).isCloud
    }

    // Get user info from Seqera Platform
    private Map callUserInfoApi(String accessToken, String apiUrl) {
        def userInfoUrl = "${apiUrl}/user-info"
        def connection = new URL(userInfoUrl).openConnection() as HttpURLConnection
        connection.requestMethod = 'GET'
        connection.setRequestProperty('Authorization', "Bearer ${accessToken}")

        if (connection.responseCode != 200) {
            def error = connection.errorStream?.text ?: "HTTP ${connection.responseCode}"
            throw new RuntimeException("Failed to get user info: ${error}")
        }

        def response = connection.inputStream.text
        def json = new groovy.json.JsonSlurper().parseText(response) as Map
        return json.user as Map
    }

    private Path getConfigFile() {
        return Const.APP_HOME_DIR.resolve('config')
    }

    private Map readConfig() {
        def configFile = getConfigFile()
        if (!Files.exists(configFile)) {
            return [:]
        }

        try {
            def configText = Files.readString(configFile)
            def config = new ConfigSlurper().parse(configText)
            return config.flatten()
        } catch (Exception e) {
            throw new RuntimeException("Failed to read config file ${configFile}: ${e.message}")
        }
    }

    private String cleanTowerConfig(String content) {
        // Remove tower scoped blocks: tower { ... }
        content = content.replaceAll(/(?ms)tower\s*\{.*?\}/, '')
        // Remove individual tower.* lines
        content = content.replaceAll(/(?m)^tower\..*$\n?/, '')
        // Remove Seqera Platform configuration comment
        content = content.replaceAll(/\/\/\s*Seqera Platform configuration\s*/, '')
        // Clean up extra whitespace
        return content.replaceAll(/\n\n+/, '\n\n').trim() + "\n\n"
    }

    private void writeConfig(Map config) {
        def configFile = getConfigFile()

        // Create directory if it doesn't exist
        if (!Files.exists(configFile.parent)) {
            Files.createDirectories(configFile.parent)
        }

        // Read existing config and clean out old tower blocks
        def configText = new StringBuilder()
        if (Files.exists(configFile)) {
            def existingContent = Files.readString(configFile)
            def cleanedContent = cleanTowerConfig(existingContent)
            configText.append(cleanedContent)
        }

        // Write tower config block
        def towerConfig = config.findAll { key, value ->
            key.toString().startsWith('tower.')
        }

        configText.append("// Seqera Platform configuration\n")
        configText.append("tower {\n")
        towerConfig.each { key, value ->
            def configKey = key.toString().substring(6) // Remove "tower." prefix
            if (configKey.endsWith('.comment')) {
                // Skip comment keys - they're handled below
                return
            }

            if (value instanceof String) {
                configText.append("    ${configKey} = '${value}'\n")
            } else {
                configText.append("    ${configKey} = ${value}\n")
            }
        }
        configText.append("}\n")

        def finalConfig = configText.toString()

        // Add workspace comment if available
        if (config.containsKey('tower.workspaceId.comment')) {
            def workspaceId = config['tower.workspaceId']
            def comment = config['tower.workspaceId.comment']
            finalConfig = finalConfig.replaceAll(
                /workspaceId = '${Pattern.quote(workspaceId.toString())}'/,
                "workspaceId = '${workspaceId}'  // ${comment}"
            )
        }

        Files.writeString(configFile, finalConfig, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)
    }

    //
    // nextflow auth login
    //
    class LoginCmd implements SubCmd {

        // Auth0 configuration per environment
        private static final Map AUTH0_CONFIG = [
            'dev': [
                domain: 'seqera-development.eu.auth0.com',
                clientId: 'Ep2LhYiYmuV9hhz0dH6dbXVq0S7s7SWZ'
            ],
            'stage': [
                domain: 'seqera-stage.eu.auth0.com',
                clientId: '60cPDjI6YhoTPjyMTIBjGtxatSUwWswB'
            ],
            'prod': [
                domain: 'seqera.eu.auth0.com',
                clientId: 'FxCM8EJ76nNeHUDidSHkZfT8VtsrhHeL'
            ]
        ]

        String apiUrl

        private Map getAuth0Config(String environment) {
            def config = AUTH0_CONFIG[environment] as Map
            if (!config) {
                throw new RuntimeException("Unknown environment: ${environment}")
            }
            return config
        }


        @Override
        String getName() { 'login' }

        @Override
        void apply(List<String> args) {
            if (args.size() > 0) {
                throw new AbortOperationException("Too many arguments for login command")
            }

            // Check if TOWER_ACCESS_TOKEN environment variable is set
            def envToken = System.getenv('TOWER_ACCESS_TOKEN')
            if (envToken) {
                println ""
                ColorUtil.printColored("WARNING: Authentication token is already configured via TOWER_ACCESS_TOKEN environment variable.", "yellow bold")
                ColorUtil.printColored("${ColorUtil.colorize('nextflow auth login', 'cyan')} sets credentials using Nextflow config files, which take precedence over the environment variable.", "dim")
                ColorUtil.printColored(" however, caution is advised to avoid confusing behaviour.", "dim")
                println ""
            }

            // Check if tower.accessToken is already set
            def config = readConfig()
            def existingToken = config['tower.accessToken']
            if (existingToken) {
                ColorUtil.printColored("Error: Authentication token is already configured in Nextflow config.", "red")
                ColorUtil.printColored("Config file: ${ColorUtil.colorize(getConfigFile().toString(), 'magenta')}", "dim")
                println " Run ${ColorUtil.colorize('nextflow auth logout', 'cyan')} to remove the current authentication."
                return
            }

            println "Nextflow authentication with ${ColorUtil.colorize('Seqera Platform', 'cyan bold')}"
            ColorUtil.printColored(" - Authentication will be saved to: ${ColorUtil.colorize(getConfigFile().toString(), 'magenta')}", "dim")

            // Use provided URL or default
            if (!apiUrl) {
                apiUrl = 'https://api.cloud.seqera.io'
            } else if (!apiUrl.startsWith('http://') && !apiUrl.startsWith('https://')) {
                apiUrl = 'https://' + apiUrl
            }

            ColorUtil.printColored(" - Seqera Platform API endpoint: ${ColorUtil.colorize(apiUrl, 'magenta')} (can be customised with ${ColorUtil.colorize('-url', 'cyan')})", "dim")

            // Check if this is a cloud endpoint or enterprise
            def endpointInfo = getCloudEndpointInfo(apiUrl)
            if (endpointInfo.isCloud) {
                try {
                    performAuth0Login(apiUrl, endpointInfo.environment as String)
                } catch (Exception e) {
                    log.debug("Authentication failed", e)
                    println ""
                    throw new AbortOperationException("${e.message}")
                }
            } else {
                // Enterprise endpoint - use PAT authentication
                handleEnterpriseAuth(apiUrl)
            }
        }

        private void performAuth0Login(String apiUrl, String environment) {
            // Get Auth0 configuration for this environment
            def auth0Config = getAuth0Config(environment)

            // Start device authorization flow
            def deviceAuth = requestDeviceAuthorization(auth0Config)

            println ""
            ColorUtil.printColored("Confirmation code: ${ColorUtil.colorize(deviceAuth.user_code as String, 'yellow')}", "cyan bold")
            def urlWithCode = "${deviceAuth.verification_uri}?user_code=${deviceAuth.user_code}"
            println "${ColorUtil.colorize('Authentication URL:', 'cyan bold')} ${ColorUtil.colorize(urlWithCode, 'magenta')}"
            ColorUtil.printColored("\n[ Press Enter to open in browser ]", "cyan bold")
            System.in.read() // Wait for Enter key

            // Try to open browser automatically
            def browserOpened = false
            try {
                // Method 1: Java Desktop API
                if (java.awt.Desktop.isDesktopSupported()) {
                    def desktop = java.awt.Desktop.getDesktop()
                    if (desktop.isSupported(java.awt.Desktop.Action.BROWSE)) {
                        desktop.browse(new URI(urlWithCode))
                        browserOpened = true
                    }
                }

                // Method 2: Platform-specific commands
                if (!browserOpened) {
                    def os = System.getProperty("os.name").toLowerCase()
                    def command = []

                    if (os.contains("mac") || os.contains("darwin")) {
                        command = ["open", urlWithCode]
                    } else if (os.contains("win")) {
                        command = ["cmd", "/c", "start", urlWithCode]
                    } else {
                        // Linux and other Unix-like systems
                        def browsers = ["xdg-open", "firefox", "google-chrome", "chromium", "safari"]
                        for (browser in browsers) {
                            try {
                                new ProcessBuilder(browser, urlWithCode).start()
                                browserOpened = true
                                break
                            } catch (Exception ignored) {
                                // Try next browser
                            }
                        }
                    }

                    if (!browserOpened && command) {
                        new ProcessBuilder(command as String[]).start()
                        browserOpened = true
                    }
                }
            } catch (Exception ignored) {
                // Will handle below
            }

            if (!browserOpened) {
                ColorUtil.printColored("Could not open browser automatically. Please copy the URL above and open it manually in your browser.", "yellow")
            }
            print("${ColorUtil.colorize('Waiting for authentication...', 'dim', true)}")

            try {
                // Poll for device token
                def tokenData = pollForDeviceToken(deviceAuth.device_code as String, deviceAuth.interval as Integer ?: 5, auth0Config)
                def accessToken = tokenData['access_token'] as String

                // Verify login by calling /user-info
                def userInfo = callUserInfoApi(accessToken, apiUrl)
                ColorUtil.printColored("\n\nAuthentication successful!", "green")
                println "Logged in to ${ColorUtil.colorize(apiUrl.replace('api.', '').replace('/api', ''), 'magenta')} as: ${ColorUtil.colorize(userInfo.userName as String, 'cyan bold')}"

                // Generate PAT
                def pat = generatePAT(accessToken, apiUrl)

                // Save to config
                saveAuthToConfig(pat, apiUrl)
                println "Seqera Platform configuration saved to ${ColorUtil.colorize(getConfigFile().toString(), 'magenta')}"

            } catch (Exception e) {
                throw new RuntimeException("Authentication failed: ${e.message}", e)
            }
        }

        private Map requestDeviceAuthorization(Map auth0Config) {
            def deviceAuthUrl = "https://${auth0Config.domain}/oauth/device/code"

            def params = [
                'client_id': auth0Config.clientId,
                'scope': 'openid profile email offline_access',
                'audience': 'platform'
            ]

            def postData = params.collect { k, v ->
                "${URLEncoder.encode(k.toString(), 'UTF-8')}=${URLEncoder.encode(v.toString(), 'UTF-8')}"
            }.join('&')

            def connection = new URL(deviceAuthUrl).openConnection() as HttpURLConnection
            connection.requestMethod = 'POST'
            connection.setRequestProperty('Content-Type', 'application/x-www-form-urlencoded')
            connection.doOutput = true

            connection.outputStream.withWriter { writer ->
                writer.write(postData)
            }

            if (connection.responseCode != 200) {
                def error = connection.errorStream?.text ?: "HTTP ${connection.responseCode}"
                throw new RuntimeException("Device authorization request failed: ${error}")
            }

            def response = connection.inputStream.text
            def json = new groovy.json.JsonSlurper().parseText(response)
            return json as Map
        }

        private Map pollForDeviceToken(String deviceCode, int intervalSeconds, Map auth0Config) {
            def tokenUrl = "https://${auth0Config.domain}/oauth/token"
            def maxRetries = 60 // 5 minutes with 5-second intervals
            def retryCount = 0

            while (retryCount < maxRetries) {
                def params = [
                    'grant_type': 'urn:ietf:params:oauth:grant-type:device_code',
                    'device_code': deviceCode,
                    'client_id': auth0Config.clientId
                ]

                def postData = params.collect { k, v ->
                    "${URLEncoder.encode(k.toString(), 'UTF-8')}=${URLEncoder.encode(v.toString(), 'UTF-8')}"
                }.join('&')

                def connection = new URL(tokenUrl).openConnection() as HttpURLConnection
                connection.requestMethod = 'POST'
                connection.setRequestProperty('Content-Type', 'application/x-www-form-urlencoded')
                connection.doOutput = true

                connection.outputStream.withWriter { writer ->
                    writer.write(postData)
                }

                if (connection.responseCode == 200) {
                    def response = connection.inputStream.text
                    def json = new groovy.json.JsonSlurper().parseText(response)
                    return json as Map
                } else {
                    def errorResponse = connection.errorStream?.text
                    if (errorResponse) {
                        def errorJson = new groovy.json.JsonSlurper().parseText(errorResponse) as Map
                        def error = errorJson.error

                        if (error == 'authorization_pending') {
                            // User hasn't completed authorization yet, continue polling
                            print "${ColorUtil.colorize('.', 'dim', true)}"
                            System.out.flush()
                        } else if (error == 'slow_down') {
                            // Increase polling interval
                            intervalSeconds += 5
                            print "${ColorUtil.colorize('.', 'dim', true)}"
                            System.out.flush()
                        } else if (error == 'expired_token') {
                            throw new RuntimeException("The device code has expired. Please try again.")
                        } else if (error == 'access_denied') {
                            throw new RuntimeException("Access denied by user")
                        } else {
                            throw new RuntimeException("Token request failed: ${error} - ${errorJson.error_description ?: ''}")
                        }
                    } else {
                        throw new RuntimeException("Token request failed: HTTP ${connection.responseCode}")
                    }
                }

                // Wait before next poll
                Thread.sleep(intervalSeconds * 1000)
                retryCount++
            }

            throw new RuntimeException("Authentication timed out. Please try again.")
        }


        private void handleEnterpriseAuth(String apiUrl) {
            println ""
            println "Please generate a Personal Access Token from your ${ColorUtil.colorize('Seqera Platform', 'cyan bold')} instance."
            println "You can create one at: ${ColorUtil.colorize(apiUrl.replace('/api', '').replace('://api.', '') + '/tokens', 'magenta')}"
            println ""

            System.out.print("Enter your Personal Access Token: ")
            System.out.flush()

            def console = System.console()
            def pat = console ?
                new String(console.readPassword()) :
                new BufferedReader(new InputStreamReader(System.in)).readLine()

            if (!pat || pat.trim().isEmpty()) {
                throw new AbortOperationException("Personal Access Token is required for Seqera Platform Enterprise authentication")
            }

            // Save to config
            saveAuthToConfig(pat.trim(), apiUrl)
            ColorUtil.printColored("Personal Access Token saved to Nextflow config", "green")
            println "Config file: ${ColorUtil.colorize(getConfigFile().toString(), 'magenta')}"
        }

        private String generatePAT(String accessToken, String apiUrl) {
            def tokensUrl = "${apiUrl}/tokens"
            def username = System.getProperty("user.name")
            def timestamp = new Date().format("yyyy-MM-dd-HH-mm")
            def tokenName = "nextflow-${username}-${timestamp}"

            def requestBody = new groovy.json.JsonBuilder([name: tokenName]).toString()

            def connection = new URL(tokensUrl).openConnection() as HttpURLConnection
            connection.requestMethod = 'POST'
            connection.setRequestProperty('Authorization', "Bearer ${accessToken}")
            connection.setRequestProperty('Content-Type', 'application/json')
            connection.doOutput = true

            connection.outputStream.withWriter { writer ->
                writer.write(requestBody)
            }

            if (connection.responseCode != 200) {
                def error = connection.errorStream?.text ?: "HTTP ${connection.responseCode}"
                throw new RuntimeException("Failed to generate PAT: ${error}")
            }

            def response = connection.inputStream.text
            def json = new groovy.json.JsonSlurper().parseText(response) as Map
            return json.accessKey as String
        }

        private void saveAuthToConfig(String accessToken, String apiUrl) {
            def config = readConfig()
            config['tower.accessToken'] = accessToken
            config['tower.endpoint'] = apiUrl
            config['tower.enabled'] = true

            writeConfig(config)
        }

        @Override
        void usage(List<String> result) {
            result << 'Authenticate with Seqera Platform'
            result << "Usage: nextflow auth $name [-u <endpoint>]".toString()
            result << ''
            result << 'Options:'
            result << '  -u, -url <endpoint>    Seqera Platform API endpoint (default: https://api.cloud.seqera.io)'
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
                throw new AbortOperationException("Too many arguments for logout command")
            }

            // Check if TOWER_ACCESS_TOKEN environment variable is set
            def envToken = System.getenv('TOWER_ACCESS_TOKEN')
            if (envToken) {
                println ""
                ColorUtil.printColored("WARNING: TOWER_ACCESS_TOKEN environment variable is set.", "yellow bold")
                println " ${ColorUtil.colorize('nextflow auth logout', 'dim cyan')}${ColorUtil.colorize(' only removes credentials from Nextflow config files.', 'dim')}"
                ColorUtil.printColored(" The environment variable will remain unaffected.", "dim")
                println ""
            }

            // Check if tower.accessToken is set
            def config = readConfig()
            def existingToken = config['tower.accessToken']
            def endpoint = config['tower.endpoint'] ?: 'https://api.cloud.seqera.io'

            if (!existingToken) {
                ColorUtil.printColored("Error: No authentication token found in Nextflow config.", "red")
                println "Config file: ${ColorUtil.colorize(getConfigFile().toString(), 'magenta')}"
                return
            }

            ColorUtil.printColored(" - Found authentication token in config file: ${ColorUtil.colorize(getConfigFile().toString(), 'magenta')}", "dim")

            // Prompt user for API URL if not already configured
            def apiUrl = endpoint as String
            if (!apiUrl || apiUrl.isEmpty()) {
                apiUrl = promptForApiUrl()
            } else {
                ColorUtil.printColored(" - Using Seqera Platform endpoint: ${ColorUtil.colorize(apiUrl, 'magenta')}", "dim")
            }

            // Validate token by calling /user-info API
            try {
                def userInfo = callUserInfoApi(existingToken as String, apiUrl)
                ColorUtil.printColored(" - Token is valid for user: ${ColorUtil.colorize(userInfo.userName as String, 'cyan bold')}", "dim")

                // Only delete PAT from platform if this is a cloud endpoint
                if (isCloudEndpoint(apiUrl)) {
                    def tokenId = decodeTokenId(existingToken as String)
                    deleteTokenViaApi(existingToken as String, apiUrl, tokenId)
                } else {
                    println " - Enterprise installation detected - PAT will not be deleted from platform."
                }

                removeAuthFromConfig()

            } catch (Exception e) {
                println "Failed to validate or delete token: ${e.message}"
                println "Removing token from config anyway..."

                // Remove from config even if API calls fail
                removeAuthFromConfig()
                ColorUtil.printColored("Token removed from Nextflow config.", "green")
            }
        }

        private String decodeTokenId(String token) {
            try {
                // Decode base64 token
                def decoded = new String(Base64.decoder.decode(token), "UTF-8")

                // Parse JSON to extract token ID
                def json = new groovy.json.JsonSlurper().parseText(decoded) as Map
                def tokenId = json.tid

                if (!tokenId) {
                    throw new RuntimeException("No token ID found in decoded token")
                }

                return tokenId.toString()
            } catch (Exception e) {
                throw new RuntimeException("Failed to decode token ID: ${e.message}")
            }
        }

        private void deleteTokenViaApi(String token, String apiUrl, String tokenId) {
            def deleteUrl = "${apiUrl}/tokens/${tokenId}"
            def connection = new URL(deleteUrl).openConnection() as HttpURLConnection
            connection.requestMethod = 'DELETE'
            connection.setRequestProperty('Authorization', "Bearer ${token}")

            if (connection.responseCode != 200 && connection.responseCode != 204) {
                def error = connection.errorStream?.text ?: "HTTP ${connection.responseCode}"
                throw new RuntimeException("Failed to delete token: ${error}")
            }

            ColorUtil.printColored("\nToken successfully deleted from Seqera Platform.", "green")
        }

        private void removeAuthFromConfig() {
            def configFile = getConfigFile()

            if (Files.exists(configFile)) {
                def existingContent = Files.readString(configFile)
                def cleanedContent = cleanTowerConfig(existingContent)
                Files.writeString(configFile, cleanedContent, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)
            }

            ColorUtil.printColored("Authentication removed from Nextflow config.", "green")
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
                throw new AbortOperationException("Too many arguments for config command")
            }

            // Check if user is authenticated
            def config = readConfig()
            def existingToken = config['tower.accessToken'] ?: System.getenv('TOWER_ACCESS_TOKEN')
            def endpoint = config['tower.endpoint'] ?: 'https://api.cloud.seqera.io'

            if (!existingToken) {
                println "No authentication found. Please run ${ColorUtil.colorize('nextflow auth login', 'cyan')} first."
                return
            }

            println "Nextflow ${ColorUtil.colorize('Seqera Platform', 'cyan bold')} configuration"
            ColorUtil.printColored(" - Config file: ${ColorUtil.colorize(getConfigFile().toString(), 'magenta')}", "dim")

            // Check if token is from environment variable
            if (!config['tower.accessToken'] && System.getenv('TOWER_ACCESS_TOKEN')) {
                ColorUtil.printColored(" - Using access token from TOWER_ACCESS_TOKEN environment variable", "dim")
            }

            try {
                // Get user info to validate token and get user ID
                def userInfo = callUserInfoApi(existingToken as String, endpoint as String)
                ColorUtil.printColored(" - Authenticated as: ${ColorUtil.colorize(userInfo.userName as String, 'cyan bold')}", "dim")
                println ""

                // Track if any changes are made
                def configChanged = false

                // Configure tower.enabled
                configChanged |= configureEnabled(config)

                // Configure workspace
                configChanged |= configureWorkspace(config, existingToken as String, endpoint as String, userInfo.id as String)

                // Save updated config only if changes were made
                if (configChanged) {
                    writeConfig(config)
                    ColorUtil.printColored("\nConfiguration saved to ${ColorUtil.colorize(getConfigFile().toString(), 'magenta')}", "green")
                }

            } catch (Exception e) {
                throw new AbortOperationException("Failed to configure settings: ${e.message}")
            }
        }

        private boolean configureEnabled(Map config) {
            def currentEnabled = config.get('tower.enabled', false)

            println "Workflow monitoring settings. Current setting: ${currentEnabled ? ColorUtil.colorize('enabled', 'green') : ColorUtil.colorize('disabled', 'red')}"
            ColorUtil.printColored("  When enabled, all workflow runs are automatically monitored by Seqera Platform", "dim")
            ColorUtil.printColored("  When disabled, you can enable per-run with the ${ColorUtil.colorize('-with-tower', 'cyan')} flag", "dim")
            println ""

            System.out.print("${ColorUtil.colorize('Enable workflow monitoring for all runs?', 'cyan bold', true)} (${currentEnabled ? ColorUtil.colorize('Y', 'green') + '/n' : 'y/' + ColorUtil.colorize('N', 'red')}): ")
            System.out.flush()

            def reader = new BufferedReader(new InputStreamReader(System.in))
            def input = reader.readLine()?.trim()?.toLowerCase()

            if (input.isEmpty()) {
                // Keep current setting if user just presses enter
                return false
            } else if (input == 'y' || input == 'yes') {
                if (!currentEnabled) {
                    config['tower.enabled'] = true
                    return true
                } else {
                    return false
                }
            } else if (input == 'n' || input == 'no') {
                if (currentEnabled) {
                    config.remove('tower.enabled') // Don't set it to false, just remove it
                    return true
                } else {
                    return false
                }
            } else {
                println "Invalid input."
                return false
            }
        }

        private boolean configureWorkspace(Map config, String accessToken, String endpoint, String userId) {
            // Get all workspaces for the user
            def workspaces = getUserWorkspaces(accessToken, endpoint, userId)

            if (!workspaces) {
                println "\nNo workspaces found for your account."
                return false
            }

            // Show current workspace (check both config and env var)
            def currentWorkspaceId = config.get('tower.workspaceId')
            def envWorkspaceId = System.getenv('TOWER_WORKFLOW_ID')
            def effectiveWorkspaceId = currentWorkspaceId ?: envWorkspaceId
            def currentWorkspace = workspaces.find { ((Map)it).workspaceId.toString() == effectiveWorkspaceId?.toString() }


            // Group by organization
            def orgWorkspaces = workspaces.groupBy { ((Map)it).orgName ?: 'Personal' }

            // If 8 or fewer total options, show all at once
            if (workspaces.size() <= 8) {
                return selectWorkspaceFromAll(config, workspaces, currentWorkspaceId, envWorkspaceId)
            } else {
                // Two-stage selection: org first, then workspace
                return selectWorkspaceByOrg(config, orgWorkspaces, currentWorkspaceId, envWorkspaceId)
            }
        }

        private boolean selectWorkspaceFromAll(Map config, List workspaces, def currentWorkspaceId, def envWorkspaceId) {
            println "\nAvailable workspaces:"
            println "  0. ${ColorUtil.colorize('Personal workspace', 'cyan', true)} ${ColorUtil.colorize('[no organization]', 'dim', true)}"

            workspaces.eachWithIndex { workspace, index ->
                def ws = workspace as Map
                def prefix = ws.orgName ? "${ColorUtil.colorize(ws.orgName as String, 'cyan', true)} / " : ""
                println "  ${index + 1}. ${prefix}${ColorUtil.colorize(ws.workspaceName as String, 'magenta', true)} ${ColorUtil.colorize('[' + (ws.workspaceFullName as String) + ']', 'dim', true)}"
            }

            // Show current workspace and prepare prompt
            def currentWorkspace = workspaces.find { ((Map)it).workspaceId.toString() == (currentWorkspaceId ?: System.getenv('TOWER_WORKFLOW_ID'))?.toString() }
            def currentWorkspaceName
            if (currentWorkspace) {
                def workspace = currentWorkspace as Map
                def source = currentWorkspaceId ? "config" : "TOWER_WORKFLOW_ID env var"
                currentWorkspaceName = "${workspace.orgName} / ${workspace.workspaceName}"
            } else if (System.getenv('TOWER_WORKFLOW_ID')) {
                currentWorkspaceName = "TOWER_WORKFLOW_ID=${System.getenv('TOWER_WORKFLOW_ID')}"
            } else {
                currentWorkspaceName = "Personal workspace"
            }

            ColorUtil.printColored("\nSelect workspace (0-${workspaces.size()}, press Enter to keep as '${currentWorkspaceName}'): ", "bold cyan")
            System.out.flush()

            def reader = new BufferedReader(new InputStreamReader(System.in))
            def input = reader.readLine()?.trim()

            if (input.isEmpty()) {
                return false
            }

            try {
                def selection = Integer.parseInt(input)
                if (selection == 0) {
                    if (envWorkspaceId) {
                        return false
                    } else {
                        def hadWorkspaceId = config.containsKey('tower.workspaceId')
                        config.remove('tower.workspaceId')
                        config.remove('tower.workspaceId.comment')
                        return hadWorkspaceId
                    }
                } else if (selection > 0 && selection <= workspaces.size()) {
                    def selectedWorkspace = workspaces[selection - 1] as Map
                    def selectedId = selectedWorkspace.workspaceId.toString()
                    if (envWorkspaceId && selectedId == envWorkspaceId) {
                        return false
                    } else {
                        def currentId = config.get('tower.workspaceId')
                        config['tower.workspaceId'] = selectedId
                        config['tower.workspaceId.comment'] = "${selectedWorkspace.orgName} / ${selectedWorkspace.workspaceName} [${selectedWorkspace.workspaceFullName}]"
                        return currentId != selectedId
                    }
                } else {
                    println "Invalid selection."
                    return false
                }
            } catch (NumberFormatException e) {
                println "Invalid input."
                return false
            }
        }

        private boolean selectWorkspaceByOrg(Map config, Map orgWorkspaces, def currentWorkspaceId, def envWorkspaceId) {
            // Get current workspace info for prompts
            def allWorkspaces = []
            orgWorkspaces.values().each { workspaceList ->
                allWorkspaces.addAll(workspaceList as List)
            }
            def currentWorkspace = allWorkspaces.find { ((Map)it).workspaceId.toString() == (currentWorkspaceId ?: envWorkspaceId)?.toString() }
            def currentOrgName
            def currentWorkspaceName
            if (currentWorkspace) {
                def workspace = currentWorkspace as Map
                currentOrgName = workspace.orgName as String
                currentWorkspaceName = workspace.workspaceName as String
            } else if (envWorkspaceId) {
                currentOrgName = "TOWER_WORKFLOW_ID=${envWorkspaceId}"
                currentWorkspaceName = null
            } else {
                currentOrgName = "Personal"
                currentWorkspaceName = null
            }

            // First, select organization
            def orgs = orgWorkspaces.keySet().toList()

            println "\nAvailable organizations:"
            // Add Personal workspace option if not already in the list
            def hasPersonal = orgs.contains('Personal')
            if (!hasPersonal) {
                println "  1. ${ColorUtil.colorize('Personal', 'cyan', true)} ${ColorUtil.colorize('[Personal workspace - no organization]', 'dim', true)}"
                orgs.eachWithIndex { orgName, index ->
                    println "  ${index + 2}. ${ColorUtil.colorize(orgName as String, 'cyan', true)}"
                }
                System.out.print("${ColorUtil.colorize("Select organization (1-${orgs.size() + 1}, Enter to keep as '${currentOrgName}'): ", 'dim', true)}")
            } else {
                orgs.eachWithIndex { orgName, index ->
                    def displayName = orgName == 'Personal' ? 'Personal [Personal workspace - no organization]' : orgName
                    println "  ${index + 1}. ${ColorUtil.colorize(displayName as String, 'cyan', true)}"
                }
                System.out.print("${ColorUtil.colorize("Select organization (1-${orgs.size()}, Enter to keep as '${currentOrgName}'): ", 'dim', true)}")
            }
            System.out.flush()

            def reader = new BufferedReader(new InputStreamReader(System.in))
            def orgInput = reader.readLine()?.trim()

            if (orgInput.isEmpty()) {
                return false
            }

            try {
                def orgSelection = Integer.parseInt(orgInput)
                def maxOrgSelection = hasPersonal ? orgs.size() : orgs.size() + 1
                if (orgSelection < 1 || orgSelection > maxOrgSelection) {
                    println "Invalid selection."
                    return false
                }

                def selectedOrgName
                if (!hasPersonal && orgSelection == 1) {
                    // Personal workspace selected
                    if (envWorkspaceId) {
                        return false
                    } else {
                        def hadWorkspaceId = config.containsKey('tower.workspaceId')
                        config.remove('tower.workspaceId')
                        config.remove('tower.workspaceId.comment')
                        return hadWorkspaceId
                    }
                } else {
                    def orgIndex = hasPersonal ? orgSelection - 1 : orgSelection - 2
                    selectedOrgName = orgs[orgIndex]
                }

                def orgWorkspaceList = orgWorkspaces[selectedOrgName] as List

                println ""
                println "Select workspace in ${selectedOrgName}:"

                if (selectedOrgName == 'Personal') {
                    println "  0. Personal workspace (default)"
                }

                orgWorkspaceList.eachWithIndex { workspace, index ->
                    def ws = workspace as Map
                    println "  ${index + 1}. ${ColorUtil.colorize(ws.workspaceName as String, 'magenta', true)} ${ColorUtil.colorize('[' + (ws.workspaceFullName as String) + ']', 'dim', true)}"
                }

                def maxSelection = orgWorkspaceList.size()
                def keepAsText = currentWorkspaceName ? "'${currentWorkspaceName}'" : "'Personal workspace'"
                System.out.print("${ColorUtil.colorize("Select workspace (${selectedOrgName == 'Personal' ? '0-' : '1-'}${maxSelection}, Enter to keep as ${keepAsText}): ", 'dim', true)}")
                System.out.flush()

                def wsInput = reader.readLine()?.trim()
                if (wsInput.isEmpty()) {
                    return false
                }

                def wsSelection = Integer.parseInt(wsInput)
                if (selectedOrgName == 'Personal' && wsSelection == 0) {
                    if (envWorkspaceId) {
                        return false
                    } else {
                        def hadWorkspaceId = config.containsKey('tower.workspaceId')
                        config.remove('tower.workspaceId')
                        config.remove('tower.workspaceId.comment')
                        return hadWorkspaceId
                    }
                } else if (wsSelection > 0 && wsSelection <= orgWorkspaceList.size()) {
                    def selectedWorkspace = orgWorkspaceList[wsSelection - 1] as Map
                    def selectedId = selectedWorkspace.workspaceId.toString()
                    if (envWorkspaceId && selectedId == envWorkspaceId) {
                        return false
                    } else {
                        def currentId = config.get('tower.workspaceId')
                        config['tower.workspaceId'] = selectedId
                        config['tower.workspaceId.comment'] = "${selectedWorkspace.orgName} / ${selectedWorkspace.workspaceName} [${selectedWorkspace.workspaceFullName}]"
                        return currentId != selectedId
                    }
                } else {
                    println "Invalid selection."
                    return false
                }

            } catch (NumberFormatException e) {
                println "Invalid input."
                return false
            }
        }

        private List getUserWorkspaces(String accessToken, String endpoint, String userId) {
            def workspacesUrl = "${endpoint}/user/${userId}/workspaces"
            def connection = new URL(workspacesUrl).openConnection() as HttpURLConnection
            connection.requestMethod = 'GET'
            connection.setRequestProperty('Authorization', "Bearer ${accessToken}")

            if (connection.responseCode != 200) {
                def error = connection.errorStream?.text ?: "HTTP ${connection.responseCode}"
                throw new RuntimeException("Failed to get workspaces: ${error}")
            }

            def response = connection.inputStream.text
            def json = new groovy.json.JsonSlurper().parseText(response) as Map
            def orgsAndWorkspaces = json.orgsAndWorkspaces as List

            // Filter to only include actual workspaces (where workspaceId is not null)
            return orgsAndWorkspaces.findAll { ((Map)it).workspaceId != null }
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
                throw new AbortOperationException("Too many arguments for status command")
            }

            def config = readConfig()

            // Collect all status information
            List<List<String>> statusRows = []

            // API endpoint
            def endpointInfo = getConfigValue(config, 'tower.endpoint', 'TOWER_API_ENDPOINT', 'https://api.cloud.seqera.io')
            statusRows.add(['API endpoint', ColorUtil.colorize(endpointInfo.value as String, 'magenta'), endpointInfo.source as String])

            // API connection check
            def apiConnectionOk = checkApiConnection(endpointInfo.value as String)
            def connectionColor = apiConnectionOk ? 'green' : 'red'
            statusRows.add(['API connection', ColorUtil.colorize(apiConnectionOk ? 'OK' : 'ERROR', connectionColor), ''])

            // Authentication check
            def tokenInfo = getConfigValue(config, 'tower.accessToken', 'TOWER_ACCESS_TOKEN')
            if (tokenInfo.value) {
                try {
                    def userInfo = callUserInfoApi(tokenInfo.value as String, endpointInfo.value as String)
                    def currentUser = userInfo.userName as String
                    statusRows.add(['Authentication', "${ColorUtil.colorize('OK', 'green')} (user: ${ColorUtil.colorize(currentUser, 'cyan')})".toString(), tokenInfo.source as String])
                } catch (Exception e) {
                    statusRows.add(['Authentication', ColorUtil.colorize('ERROR', 'red'), 'failed'])
                }
            } else {
                statusRows.add(['Authentication', "${ColorUtil.colorize('ERROR', 'red')} ${ColorUtil.colorize('(no token)', 'dim')}".toString(), 'not set'])
            }

            // Monitoring enabled
            def enabledInfo = getConfigValue(config, 'tower.enabled', null, 'false')
            def enabledValue = enabledInfo.value?.toString()?.toLowerCase() in ['true', '1', 'yes'] ? 'Yes' : 'No'
            def enabledColor = enabledValue == 'Yes' ? 'green' : 'yellow'
            statusRows.add(['Workflow monitoring', ColorUtil.colorize(enabledValue, enabledColor), (enabledInfo.source ?: 'default') as String])

            // Default workspace
            def workspaceInfo = getConfigValue(config, 'tower.workspaceId', 'TOWER_WORKFLOW_ID')
            if (workspaceInfo.value) {
                // Try to get workspace name from API if we have a token
                def workspaceDetails = null
                if (tokenInfo.value) {
                    workspaceDetails = getWorkspaceDetailsFromApi(tokenInfo.value as String, endpointInfo.value as String, workspaceInfo.value as String)
                }

                if (workspaceDetails) {
                    // Add workspace ID row
                    statusRows.add(['Default workspace ID', ColorUtil.colorize(workspaceInfo.value as String, 'blue'), workspaceInfo.source as String])
                    // Add org/name row
                    statusRows.add([' - workspace name', "${ColorUtil.colorize(workspaceDetails.orgName as String, 'cyan bold')} / ${ColorUtil.colorize(workspaceDetails.workspaceName as String, 'cyan')}".toString(), ''])
                    // Add full name row (truncate if too long)
                    def fullName = workspaceDetails.workspaceFullName as String
                    def truncatedFullName = fullName.length() > 50 ? fullName.substring(0, 47) + '...' : fullName
                    statusRows.add([' - workspace full name', ColorUtil.colorize(truncatedFullName, 'cyan dim'), ''])
                } else {
                    statusRows.add(['Default workspace', ColorUtil.colorize(workspaceInfo.value as String, 'blue', true), workspaceInfo.source as String])
                }
            } else {
                statusRows.add(['Default workspace', ColorUtil.colorize('Personal workspace', 'cyan', true), 'default'])
            }

            // Print table
            println ""
            printStatusTable(statusRows)
        }

        private void printStatusTable(List<List<String>> rows) {
            if (!rows) return

            // Calculate column widths (accounting for ANSI codes)
            def col1Width = rows.collect { stripAnsiCodes(it[0]).length() }.max()
            def col2Width = rows.collect { stripAnsiCodes(it[1]).length() }.max()
            def col3Width = rows.collect { stripAnsiCodes(it[2]).length() }.max()

            // Add some padding
            col1Width = Math.max(col1Width, 15) + 2
            col2Width = Math.max(col2Width, 15) + 2
            col3Width = Math.max(col3Width, 10) + 2

            // Print table header
            ColorUtil.printColored("${'Setting'.padRight(col1Width)} ${'Value'.padRight(col2Width)} Source", "cyan bold")
            println "${'-' * col1Width} ${'-' * col2Width} ${'-' * col3Width}"

            // Print rows
            rows.each { row ->
                def paddedCol1 = padStringWithAnsi(row[0], col1Width)
                def paddedCol2 = padStringWithAnsi(row[1], col2Width)
                def paddedCol3 = ColorUtil.colorize(row[2], 'dim', true)
                println "${paddedCol1} ${paddedCol2} ${paddedCol3}"
            }
        }

        private String stripAnsiCodes(String text) {
            return text?.replaceAll(/\u001B\[[0-9;]*m/, '') ?: ''
        }

        private String padStringWithAnsi(String text, int width) {
            def plainText = stripAnsiCodes(text)
            def padding = width - plainText.length()
            return padding > 0 ? text + (' ' * padding) : text
        }

        private String shortenPath(String path) {
            def userHome = System.getProperty('user.home')
            if (path.startsWith(userHome)) {
                return '~' + path.substring(userHome.length())
            }
            return path
        }

        private Map getConfigValue(Map config, String configKey, String envVarName, String defaultValue = null) {
            def configValue = config[configKey]
            def envValue = envVarName ? System.getenv(envVarName) : null
            def effectiveValue = configValue ?: envValue ?: defaultValue

            def source = null
            if (configValue) {
                source = shortenPath(getConfigFile().toString())
            } else if (envValue) {
                source = "env var \$${envVarName}"
            } else if (defaultValue) {
                source = "default"
            }

            return [
                value: effectiveValue,
                source: source,
                fromConfig: configValue != null,
                fromEnv: envValue != null,
                isDefault: !configValue && !envValue
            ]
        }

        private String getWorkspaceNameFromApi(String accessToken, String endpoint, String workspaceId) {
            try {
                // Get user info to get user ID
                def userInfo = callUserInfoApi(accessToken, endpoint)
                def userId = userInfo.id as String

                // Get workspaces for the user
                def workspacesUrl = "${endpoint}/user/${userId}/workspaces"
                def connection = new URL(workspacesUrl).openConnection() as HttpURLConnection
                connection.requestMethod = 'GET'
                connection.connectTimeout = 10000 // 10 second timeout
                connection.readTimeout = 10000
                connection.setRequestProperty('Authorization', "Bearer ${accessToken}")

                if (connection.responseCode != 200) {
                    return null
                }

                def response = connection.inputStream.text
                def json = new groovy.json.JsonSlurper().parseText(response) as Map
                def orgsAndWorkspaces = json.orgsAndWorkspaces as List

                // Find the workspace with matching ID
                def workspace = orgsAndWorkspaces.find { ((Map)it).workspaceId?.toString() == workspaceId }
                if (workspace) {
                    def ws = workspace as Map
                    return "${ws.orgName} / ${ws.workspaceName} [${ws.workspaceFullName}]"
                }

                return null
            } catch (Exception e) {
                return null
            }
        }

        private Map getWorkspaceDetailsFromApi(String accessToken, String endpoint, String workspaceId) {
            try {
                // Get user info to get user ID
                def userInfo = callUserInfoApi(accessToken, endpoint)
                def userId = userInfo.id as String

                // Get workspaces for the user
                def workspacesUrl = "${endpoint}/user/${userId}/workspaces"
                def connection = new URL(workspacesUrl).openConnection() as HttpURLConnection
                connection.requestMethod = 'GET'
                connection.connectTimeout = 10000 // 10 second timeout
                connection.readTimeout = 10000
                connection.setRequestProperty('Authorization', "Bearer ${accessToken}")

                if (connection.responseCode != 200) {
                    return null
                }

                def response = connection.inputStream.text
                def json = new groovy.json.JsonSlurper().parseText(response) as Map
                def orgsAndWorkspaces = json.orgsAndWorkspaces as List

                // Find the workspace with matching ID
                def workspace = orgsAndWorkspaces.find { ((Map)it).workspaceId?.toString() == workspaceId }
                if (workspace) {
                    def ws = workspace as Map
                    return [
                        orgName: ws.orgName,
                        workspaceName: ws.workspaceName,
                        workspaceFullName: ws.workspaceFullName
                    ]
                }

                return null
            } catch (Exception e) {
                return null
            }
        }

        private boolean checkApiConnection(String endpoint) {
            try {
                def serviceInfoUrl = "${endpoint}/service-info"
                def connection = new URL(serviceInfoUrl).openConnection() as HttpURLConnection
                connection.requestMethod = 'GET'
                connection.connectTimeout = 10000 // 10 second timeout
                connection.readTimeout = 10000

                return connection.responseCode == 200
            } catch (Exception e) {
                return false
            }
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
