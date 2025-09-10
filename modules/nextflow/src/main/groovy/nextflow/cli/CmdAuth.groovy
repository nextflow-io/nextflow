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
import nextflow.exception.AbortOperationException

import java.awt.Desktop
import java.net.ServerSocket
import java.net.URI
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardOpenOption
import java.security.MessageDigest
import java.security.SecureRandom
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

    // Shared methods
    private String promptForApiUrl() {
        System.out.print("Seqera Platform API endpoint [Default https://api.cloud.seqera.io]: ")
        System.out.flush()

        def reader = new BufferedReader(new InputStreamReader(System.in))
        def input = reader.readLine()?.trim()

        return input?.isEmpty() || input == null ? 'https://api.cloud.seqera.io' : input
    }

    private boolean isCloudEndpoint(String apiUrl) {
        return apiUrl == 'https://api.cloud.seqera.io' ||
               apiUrl == 'https://api.cloud.stage-seqera.io' ||
               apiUrl == 'https://api.cloud.dev-seqera.io' ||
               apiUrl == 'https://cloud.seqera.io/api' ||
               apiUrl == 'https://cloud.stage-seqera.io/api' ||
               apiUrl == 'https://cloud.dev-seqera.io/api'
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

        private static final String AUTH0_DOMAIN = "seqera-development.eu.auth0.com"
        private static final String AUTH0_CLIENT_ID = "Ep2LhYiYmuV9hhz0dH6dbXVq0S7s7SWZ"
        private static final int CALLBACK_PORT = 8085

        String apiUrl

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
                println "WARNING: Authentication token is already configured via TOWER_ACCESS_TOKEN environment variable."
                println "nextflow auth login' sets credentials using Nextflow config files, which take precedence over the environment variable."
                println "however, caution is advised to avoid confusing behaviour."
                println ""
            }

            // Check if tower.accessToken is already set
            def config = readConfig()
            def existingToken = config['tower.accessToken']
            if (existingToken) {
                println "Authentication token is already configured in Nextflow config."
                println "Config file: ${getConfigFile()}"
                println "Run 'nextflow auth logout' to remove the current authentication."
                return
            }

            println "Nextflow authentication with Seqera Platform"
            println " - Authentication will be saved to: ${getConfigFile()}"

            // Use provided URL or default
            if (!apiUrl) {
                apiUrl = 'https://api.cloud.seqera.io'
            } else if (!apiUrl.startsWith('http://') && !apiUrl.startsWith('https://')) {
                apiUrl = 'https://' + apiUrl
            }

            println " - Seqera Platform API endpoint: ${apiUrl}"
            println "   (can be customised with -u option)"
            Thread.sleep(2000)

            // Check if this is a cloud endpoint or enterprise
            if (isCloudEndpoint(apiUrl)) {
                try {
                    performAuth0Login(apiUrl)
                } catch (Exception e) {
                    log.debug("Authentication failed", e)
                    throw new AbortOperationException("Authentication failed: ${e.message}")
                }
            } else {
                // Enterprise endpoint - use PAT authentication
                handleEnterpriseAuth(apiUrl)
            }
        }

        private void performAuth0Login(String apiUrl) {
            println " - Opening browser for authentication..."

            // Generate PKCE parameters
            def codeVerifier = generateCodeVerifier()
            def codeChallenge = generateCodeChallenge(codeVerifier)
            def state = generateState()

            // Start local server for callback
            def callbackFuture = startCallbackServer(state)

            // Build authorization URL
            def authUrl = buildAuthUrl(codeChallenge, state)

            // Open browser
            openBrowser(authUrl)

            try {
                // Wait for callback with timeout
                def authCode = callbackFuture.get(5, TimeUnit.MINUTES)

                // Exchange code for token
                def tokenData = exchangeCodeForToken(authCode, codeVerifier)
                def accessToken = tokenData['access_token'] as String

                // Verify login by calling /user-info
                def userInfo = callUserInfoApi(accessToken, apiUrl)
                println " - Authentication successful! Logged in as: ${userInfo.userName}"

                // Generate PAT
                def pat = generatePAT(accessToken, apiUrl)

                // Save to config
                saveAuthToConfig(pat, apiUrl)
                println " - Seqera Platform configuration saved to ${getConfigFile()}"

            } catch (Exception e) {
                throw new RuntimeException("Authentication timeout or failed: ${e.message}", e)
            }
        }

        private String generateCodeVerifier() {
            def random = new SecureRandom()
            def bytes = new byte[32]
            random.nextBytes(bytes)
            return Base64.urlEncoder.withoutPadding().encodeToString(bytes)
        }

        private String generateCodeChallenge(String codeVerifier) {
            def digest = MessageDigest.getInstance("SHA-256")
            def hash = digest.digest(codeVerifier.getBytes("UTF-8"))
            return Base64.urlEncoder.withoutPadding().encodeToString(hash)
        }

        private String generateState() {
            def random = new SecureRandom()
            def bytes = new byte[16]
            random.nextBytes(bytes)
            return Base64.urlEncoder.withoutPadding().encodeToString(bytes)
        }

        private CompletableFuture<String> startCallbackServer(String expectedState) {
            def future = new CompletableFuture<String>()

            Thread.start {
                try {
                    def server = new ServerSocket(CALLBACK_PORT)
                    server.soTimeout = 300000 // 5 minutes

                    def socket = server.accept()
                    def input = new BufferedReader(new InputStreamReader(socket.inputStream))
                    def output = new PrintWriter(socket.outputStream)

                    def requestLine = input.readLine()
                    if (requestLine?.startsWith("GET /callback")) {
                        def query = requestLine.split("\\?")[1]?.split(" ")[0]
                        def params = parseQueryParams(query)

                        if (params['state'] != expectedState) {
                            throw new RuntimeException("Invalid state parameter")
                        }

                        if (params['error']) {
                            throw new RuntimeException("Auth error: ${params['error']} - ${params['error_description']}")
                        }

                        def code = params['code']
                        if (!code) {
                            throw new RuntimeException("No authorization code received")
                        }

                        // Send success response
                        output.println("HTTP/1.1 200 OK")
                        output.println("Content-Type: text/html")
                        output.println()
                        output.println("""<html>
<head>
    <style>
        body {
            font-family: Inter, sans-serif;
            margin: 0;
            padding: 0;
            display: flex;
            justify-content: center;
            align-items: center;
            min-height: 100vh;
        }
        .container {
            padding: 2rem;
            background: white;
            border: 1px solid rgba(69, 63, 81, .2);
        }
        h1 {
            color: #1e293b;
            margin: 0 0 1rem;
            font-size: 21px;
            font-weight: 600;
        }
        p { color: #64748b; margin: 0; }
    </style>
</head>
<body>
    <div class="container">
        <h1>Nextflow authentication successful!</h1>
        <p>You can now close this window.</p>
    </div>
</body>
</html>
""")
                        output.flush()

                        future.complete(code)
                    } else {
                        future.completeExceptionally(new RuntimeException("Invalid callback request"))
                    }

                    socket.close()
                    server.close()

                } catch (Exception e) {
                    future.completeExceptionally(e)
                }
            }

            return future
        }

        private Map<String, String> parseQueryParams(String query) {
            Map<String, String> params = [:]
            if (query) {
                query.split('&').each { param ->
                    def parts = param.split('=', 2)
                    if (parts.length == 2) {
                        params[URLDecoder.decode(parts[0], "UTF-8")] = URLDecoder.decode(parts[1], "UTF-8")
                    }
                }
            }
            return params
        }

        private String buildAuthUrl(String codeChallenge, String state) {
            def params = [
                'response_type': 'code',
                'client_id': AUTH0_CLIENT_ID,
                'redirect_uri': "http://localhost:${CALLBACK_PORT}/callback",
                'scope': 'openid profile email offline_access',
                'audience': 'platform',
                'code_challenge': codeChallenge,
                'code_challenge_method': 'S256',
                'state': state,
                'prompt': 'login'
            ]

            def query = params.collect { k, v ->
                "${URLEncoder.encode(k.toString(), 'UTF-8')}=${URLEncoder.encode(v.toString(), 'UTF-8')}"
            }.join('&')

            return "https://${AUTH0_DOMAIN}/authorize?${query}"
        }

        private void openBrowser(String url) {
            try {
                if (Desktop.isDesktopSupported()) {
                    def desktop = Desktop.desktop
                    if (desktop.isSupported(Desktop.Action.BROWSE)) {
                        desktop.browse(new URI(url))
                        return
                    }
                }

                // Fallback: try OS-specific commands
                def os = System.getProperty("os.name").toLowerCase()
                if (os.contains("mac")) {
                    Runtime.runtime.exec(["open", url] as String[])
                    return
                } else if (os.contains("linux")) {
                    Runtime.runtime.exec(["xdg-open", url] as String[])
                    return
                } else if (os.contains("windows")) {
                    Runtime.runtime.exec(["rundll32", "url.dll,FileProtocolHandler", url] as String[])
                    return
                }

                // If all else fails, show URL
                println "Could not open browser automatically. Please visit: ${url}"

            } catch (Exception e) {
                log.debug("Failed to open browser", e)
                println "Failed to open browser automatically. Please visit: ${url}"
            }
        }

        private Map exchangeCodeForToken(String authCode, String codeVerifier) {
            def tokenUrl = "https://${AUTH0_DOMAIN}/oauth/token"

            def params = [
                'grant_type': 'authorization_code',
                'client_id': AUTH0_CLIENT_ID,
                'code': authCode,
                'redirect_uri': "http://localhost:${CALLBACK_PORT}/callback",
                'code_verifier': codeVerifier
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

            if (connection.responseCode != 200) {
                def error = connection.errorStream?.text ?: "HTTP ${connection.responseCode}"
                throw new RuntimeException("Token exchange failed: ${error}")
            }

            def response = connection.inputStream.text
            def json = new groovy.json.JsonSlurper().parseText(response)
            return json as Map
        }

        private void handleEnterpriseAuth(String apiUrl) {
            println ""
            println "Please generate a Personal Access Token from your Seqera Platform instance."
            println "You can create one at: ${apiUrl.replace('/api', '').replace('://api.', '')}/tokens"
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
            println "Personal Access Token saved to Nextflow config"
            println "Config file: ${getConfigFile()}"
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
            result << '  1. Open browser for OAuth2 authentication (Cloud) or prompt for PAT (Enterprise)'
            result << '  2. Generate and save access token to home-directory Nextflow config'
            result << '  3. Configure tower.accessToken, tower.endpoint, and tower.enabled settings'
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

            println "Nextflow authentication logout"
            println ""

            // Check if TOWER_ACCESS_TOKEN environment variable is set
            def envToken = System.getenv('TOWER_ACCESS_TOKEN')
            if (envToken) {
                println "WARNING: TOWER_ACCESS_TOKEN environment variable is set."
                println "'nextflow auth logout' only removes credentials from Nextflow config files."
                println "The environment variable will remain unaffected."
                println ""
            }

            // Check if tower.accessToken is set
            def config = readConfig()
            def existingToken = config['tower.accessToken']
            def endpoint = config['tower.endpoint'] ?: 'https://api.cloud.seqera.io'

            if (!existingToken) {
                println "No authentication token found in Nextflow config."
                println "Config file: ${getConfigFile()}"
                return
            }

            println " - Found authentication token in config file: ${getConfigFile()}"

            // Prompt user for API URL if not already configured
            def apiUrl = endpoint as String
            if (!apiUrl || apiUrl.isEmpty()) {
                apiUrl = promptForApiUrl()
            } else {
                println " - Using Seqera Platform endpoint: ${apiUrl}"
            }

            // Validate token by calling /user-info API
            try {
                def userInfo = callUserInfoApi(existingToken as String, apiUrl)
                println " - Token is valid for user: ${userInfo.userName}"

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
                println "Token removed from Nextflow config."
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

            println " - Token successfully deleted from Seqera Platform."
        }

        private void removeAuthFromConfig() {
            def configFile = getConfigFile()

            if (Files.exists(configFile)) {
                def existingContent = Files.readString(configFile)
                def cleanedContent = cleanTowerConfig(existingContent)
                Files.writeString(configFile, cleanedContent, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)
            }

            println " - Authentication removed from Nextflow config."
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
                println "No authentication found. Please run 'nextflow auth login' first."
                return
            }

            println "Nextflow Seqera Platform configuration"
            println " - Config file: ${getConfigFile()}"

            try {
                // Get user info to validate token and get user ID
                def userInfo = callUserInfoApi(existingToken as String, endpoint as String)
                println " - Authenticated as: ${userInfo.userName}"
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
                    println " - Configuration saved to ${getConfigFile()}"
                }

            } catch (Exception e) {
                throw new AbortOperationException("Failed to configure settings: ${e.message}")
            }
        }

        private boolean configureEnabled(Map config) {
            def currentEnabled = config.get('tower.enabled', false)

            println "Workflow monitoring settings:"
            println "  Current: ${currentEnabled ? 'enabled' : 'disabled'}"
            println "  When enabled, all workflow runs are automatically monitored by Seqera Platform"
            println "  When disabled, you can enable per-run with the -with-tower flag"
            println ""

            System.out.print("Enable workflow monitoring for all runs? (${currentEnabled ? 'Y/n' : 'y/N'}): ")
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

            println "Default workspace settings:"
            if (currentWorkspace) {
                def workspace = currentWorkspace as Map
                def source = currentWorkspaceId ? "config" : (envWorkspaceId ? "TOWER_WORKFLOW_ID env var" : "config")
                println "  Current: ${workspace.orgName} / ${workspace.workspaceFullName} (from ${source})"
            } else if (envWorkspaceId) {
                println "  Current: TOWER_WORKFLOW_ID=${envWorkspaceId} (workspace not found in available workspaces)"
            } else {
                println "  Current: Personal workspace (default)"
            }
            println ""

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
            println "Select default workspace:"
            println "  0. Personal workspace (no organization)"

            workspaces.eachWithIndex { workspace, index ->
                def ws = workspace as Map
                def prefix = ws.orgName ? "${ws.orgName} / " : ""
                println "  ${index + 1}. ${prefix}${ws.workspaceFullName}"
            }

            System.out.print("Select workspace (0-${workspaces.size()}, Enter to keep current): ")
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
                        config['tower.workspaceId.comment'] = "${selectedWorkspace.orgName} / ${selectedWorkspace.workspaceFullName}"
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
            // First, select organization
            def orgs = orgWorkspaces.keySet().toList()

            println "Select organization:"
            orgs.eachWithIndex { orgName, index ->
                println "  ${index + 1}. ${orgName}"
            }

            System.out.print("Select organization (1-${orgs.size()}, Enter to keep current): ")
            System.out.flush()

            def reader = new BufferedReader(new InputStreamReader(System.in))
            def orgInput = reader.readLine()?.trim()

            if (orgInput.isEmpty()) {
                return false
            }

            try {
                def orgSelection = Integer.parseInt(orgInput)
                if (orgSelection < 1 || orgSelection > orgs.size()) {
                    println "Invalid selection."
                    return false
                }

                def selectedOrgName = orgs[orgSelection - 1]
                def orgWorkspaceList = orgWorkspaces[selectedOrgName] as List

                println ""
                println "Select workspace in ${selectedOrgName}:"

                if (selectedOrgName == 'Personal') {
                    println "  0. Personal workspace (default)"
                }

                orgWorkspaceList.eachWithIndex { workspace, index ->
                    def ws = workspace as Map
                    println "  ${index + 1}. ${ws.workspaceFullName}"
                }

                def maxSelection = orgWorkspaceList.size()
                System.out.print("Select workspace (${selectedOrgName == 'Personal' ? '0-' : '1-'}${maxSelection}, Enter to keep current): ")
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
                        config['tower.workspaceId.comment'] = "${selectedWorkspace.orgName} / ${selectedWorkspace.workspaceFullName}"
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

            // API endpoint
            def endpointInfo = getConfigValue(config, 'tower.endpoint', 'TOWER_API_ENDPOINT', 'https://api.cloud.seqera.io')
            println "API endpoint: ${endpointInfo.value} (${endpointInfo.source})"

            // API connection check
            def apiConnectionOk = checkApiConnection(endpointInfo.value as String)
            println "API connection check: ${apiConnectionOk ? 'OK' : 'ERROR'}"

            // Access token status
            def tokenInfo = getConfigValue(config, 'tower.accessToken', 'TOWER_ACCESS_TOKEN')
            if (tokenInfo.value) {
                println "Access token: Configured (${tokenInfo.source})"
            } else {
                println "Access token: Not configured"
            }

            // Authentication check
            if (tokenInfo.value) {
                try {
                    def userInfo = callUserInfoApi(tokenInfo.value as String, endpointInfo.value as String)
                    def currentUser = userInfo.userName
                    println "Authentication: Success (${currentUser})"
                } catch (Exception e) {
                    println "Authentication: Error (${e})"
                }
            } else {
                println "Authentication: ERROR (no token)"
            }

            // Monitoring enabled
            def enabledInfo = getConfigValue(config, 'tower.enabled', null, 'false')
            def enabledValue = enabledInfo.value?.toString()?.toLowerCase() in ['true', '1', 'yes'] ? 'Yes' : 'No'
            println "Local workflow monitoring enabled: ${enabledValue} (${enabledInfo.source ?: 'default'})"

            // Default workspace
            def workspaceInfo = getConfigValue(config, 'tower.workspaceId', 'TOWER_WORKFLOW_ID')
            if (workspaceInfo.value) {
                // Try to get workspace name from API if we have a token
                def workspaceName = null
                if (tokenInfo.value) {
                    workspaceName = getWorkspaceNameFromApi(tokenInfo.value as String, endpointInfo.value as String, workspaceInfo.value as String)
                }

                if (workspaceName) {
                    println "Default workspace: '${workspaceName}' [${workspaceInfo.value}] (${workspaceInfo.source})"
                } else {
                    println "Default workspace: ${workspaceInfo.value} (${workspaceInfo.source})"
                }
            } else {
                println "Default workspace: Personal workspace (default)"
            }
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
                    return "${ws.orgName} / ${ws.workspaceFullName}"
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
