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

    CmdAuth() {
        commands.add(new LoginCmd())
        commands.add(new LogoutCmd())
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
            getCmd(args).apply(args.drop(1))
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
            if (value instanceof String) {
                configText.append("    ${configKey} = '${value}'\n")
            } else {
                configText.append("    ${configKey} = ${value}\n")
            }
        }
        configText.append("}\n")

        Files.writeString(configFile, configText.toString(), StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)
    }

    //
    // nextflow auth login
    //
    class LoginCmd implements SubCmd {

        private static final String AUTH0_DOMAIN = "seqera-development.eu.auth0.com"
        private static final String AUTH0_CLIENT_ID = "Ep2LhYiYmuV9hhz0dH6dbXVq0S7s7SWZ"
        private static final int CALLBACK_PORT = 8085

        @Override
        String getName() { 'login' }

        @Override
        void apply(List<String> args) {
            if (args.size() > 0) {
                throw new AbortOperationException("Too many arguments for login command")
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
            println ""

            // Prompt user for API URL
            def apiUrl = promptForApiUrl()

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
            font-family: Inter, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            margin: 0;
            padding: 0;
            display: flex;
            justify-content: center;
            align-items: center;
            min-height: 100vh;
        }
        .container {
            text-align: center;
            padding: 2rem;
            background: white;
            border-radius: 8px;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
        }
        h1 { color: #1e293b; margin: 0 0 0.5rem; }
        p { color: #64748b; margin: 0; }
    </style>
</head>
<body>
    <div class="container">
        <h1>Nextflow authentication successful!</h1>
        <p>You can now close this window.</p>
    </div>
</body>
</html>""")
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

            // Ask user if they want to enable workflow monitoring by default
            System.out.print("Enable workflow monitoring for all runs? (Y/n): ")
            System.out.flush()

            def reader = new BufferedReader(new InputStreamReader(System.in))
            def input = reader.readLine()?.trim()?.toLowerCase()

            if (input == 'n' || input == 'no') {
                println "Workflow monitoring not enabled by default. You can enable it per-run with -with-tower."
            } else {
                config['tower.enabled'] = true
            }

            writeConfig(config)
        }

        @Override
        void usage(List<String> result) {
            result << 'Authenticate with Seqera Platform'
            result << "Usage: nextflow auth $name".toString()
            result << ''
            result << 'This command will:'
            result << '  1. Prompt for Seqera Platform API endpoint'
            result << '  2. Open browser for OAuth2 authentication (Cloud) or prompt for PAT (Enterprise)'
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

            println "Nextflow authentication logout"
            println ""

            // Check if tower.accessToken is set
            def config = readConfig()
            def existingToken = config['tower.accessToken']
            def endpoint = config['tower.endpoint'] ?: 'https://api.cloud.seqera.io'

            if (!existingToken) {
                println "No authentication token found in Nextflow config. Already logged out."
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
}
