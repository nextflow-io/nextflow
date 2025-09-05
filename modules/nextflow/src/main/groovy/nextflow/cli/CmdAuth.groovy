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
import nextflow.exception.AbortOperationException

import java.awt.Desktop
import java.net.ServerSocket
import java.net.URI
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
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

    // Shared helper methods for both login and logout
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

    private String findTokenInShellConfig(String token) {
        def shellConfigInfo = detectShellConfig()

        if (shellConfigInfo.configFile) {
            def content = Files.readString(Paths.get(shellConfigInfo.configFile as String))
            if (content.contains("export TOWER_ACCESS_TOKEN=${token}") || content.contains("export TOWER_ACCESS_TOKEN=\"${token}\"")) {
                return shellConfigInfo.configFile
            }
        }

        return null
    }

    private Map detectShellConfig() {
        def shell = System.getenv("SHELL")
        def homeDir = System.getProperty("user.home")
        def configFile
        def shellName

        if (shell?.contains("zsh")) {
            shellName = "zsh"
            configFile = "${homeDir}/.zshrc"
        } else if (shell?.contains("fish")) {
            shellName = "fish"
            configFile = "${homeDir}/.config/fish/config.fish"
        } else if (shell?.contains("bash")) {
            shellName = "bash"
            // Check for .bash_profile first, then .bashrc
            def bashProfile = "${homeDir}/.bash_profile"
            def bashrc = "${homeDir}/.bashrc"
            configFile = Files.exists(Paths.get(bashProfile)) ? bashProfile : bashrc
        } else {
            // Unrecognized shell
            shellName = "unknown"
            configFile = null
        }

        return [shell: shellName, configFile: configFile]
    }

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

            // Check if TOWER_ACCESS_TOKEN is already set
            def existingToken = System.getenv("TOWER_ACCESS_TOKEN")
            if (existingToken) {
                println "TOWER_ACCESS_TOKEN environment variable is already set."

                // Try to find the token in shell config files
                def configFile = findTokenInShellConfig(existingToken)
                if (configFile) {
                    println "Token found in: ${configFile}"
                }

                println "Run `nextflow auth status` to view details."
                return
            }

            println "Nextflow authentication with Seqera Platform"
            println ""

            // Detect shell and config file
            def shellConfigInfo = detectShellConfig()

            if (shellConfigInfo.shell == "unknown") {
                println "Unrecognized shell detected. After authentication, you will need to manually add the TOWER_ACCESS_TOKEN export to your shell configuration file."
                println ""
            } else {
                // Check if config file exists before proceeding
                def configFile = Paths.get(shellConfigInfo.configFile as String)
                if (!Files.exists(configFile)) {
                    throw new AbortOperationException("Shell ${shellConfigInfo.shell} was detected but config file ${shellConfigInfo.configFile} was not found. Please create this file first.")
                }

                println "A Personal Access Token will be generated and saved to: ${shellConfigInfo.configFile}"
                println ""
            }

            // Prompt user for API URL
            def apiUrl = promptForApiUrl()

            // Check if this is a cloud endpoint or enterprise
            if (isCloudEndpoint(apiUrl)) {
                try {
                    performAuth0Login(apiUrl, shellConfigInfo)
                } catch (Exception e) {
                    log.debug("Authentication failed", e)
                    throw new AbortOperationException("Authentication failed: ${e.message}")
                }
            } else {
                // Enterprise endpoint - use PAT authentication
                handleEnterpriseAuth(apiUrl)
            }
        }

        private String promptForApiUrl() {
            System.out.print("Seqera Platform API URL [Default Seqera Cloud, https://api.cloud.seqera.io]: ")
            System.out.flush()

            def reader = new BufferedReader(new InputStreamReader(System.in))
            def input = reader.readLine()?.trim()

            return input?.isEmpty() || input == null ? 'https://api.cloud.seqera.io' : input
        }

        private void performAuth0Login(String apiUrl, Map shellConfigInfo) {
            println "- Opening browser for authentication..."

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

            println "- Waiting for authentication to complete..."

            try {
                // Wait for callback with timeout
                def authCode = callbackFuture.get(5, TimeUnit.MINUTES)

                // Exchange code for token
                def tokenData = exchangeCodeForToken(authCode, codeVerifier)
                def accessToken = tokenData['access_token'] as String

                // Verify login by calling /user-info
                def userInfo = callUserInfoApi(accessToken, apiUrl)
                println "Authentication successful! Logged in as: ${userInfo.userName}"

                // Generate PAT
                def pat = generatePAT(accessToken, apiUrl)

                // Add to shell config or display for manual addition
                if (shellConfigInfo.shell == "unknown") {
                    displayTokenForManualSetup(pat)
                } else {
                    addTokenToShellConfig(pat, shellConfigInfo)
                    println "- Personal Access Token generated and added to ${shellConfigInfo.configFile}"
                    println "Please restart your terminal or run: source ${shellConfigInfo.configFile}"
                }

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
                        output.println("<html><body><h1>Authentication successful!</h1><p>You can close this window.</p></body></html>")
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
                    Desktop.desktop.browse(new URI(url))
                } else {
                    // Fallback for systems without Desktop support
                    def os = System.getProperty("os.name").toLowerCase()
                    if (os.contains("win")) {
                        Runtime.runtime.exec("rundll32 url.dll,FileProtocolHandler ${url}")
                    } else if (os.contains("mac")) {
                        Runtime.runtime.exec("open ${url}")
                    } else {
                        Runtime.runtime.exec("xdg-open ${url}")
                    }
                }
            } catch (Exception e) {
                println "Could not open browser automatically. Please visit: ${url}"
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


        private boolean isCloudEndpoint(String apiUrl) {
            return apiUrl == 'https://api.cloud.seqera.io' ||
                   apiUrl == 'https://api.cloud.stage-seqera.io' ||
                   apiUrl == 'https://api.cloud.dev-seqera.io' ||
                   apiUrl == 'https://cloud.seqera.io/api' ||
                   apiUrl == 'https://cloud.stage-seqera.io/api' ||
                   apiUrl == 'https://cloud.dev-seqera.io/api'

        }

        private void handleEnterpriseAuth(String apiUrl) {
            println ""
            println "Please generate a Personal Access Token from your Seqera Platform instance."
            println "You can create one at: ${apiUrl.replace('/api', '').replace('api.', '')}/tokens"
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

            // Add PAT to shell config or display for manual addition
            def shellConfigInfo = detectShellConfig()
            if (shellConfigInfo.shell == "unknown") {
                displayTokenForManualSetup(pat.trim())
            } else {
                addTokenToShellConfig(pat.trim(), shellConfigInfo)
                println "Personal Access Token added to ${shellConfigInfo.configFile}"
                println "Please restart your terminal or run: source ${shellConfigInfo.configFile}"
            }
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

        private void displayTokenForManualSetup(String pat) {
            println ""
            println "=".repeat(60)
            println "MANUAL SETUP REQUIRED"
            println "=".repeat(60)
            println ""
            println "Please add the following line to your shell configuration file:"
            println ""
            println "export TOWER_ACCESS_TOKEN=${pat}"
            println ""
            println "Common shell config files:"
            println "  ~/.bashrc (bash)"
            println "  ~/.zshrc (zsh)"
            println "  ~/.config/fish/config.fish (fish)"
            println ""
            println "After adding the export, restart your terminal or run 'source <config-file>'"
            println "=".repeat(60)
        }

        private void addTokenToShellConfig(String pat, Map shellConfigInfo) {
            def configFile = Paths.get(shellConfigInfo.configFile as String)
            def commentLine = "# Seqera Platform access token"
            def exportLine = "export TOWER_ACCESS_TOKEN=\"${pat}\""

            // Append to existing file (we already checked it exists at command start)
            Files.writeString(configFile, "\n${commentLine}\n${exportLine}\n", java.nio.file.StandardOpenOption.APPEND)
        }

        @Override
        void usage(List<String> result) {
            result << 'Authenticate with Seqera Platform'
            result << "Usage: nextflow auth $name".toString()
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

            // Check if TOWER_ACCESS_TOKEN is set
            def existingToken = System.getenv("TOWER_ACCESS_TOKEN")
            if (!existingToken) {
                println "No TOWER_ACCESS_TOKEN environment variable found. Already logged out."
                return
            }

            // Find token in shell config
            def configFilePath = findTokenInShellConfig(existingToken)
            if (!configFilePath) {
                println "Token found in environment but not in shell config files."
                println "You may need to manually remove: export TOWER_ACCESS_TOKEN=${existingToken}"
                return
            }

            println "Found token in: ${configFilePath}"

            // Prompt user for API URL
            def apiUrl = promptForApiUrl()

            // Validate token by calling /user-info API
            try {
                def userInfo = callUserInfoApi(existingToken, apiUrl)
                println "Token is valid for user: ${userInfo.userName}"

                // Decode token to get token ID
                def tokenId = decodeTokenId(existingToken)
                println "Token ID: ${tokenId}"

                // Delete token via API
                deleteTokenViaApi(existingToken, apiUrl, tokenId)

                // Remove from shell config
                removeTokenFromShellConfig(existingToken, configFilePath)

                println "Successfully logged out and removed token."

            } catch (Exception e) {
                println "Failed to validate or delete token: ${e.message}"
                println "You may need to manually remove the export line from ${configFilePath}"
            }
        }

        private String promptForApiUrl() {
            System.out.print("Seqera Platform API URL [Default Seqera Cloud, https://api.cloud.seqera.io]: ")
            System.out.flush()

            def reader = new BufferedReader(new InputStreamReader(System.in))
            def input = reader.readLine()?.trim()

            return input?.isEmpty() || input == null ? 'https://api.cloud.seqera.io' : input
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

            println "Token successfully deleted from Seqera Platform."
        }

        private void removeTokenFromShellConfig(String token, String configFilePath) {
            def configFile = Paths.get(configFilePath)
            def content = Files.readString(configFile)

            // Look for the export line and optional comment above it
            def exportLine = "export TOWER_ACCESS_TOKEN=\"${token}\""
            def commentLine = "# Seqera Platform access token"

            def lines = content.split('\n').toList()
            def newLines = []
            def i = 0

            while (i < lines.size()) {
                def currentLine = lines[i]

                // Check if this is the export line we want to remove
                if (currentLine.trim() == exportLine.trim()) {
                    // Check if previous line is our comment
                    if (i > 0 && lines[i-1].trim() == commentLine.trim()) {
                        // Remove the comment line too (it was just added)
                        newLines.removeLast()
                    }
                    // Skip the export line (don't add it to newLines)
                } else {
                    newLines.add(currentLine)
                }
                i++
            }

            // Write back the modified content
            Files.writeString(configFile, newLines.join('\n'))
            println "Removed token export from ${configFilePath}"
        }

        @Override
        void usage(List<String> result) {
            result << 'Log out and remove Seqera Platform authentication'
            result << "Usage: nextflow auth $name".toString()
        }
    }
}
