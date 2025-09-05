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
import groovy.yaml.YamlBuilder
import groovy.yaml.YamlSlurper
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

            println "Nextflow authentication with Seqera Platform"
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

        private String promptForApiUrl() {
            System.out.print("Seqera Platform API URL [Default Seqera Cloud, https://api.cloud.seqera.io]: ")
            System.out.flush()

            def reader = new BufferedReader(new InputStreamReader(System.in))
            def input = reader.readLine()?.trim()

            return input?.isEmpty() || input == null ? 'https://api.cloud.seqera.io' : input
        }

        private void performAuth0Login(String apiUrl) {
            println "Opening browser for authentication..."

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

            println "Waiting for authentication to complete..."

            try {
                // Wait for callback with timeout
                def authCode = callbackFuture.get(5, TimeUnit.MINUTES)

                // Exchange code for token
                def tokenData = exchangeCodeForToken(authCode, codeVerifier)

                // Save credentials
                saveCredentials('oauth', tokenData['access_token'], apiUrl)

                println "Authentication successful! Credentials saved to ${getCredentialsPath()}"

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

        private void saveCredentials(String type, String token, String apiUrl) {
            def xdgConfigHome = System.getenv("XDG_CONFIG_HOME")
            def credentialsDir = xdgConfigHome ?
                Paths.get(xdgConfigHome, "seqera") :
                Paths.get(System.getProperty("user.home"), ".config", "seqera")
            Files.createDirectories(credentialsDir)

            def credentialsFile = credentialsDir.resolve("credentials.yml")

            def config = [:]
            if (Files.exists(credentialsFile)) {
                try {
                    def yaml = new YamlSlurper()
                    config = yaml.parse(credentialsFile.toFile()) as Map ?: [:]
                } catch (Exception e) {
                    log.debug("Could not parse existing credentials file", e)
                    config = [:]
                }
            }

            config['default'] = [
                'type': type,
                'token': token,
                'endpoint': apiUrl
            ]

            def yaml = new YamlBuilder()
            yaml(config)

            Files.write(credentialsFile, yaml.toString().getBytes("UTF-8"))

            // Set file permissions
            try {
                def perms = java.nio.file.attribute.PosixFilePermissions.fromString("rw-r--r--")
                Files.setPosixFilePermissions(credentialsFile, perms)
            } catch (Exception e) {
                log.debug("Could not set file permissions", e)
            }
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

            // Save PAT credentials
            saveCredentials('pat', pat.trim(), apiUrl)
            println "Authentication successful! Credentials saved to ${getCredentialsPath()}"
        }


        private String getCredentialsPath() {
            def xdgConfigHome = System.getenv("XDG_CONFIG_HOME")
            def credentialsDir = xdgConfigHome ?
                Paths.get(xdgConfigHome, "seqera") :
                Paths.get(System.getProperty("user.home"), ".config", "seqera")
            return credentialsDir.resolve("credentials.yml").toString()
        }

        @Override
        void usage(List<String> result) {
            result << 'Authenticate with Seqera Platform'
            result << "Usage: nextflow auth $name".toString()
        }
    }
}
