package io.seqera.tower.plugin.cli

import groovy.json.JsonBuilder
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.cli.CmdAuth
import nextflow.cli.ColorUtil
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException

import java.awt.Desktop
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardOpenOption

@Slf4j
@CompileStatic
class AuthCommandImpl implements CmdAuth.AuthCommand {

    static final int API_TIMEOUT_MS = 10000
    static final int AUTH_POLL_TIMEOUT_RETRIES = 60
    static final int AUTH_POLL_INTERVAL_SECONDS = 5
    static final int WORKSPACE_SELECTION_THRESHOLD = 8  // Max workspaces to show in single list; above this uses org-first selection
    private static final String DEFAULT_API_ENDPOINT = 'https://api.cloud.seqera.io'

    private static final Map SEQERA_API_TO_AUTH0 = [
        'https://api.cloud.dev-seqera.io'  : [
            domain  : 'seqera-development.eu.auth0.com',
            clientId: 'Ep2LhYiYmuV9hhz0dH6dbXVq0S7s7SWZ'
        ],
        'https://api.cloud.stage-seqera.io': [
            domain  : 'seqera-stage.eu.auth0.com',
            clientId: '60cPDjI6YhoTPjyMTIBjGtxatSUwWswB'
        ],
        'https://api.cloud.seqera.io' : [
            domain  : 'seqera.eu.auth0.com',
            clientId: 'FxCM8EJ76nNeHUDidSHkZfT8VtsrhHeL'
        ]
    ]

    @Override
    void login(String apiUrl) {

        // Check if TOWER_ACCESS_TOKEN environment variable is set
        def envToken = System.getenv('TOWER_ACCESS_TOKEN')
        if( envToken ) {
            println ""
            ColorUtil.printColored("WARNING: Authentication token is already configured via TOWER_ACCESS_TOKEN environment variable.", "yellow bold")
            ColorUtil.printColored("${ColorUtil.colorize('nextflow auth login', 'cyan')} sets credentials using Nextflow config files, which take precedence over the environment variable.", "dim")
            ColorUtil.printColored(" however, caution is advised to avoid confusing behaviour.", "dim")
            println ""
        }

        // Check if .login file already exists
        def loginFile = getLoginFile()
        if( Files.exists(loginFile) ) {
            ColorUtil.printColored("Error: Authentication token is already configured in Nextflow config.", "red")
            ColorUtil.printColored("Login file: ${ColorUtil.colorize(loginFile.toString(), 'magenta')}", "dim")
            println " Run ${ColorUtil.colorize('nextflow auth logout', 'cyan')} to remove the current authentication."
            return
        }

        ColorUtil.printColored("Nextflow authentication with Seqera Platform", "cyan bold")
        ColorUtil.printColored(" - Authentication will be saved to: ${ColorUtil.colorize(getLoginFile().toString(), 'magenta')}", "dim")

        apiUrl = normalizeApiUrl(apiUrl)
        ColorUtil.printColored(" - Seqera Platform API endpoint: ${ColorUtil.colorize(apiUrl, 'magenta')} (can be customised with ${ColorUtil.colorize('-url', 'cyan')})", "dim")

        // Check if this is a cloud endpoint or enterprise
        def endpointInfo = getCloudEndpointInfo(apiUrl)
        if( endpointInfo.isCloud ) {
            try {
                performAuth0Login(endpointInfo.endpoint as String, endpointInfo.auth as Map)
            } catch( Exception e ) {
                log.debug("Authentication failed", e)
                println ""
                throw new AbortOperationException("${e.message}")
            }
        } else {
            // Enterprise endpoint - use PAT authentication
            handleEnterpriseAuth(apiUrl)
        }
    }

    private void performAuth0Login(String apiUrl, Map auth0Config) {

        // Start device authorization flow
        def deviceAuth = requestDeviceAuthorization(auth0Config)

        println ""
        ColorUtil.printColored("Confirmation code: ${ColorUtil.colorize(deviceAuth.user_code as String, 'yellow')}", "cyan bold")
        def urlWithCode = "${deviceAuth.verification_uri}?user_code=${deviceAuth.user_code}"
        println "${ColorUtil.colorize('Authentication URL:', 'cyan bold')} ${ColorUtil.colorize(urlWithCode, 'magenta')}"
        ColorUtil.printColored("\n[ Press Enter to open in browser ]", "cyan bold")
        System.in.read() // Wait for Enter key

        // Try to open browser automatically
        boolean browserOpened = openBrowser(urlWithCode)

        if( !browserOpened ) {
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

            // Generate PAT
            def pat = generatePAT(accessToken, apiUrl)

            // Save to config
            saveAuthToConfig(pat, apiUrl)

            // Automatically run configuration
            runConfiguration()

        } catch( Exception e ) {
            throw new RuntimeException("Authentication failed: ${e.message}", e)
        }
    }

    private boolean openBrowser(GString urlWithCode) {
        def browserOpened = false
        try {
            // Method 1: Java Desktop API
            if( Desktop.isDesktopSupported() ) {
                def desktop = Desktop.getDesktop()
                if( desktop.isSupported(Desktop.Action.BROWSE) ) {
                    desktop.browse(new URI(urlWithCode))
                    browserOpened = true
                }
            }

            // Method 2: Platform-specific commands
            if( !browserOpened ) {
                def os = System.getProperty("os.name").toLowerCase()
                def command = []

                if( os.contains("mac") || os.contains("darwin") ) {
                    command = ["open", urlWithCode]
                } else if( os.contains("win") ) {
                    command = ["cmd", "/c", "start", urlWithCode]
                } else {
                    // Linux and other Unix-like systems
                    def browsers = ["xdg-open", "firefox", "google-chrome", "chromium", "safari"]
                    for( browser in browsers ) {
                        try {
                            new ProcessBuilder(browser, urlWithCode).start()
                            browserOpened = true
                            break
                        } catch( Exception ignored ) {
                            // Try next browser
                        }
                    }
                }

                if( !browserOpened && command ) {
                    new ProcessBuilder(command as String[]).start()
                    browserOpened = true
                }
            }
        } catch( Exception e ) {
            log.debug("Exception opening browser", e)
        }
        return browserOpened
    }

    private void runConfiguration() {
        try {
            println ""
            // Just run the existing config command
            config()
        } catch( Exception e ) {
            ColorUtil.printColored("Configuration setup failed: ${e.message}", "red")
            ColorUtil.printColored("You can run 'nextflow auth config' later to set up your configuration.", "dim")
        }
    }

    private Map requestDeviceAuthorization(Map auth0Config) {
        def params = [
            'client_id': auth0Config.clientId,
            'scope'    : 'openid profile email offline_access',
            'audience' : 'platform'
        ]
        return performAuth0Request("https://${auth0Config.domain}/oauth/device/code", params)
    }

    private Map pollForDeviceToken(String deviceCode, int intervalSeconds, Map auth0Config) {
        def tokenUrl = "https://${auth0Config.domain}/oauth/token"
        def retryCount = 0

        while( retryCount < AUTH_POLL_TIMEOUT_RETRIES ) {
            def params = [
                'grant_type' : 'urn:ietf:params:oauth:grant-type:device_code',
                'device_code': deviceCode,
                'client_id'  : auth0Config.clientId
            ]

            try {
                def result = performAuth0Request(tokenUrl, params)
                return result
            } catch( RuntimeException e ) {
                def message = e.message
                if( message.contains('authorization_pending') ) {
                    print "${ColorUtil.colorize('.', 'dim', true)}"
                    System.out.flush()
                } else if( message.contains('slow_down') ) {
                    intervalSeconds += 5
                    print "${ColorUtil.colorize('.', 'dim', true)}"
                    System.out.flush()
                } else if( message.contains('expired_token') ) {
                    throw new RuntimeException("The device code has expired. Please try again.")
                } else if( message.contains('access_denied') ) {
                    throw new RuntimeException("Access denied by user")
                } else {
                    throw e
                }
            }

            Thread.sleep(intervalSeconds * 1000)
            retryCount++
        }

        throw new RuntimeException("Authentication timed out. Please try again.")
    }


    private void handleEnterpriseAuth(String apiUrl) {
        println ""
        ColorUtil.printColored("Please generate a Personal Access Token from your Seqera Platform instance.", "cyan bold")
        println "You can create one at: ${ColorUtil.colorize(apiUrl.replace('/api', '').replace('://api.', '') + '/tokens', 'magenta')}"
        println ""

        System.out.print("Enter your Personal Access Token: ")
        System.out.flush()

        def console = System.console()
        def pat = console ?
            new String(console.readPassword()) :
            new BufferedReader(new InputStreamReader(System.in)).readLine()

        if( !pat || pat.trim().isEmpty() ) {
            throw new AbortOperationException("Personal Access Token is required for Seqera Platform Enterprise authentication")
        }

        // Save to config
        saveAuthToConfig(pat.trim(), apiUrl)
        ColorUtil.printColored("Personal Access Token saved to Nextflow login config (${getLoginFile().toString()})", "green")
    }

    private String generatePAT(String accessToken, String apiUrl) {
        def tokensUrl = "${apiUrl}/tokens"
        def username = System.getProperty("user.name")
        def timestamp = new Date().format("yyyy-MM-dd-HH-mm")
        def tokenName = "nextflow-auth-${username}-${timestamp}"

        def requestBody = new JsonBuilder([name: tokenName]).toString()

        def connection = createHttpConnection(tokensUrl, 'POST', accessToken)
        connection.setRequestProperty('Content-Type', 'application/json')
        connection.doOutput = true

        connection.outputStream.withWriter { writer ->
            writer.write(requestBody)
        }

        if( connection.responseCode != 200 ) {
            def error = connection.errorStream?.text ?: "HTTP ${connection.responseCode}"
            throw new RuntimeException("Failed to generate PAT: ${error}")
        }

        def response = connection.inputStream.text
        def json = new JsonSlurper().parseText(response) as Map
        return json.accessKey as String
    }

    private String normalizeApiUrl(String url) {
        if( !url ) {
            return DEFAULT_API_ENDPOINT
        }
        if( !url.startsWith('http://') && !url.startsWith('https://') ) {
            return 'https://' + url
        }
        return url
    }

    private Map performAuth0Request(String url, Map params) {
        def postData = params.collect { k, v ->
            "${URLEncoder.encode(k.toString(), 'UTF-8')}=${URLEncoder.encode(v.toString(), 'UTF-8')}"
        }.join('&')

        def connection = new URL(url).openConnection() as HttpURLConnection
        connection.requestMethod = 'POST'
        connection.connectTimeout = API_TIMEOUT_MS
        connection.readTimeout = API_TIMEOUT_MS
        connection.setRequestProperty('Content-Type', 'application/x-www-form-urlencoded')
        connection.doOutput = true

        connection.outputStream.withWriter { writer ->
            writer.write(postData)
        }

        if( connection.responseCode == 200 ) {
            def response = connection.inputStream.text
            def json = new JsonSlurper().parseText(response)
            return json as Map
        } else {
            def errorResponse = connection.errorStream?.text
            if( errorResponse ) {
                def errorJson = new JsonSlurper().parseText(errorResponse) as Map
                def error = errorJson.error
                throw new RuntimeException("${error}: ${errorJson.error_description ?: ''}")
            } else {
                throw new RuntimeException("Request failed: HTTP ${connection.responseCode}")
            }
        }
    }

    private void saveAuthToConfig(String accessToken, String apiUrl) {
        def config = [:]
        config['tower.accessToken'] = accessToken
        config['tower.endpoint'] = apiUrl
        config['tower.enabled'] = true

        writeConfig(config, null)
    }

    @Override
    void logout() {
        // Check if .login file exists
        def loginFile = getLoginFile()
        if( !Files.exists(loginFile) ) {
            ColorUtil.printColored("No previous login found.", "green")
            return
        }
        // Read token from .login file
        def loginConfig = readLoginFile()
        def existingToken = loginConfig['tower.accessToken']
        def apiUrl = loginConfig['tower.endpoint'] as String ?: DEFAULT_API_ENDPOINT

        if( !existingToken ) {
            ColorUtil.printColored("WARN: No authentication token found in login file.", "yellow bold")
            println "Removing file: ${ColorUtil.colorize(loginFile.toString(), 'magenta')}"
            removeAuthFromConfig()
            return
        }

        // Check if TOWER_ACCESS_TOKEN environment variable is set
        def envToken = System.getenv('TOWER_ACCESS_TOKEN')
        if( envToken ) {
            println ""
            ColorUtil.printColored("WARNING: TOWER_ACCESS_TOKEN environment variable is set.", "yellow bold")
            println " ${ColorUtil.colorize('nextflow auth logout', 'dim cyan')}${ColorUtil.colorize(' only removes credentials from Nextflow config files.', 'dim')}"
            ColorUtil.printColored(" The environment variable will remain unaffected.", "dim")
            println ""
        }

        ColorUtil.printColored(" - Found authentication token in login file: ${ColorUtil.colorize(loginFile.toString(), 'magenta')}", "dim")
        ColorUtil.printColored(" - Using Seqera Platform endpoint: ${ColorUtil.colorize(apiUrl, 'magenta')}", "dim")

        // Validate token by calling /user-info API
        try {
            def userInfo = callUserInfoApi(existingToken as String, apiUrl)
            ColorUtil.printColored(" - Token is valid for user: ${ColorUtil.colorize(userInfo.userName as String, 'cyan bold')}", "dim")

            // Only delete PAT from platform if this is a cloud endpoint
            if( isCloudEndpoint(apiUrl) ) {
                def tokenId = decodeTokenId(existingToken as String)
                deleteTokenViaApi(existingToken as String, apiUrl, tokenId)
            } else {
                println " - Enterprise installation detected - PAT will not be deleted from platform."
            }
            removeAuthFromConfig()

        } catch( Exception e ) {
            println "Failed to validate or delete token: ${e.message}"
            println "Removing token from config anyway..."

            // Remove from config even if API calls fail
            removeAuthFromConfig()
        }
    }

    private String decodeTokenId(String token) {
        try {
            // Decode base64 token
            def decoded = new String(Base64.decoder.decode(token), "UTF-8")

            // Parse JSON to extract token ID
            def json = new JsonSlurper().parseText(decoded) as Map
            def tokenId = json.tid

            if( !tokenId ) {
                throw new RuntimeException("No token ID found in decoded token")
            }

            return tokenId.toString()
        } catch( Exception e ) {
            throw new RuntimeException("Failed to decode token ID: ${e.message}")
        }
    }

    private void deleteTokenViaApi(String token, String apiUrl, String tokenId) {
        def connection = createHttpConnection("${apiUrl}/tokens/${tokenId}", 'DELETE', token)

        if( connection.responseCode != 200 && connection.responseCode != 204 ) {
            def error = connection.errorStream?.text ?: "HTTP ${connection.responseCode}"
            throw new RuntimeException("Failed to delete token: ${error}")
        }

        ColorUtil.printColored("\nToken successfully deleted from Seqera Platform.", "green")
    }

    private void removeAuthFromConfig() {
        def configFile = getConfigFile()
        def loginFile = getLoginFile()

        // Remove includeConfig line from main config file
        if( Files.exists(configFile) ) {
            def existingContent = Files.readString(configFile)
            def updatedContent = removeIncludeConfigLine(existingContent)
            Files.writeString(configFile, updatedContent.toString(), StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)
        }

        // Delete .login file
        if( Files.exists(loginFile) ) {
            Files.delete(loginFile)
        }

        ColorUtil.printColored("Authentication removed from Nextflow config.", "green")
    }

    @Override
    void config() {
        // Read from both main config and .login file
        def config = readConfig()

        // Token can come from .login file or environment variable
        def existingToken = config['tower.accessToken'] ?: System.getenv('TOWER_ACCESS_TOKEN')
        def endpoint = config['tower.endpoint'] ?: DEFAULT_API_ENDPOINT

        if( !existingToken ) {
            println "No authentication found. Please run ${ColorUtil.colorize('nextflow auth login', 'cyan')} first."
            return
        }

        ColorUtil.printColored("Nextflow Seqera Platform configuration", "cyan bold")
        ColorUtil.printColored(" - Config file: ${ColorUtil.colorize(getLoginFile().toString(), 'magenta')}", "dim")

        // Check if token is from environment variable
        if( !config['tower.accessToken'] && System.getenv('TOWER_ACCESS_TOKEN') ) {
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
            def workspaceResult = configureWorkspace(config, existingToken as String, endpoint as String, userInfo.id as String)
            configChanged = configChanged || (workspaceResult.changed as boolean)

            // Save updated config only if changes were made
            if( configChanged ) {
                writeConfig(config, workspaceResult.metadata as Map)

                // Show the new configuration
                println "\nNew configuration:"
                showCurrentConfig(config, existingToken as String, endpoint as String)

                ColorUtil.printColored("\nConfiguration saved to ${ColorUtil.colorize(getLoginFile().toString(), 'magenta')}", "green")
            } else {
                ColorUtil.printColored("\nNo configuration changes were made.", "dim")
            }

        } catch( Exception e ) {
            throw new AbortOperationException("Failed to configure settings: ${e.message}")
        }
    }

    private void showCurrentConfig(Map config, String accessToken, String endpoint) {
        // Show workflow monitoring status
        def monitoringEnabled = config.get('tower.enabled', false)
        println "  ${ColorUtil.colorize('Workflow monitoring:', 'cyan')} ${monitoringEnabled ? ColorUtil.colorize('enabled', 'green') : ColorUtil.colorize('disabled', 'red')}"

        // Show workspace setting
        def workspaceId = config.get('tower.workspaceId')
        if( workspaceId ) {
            // Try to get workspace details from API for display
            def workspaceDetails = getWorkspaceDetailsFromApi(accessToken, endpoint, workspaceId as String)
            if( workspaceDetails ) {
                def details = workspaceDetails as Map
                println "  ${ColorUtil.colorize('Default workspace:', 'cyan')} ${ColorUtil.colorize(workspaceId as String, 'magenta')} ${ColorUtil.colorize('[' + details.orgName + ' / ' + details.workspaceName + ']', 'dim', true)}"
            } else {
                println "  ${ColorUtil.colorize('Default workspace:', 'cyan')} ${ColorUtil.colorize(workspaceId as String, 'magenta')}"
            }
        } else {
            println "  ${ColorUtil.colorize('Default workspace:', 'cyan')} ${ColorUtil.colorize('None (Personal workspace)', 'magenta')}"
        }

        // Show API endpoint
        def apiEndpoint = config.get('tower.endpoint') ?: endpoint
        println "  ${ColorUtil.colorize('API endpoint:', 'cyan', true)} $apiEndpoint"
    }

    private boolean configureEnabled(Map config) {
        def currentEnabled = config.get('tower.enabled', false)

        println "Workflow monitoring settings. Current setting: ${currentEnabled ? ColorUtil.colorize('enabled', 'green') : ColorUtil.colorize('disabled', 'red')}"
        ColorUtil.printColored("  When enabled, all workflow runs are automatically monitored by Seqera Platform", "dim")
        ColorUtil.printColored("  When disabled, you can enable per-run with the ${ColorUtil.colorize('-with-tower', 'cyan')} flag", "dim")
        println ""

        def promptText = "${ColorUtil.colorize('Enable workflow monitoring for all runs?', 'cyan bold', true)} (${currentEnabled ? ColorUtil.colorize('Y', 'green') + '/n' : 'y/' + ColorUtil.colorize('N', 'red')}): "
        def input = promptForYesNo(promptText, currentEnabled)

        if( input == null ) {
            return false // No change
        } else if( input && !currentEnabled ) {
            config['tower.enabled'] = true
            return true
        } else if( !input && currentEnabled ) {
            config.remove('tower.enabled')
            return true
        }
        return false
    }

    private Map configureWorkspace(Map config, String accessToken, String endpoint, String userId) {
        // Check if TOWER_WORKFLOW_ID environment variable is set
        def envWorkspaceId = System.getenv('TOWER_WORKFLOW_ID')
        if( envWorkspaceId ) {
            println "\nDefault workspace: ${ColorUtil.colorize('TOWER_WORKFLOW_ID environment variable is set', 'yellow')}"
            ColorUtil.printColored("  Not prompting for default workspace configuration as environment variable takes precedence", "dim")
            return [changed: false, metadata: null]
        }

        // Get all workspaces for the user
        def workspaces = getUserWorkspaces(accessToken, endpoint, userId)

        if( !workspaces ) {
            println "\nNo workspaces found for your account."
            return [changed: false, metadata: null]
        }

        // Show current workspace setting
        def currentWorkspaceId = config.get('tower.workspaceId')
        def currentWorkspace = workspaces.find { ((Map) it).workspaceId.toString() == currentWorkspaceId?.toString() }

        def currentSetting
        if( currentWorkspace ) {
            def workspace = currentWorkspace as Map
            currentSetting = "${workspace.orgName} / ${workspace.workspaceName}"
        } else {
            currentSetting = "None (Personal workspace)"
        }

        println "\nDefault workspace. Current setting: ${ColorUtil.colorize(currentSetting as String, 'cyan', true)}"
        ColorUtil.printColored("  Workflow runs use this workspace by default", "dim")
        // Group by organization
        def orgWorkspaces = workspaces.groupBy { ((Map) it).orgName ?: 'Personal' }

        // If threshold or fewer total options, show all at once
        if( workspaces.size() <= WORKSPACE_SELECTION_THRESHOLD ) {
            return selectWorkspaceFromAll(config, workspaces, currentWorkspaceId)
        } else {
            // Two-stage selection: org first, then workspace
            return selectWorkspaceByOrg(config, orgWorkspaces, currentWorkspaceId)
        }
    }

    private Map selectWorkspaceFromAll(Map config, List workspaces, def currentWorkspaceId) {
        println "\nAvailable workspaces:"
        println "  0. ${ColorUtil.colorize('None (Personal workspace)', 'cyan', true)} ${ColorUtil.colorize('[no organization]', 'dim', true)}"

        workspaces.eachWithIndex { workspace, index ->
            def ws = workspace as Map
            def prefix = ws.orgName ? "${ColorUtil.colorize(ws.orgName as String, 'cyan', true)} / " : ""
            println "  ${index + 1}. ${prefix}${ColorUtil.colorize(ws.workspaceName as String, 'magenta', true)} ${ColorUtil.colorize('[' + (ws.workspaceFullName as String) + ']', 'dim', true)}"
        }

        // Show current workspace and prepare prompt
        def currentWorkspace = workspaces.find { ((Map) it).workspaceId.toString() == currentWorkspaceId?.toString() }
        def currentWorkspaceName
        if( currentWorkspace ) {
            def workspace = currentWorkspace as Map
            currentWorkspaceName = "${workspace.orgName} / ${workspace.workspaceName}"
        } else {
            currentWorkspaceName = "None (Personal workspace)"
        }

        def prompt = "\nSelect workspace (0-${workspaces.size()}, press Enter to keep as '${currentWorkspaceName}'): "
        def selection = promptForNumber(prompt, 0, workspaces.size(), true)

        if( selection == null ) {
            return [changed: false, metadata: null]
        }

        if( selection == 0 ) {
            def hadWorkspaceId = config.containsKey('tower.workspaceId')
            config.remove('tower.workspaceId')
            return [changed: hadWorkspaceId, metadata: null]
        } else {
            def selectedWorkspace = workspaces[selection - 1] as Map
            def selectedId = selectedWorkspace.workspaceId.toString()
            def currentId = config.get('tower.workspaceId')
            config['tower.workspaceId'] = selectedId
            def metadata = [
                orgName          : selectedWorkspace.orgName,
                workspaceName    : selectedWorkspace.workspaceName,
                workspaceFullName: selectedWorkspace.workspaceFullName
            ]
            return [changed: currentId != selectedId, metadata: metadata]
        }
    }

    private Map selectWorkspaceByOrg(Map config, Map orgWorkspaces, def currentWorkspaceId) {
        // Get current workspace info for prompts
        def allWorkspaces = []
        orgWorkspaces.values().each { workspaceList ->
            allWorkspaces.addAll(workspaceList as List)
        }
        def currentWorkspace = allWorkspaces.find { ((Map) it).workspaceId.toString() == currentWorkspaceId?.toString() }
        def currentWorkspaceName
        def currentWorkspaceDisplay
        if( currentWorkspace ) {
            def workspace = currentWorkspace as Map
            currentWorkspaceName = workspace.workspaceName as String
            currentWorkspaceDisplay = "${workspace.orgName} / ${workspace.workspaceName}"
        } else {
            currentWorkspaceName = null
            currentWorkspaceDisplay = "None"
        }

        // First, select organization
        def orgs = orgWorkspaces.keySet().toList()

        // Always add Personal as first option (it's never returned by the API but should always be available)
        orgs.add(0, 'Personal')

        println "\nAvailable organizations:"
        orgs.eachWithIndex { orgName, index ->
            def displayName = orgName == 'Personal' ? 'None [Personal workspace]' : orgName
            println "  ${index + 1}. ${ColorUtil.colorize(displayName as String, 'cyan', true)}"
        }
        System.out.print("${ColorUtil.colorize("Select organization (1-${orgs.size()}, leave blank to keep as '${currentWorkspaceDisplay}'): ", 'dim', true)}")
        System.out.flush()

        def reader = new BufferedReader(new InputStreamReader(System.in))
        def orgSelection
        while( true ) {
            def orgInput = reader.readLine()?.trim()

            if( orgInput.isEmpty() ) {
                return [changed: false, metadata: null]
            }

            try {
                orgSelection = Integer.parseInt(orgInput)
                if( orgSelection >= 1 && orgSelection <= orgs.size() ) {
                    break
                }
            } catch( NumberFormatException ignored ) {
                // Fall through to error message
            }
            ColorUtil.printColored("Invalid input. Please enter a number between 1 and ${orgs.size()}.", "red")
            System.out.print("${ColorUtil.colorize("Select organization (1-${orgs.size()}, leave blank to keep as '${currentWorkspaceDisplay}'): ", 'dim', true)}")
            System.out.flush()
        }

        def selectedOrgName = orgs[orgSelection - 1]

        // If Personal was selected, remove workspace ID (use personal workspace)
        if( selectedOrgName == 'Personal' ) {
            def hadWorkspaceId = config.containsKey('tower.workspaceId')
            config.remove('tower.workspaceId')
            return [changed: hadWorkspaceId, metadata: null]
        }

        def orgWorkspaceList = orgWorkspaces[selectedOrgName] as List

        println ""
        println "Select workspace in ${selectedOrgName}:"

        orgWorkspaceList.eachWithIndex { workspace, index ->
            def ws = workspace as Map
            println "  ${index + 1}. ${ColorUtil.colorize(ws.workspaceName as String, 'magenta', true)} ${ColorUtil.colorize('[' + (ws.workspaceFullName as String) + ']', 'dim', true)}"
        }

        def maxSelection = orgWorkspaceList.size()
        def wsSelection
        while( true ) {
            System.out.print("${ColorUtil.colorize("Select workspace (1-${maxSelection}): ", 'dim', true)}")
            System.out.flush()

            def wsInput = reader.readLine()?.trim()
            if( wsInput.isEmpty() ) {
                ColorUtil.printColored("Please enter a selection.", "red")
                continue
            }

            try {
                wsSelection = Integer.parseInt(wsInput)
                if( wsSelection >= 1 && wsSelection <= maxSelection ) {
                    break
                }
            } catch( NumberFormatException ignored ) {
                // Fall through to error message
            }
            ColorUtil.printColored("Invalid input. Please enter a number between 1 and ${maxSelection}.", "red")
        }

        def selectedWorkspace = orgWorkspaceList[wsSelection - 1] as Map
        def selectedId = selectedWorkspace.workspaceId.toString()
        def currentId = config.get('tower.workspaceId')
        config['tower.workspaceId'] = selectedId
        def metadata = [
            orgName          : selectedWorkspace.orgName,
            workspaceName    : selectedWorkspace.workspaceName,
            workspaceFullName: selectedWorkspace.workspaceFullName
        ]
        return [changed: currentId != selectedId, metadata: metadata]
    }

    private List getUserWorkspaces(String accessToken, String endpoint, String userId) {
        def connection = createHttpConnection("${endpoint}/user/${userId}/workspaces", 'GET', accessToken)

        if( connection.responseCode != 200 ) {
            def error = connection.errorStream?.text ?: "HTTP ${connection.responseCode}"
            throw new RuntimeException("Failed to get workspaces: ${error}")
        }

        def response = connection.inputStream.text
        def json = new JsonSlurper().parseText(response) as Map
        def orgsAndWorkspaces = json.orgsAndWorkspaces as List

        return orgsAndWorkspaces.findAll { ((Map) it).workspaceId != null }
    }

    private Boolean promptForYesNo(String prompt, Boolean defaultValue) {
        while( true ) {
            System.out.print(prompt)
            System.out.flush()
            def input = readUserInput()?.toLowerCase()

            if( input?.isEmpty() ) {
                return defaultValue
            } else if( input in ['y', 'yes'] ) {
                return true
            } else if( input in ['n', 'no'] ) {
                return false
            } else {
                ColorUtil.printColored("Invalid input. Please enter 'y', 'n', or press Enter to keep current setting.", "red")
            }
        }
    }

    private Integer promptForNumber(String prompt, int min, int max, boolean allowEmpty = false) {
        while( true ) {
            ColorUtil.printColored(prompt, "bold cyan")
            System.out.flush()
            def input = readUserInput()

            if( input?.isEmpty() && allowEmpty ) {
                return null
            }

            try {
                def number = Integer.parseInt(input)
                if( number >= min && number <= max ) {
                    return number
                }
            } catch( NumberFormatException e ) {
                // Fall through to error message
            }
            ColorUtil.printColored("Invalid input. Please enter a number between ${min} and ${max}.", "red")
        }
    }

    @Override
    void status() {
        def config = readConfig()
        def loginConfig = readLoginFile()
        // Collect all status information
        List<List<String>> statusRows = []
        def workspaceDisplayInfo = null

        // API endpoint
        def endpointInfo = getConfigValue(config, loginConfig, 'tower.endpoint', 'TOWER_API_ENDPOINT', DEFAULT_API_ENDPOINT)
        statusRows.add(['API endpoint', ColorUtil.colorize(endpointInfo.value as String, 'magenta'), endpointInfo.source as String])

        // API connection check
        def apiConnectionOk = checkApiConnection(endpointInfo.value as String)
        def connectionColor = apiConnectionOk ? 'green' : 'red'
        statusRows.add(['API connection', ColorUtil.colorize(apiConnectionOk ? 'OK' : 'ERROR', connectionColor), ''])

        // Authentication check
        def tokenInfo = getConfigValue(config, loginConfig, 'tower.accessToken', 'TOWER_ACCESS_TOKEN')
        if( tokenInfo.value ) {
            try {
                def userInfo = callUserInfoApi(tokenInfo.value as String, endpointInfo.value as String)
                def currentUser = userInfo.userName as String
                statusRows.add(['Authentication', "${ColorUtil.colorize('OK', 'green')} (user: ${ColorUtil.colorize(currentUser, 'cyan')})".toString(), tokenInfo.source as String])
            } catch( Exception e ) {
                statusRows.add(['Authentication', ColorUtil.colorize('ERROR', 'red'), 'failed'])
            }
        } else {
            statusRows.add(['Authentication', "${ColorUtil.colorize('ERROR', 'red')} ${ColorUtil.colorize('(no token)', 'dim')}".toString(), 'not set'])
        }

        // Monitoring enabled
        def enabledInfo = getConfigValue(config, loginConfig,'tower.enabled', null, 'false')
        def enabledValue = enabledInfo.value?.toString()?.toLowerCase() in ['true', '1', 'yes'] ? 'Yes' : 'No'
        def enabledColor = enabledValue == 'Yes' ? 'green' : 'yellow'
        statusRows.add(['Workflow monitoring', ColorUtil.colorize(enabledValue, enabledColor), (enabledInfo.source ?: 'default') as String])

        // Default workspace
        def workspaceInfo = getConfigValue(config, loginConfig, 'tower.workspaceId', 'TOWER_WORKFLOW_ID')
        if( workspaceInfo.value ) {
            // Try to get workspace name from API if we have a token
            def workspaceDetails = null
            if( tokenInfo.value ) {
                workspaceDetails = getWorkspaceDetailsFromApi(tokenInfo.value as String, endpointInfo.value as String, workspaceInfo.value as String)
            }

            if( workspaceDetails ) {
                // Add workspace ID row
                statusRows.add(['Default workspace', ColorUtil.colorize(workspaceInfo.value as String, 'blue'), workspaceInfo.source as String])
                // Store workspace details for display after the table
                workspaceDisplayInfo = workspaceDetails
            } else {
                statusRows.add(['Default workspace', ColorUtil.colorize(workspaceInfo.value as String, 'blue', true), workspaceInfo.source as String])
            }
        } else {
            statusRows.add(['Default workspace', ColorUtil.colorize('None (Personal workspace)', 'cyan', true), 'default'])
        }

        // Print table
        println ""
        printStatusTable(statusRows)

        // Print workspace details if available
        if( workspaceDisplayInfo ) {
            println "${' ' * 22}${ColorUtil.colorize(workspaceDisplayInfo.orgName as String, 'cyan')} / ${ColorUtil.colorize(workspaceDisplayInfo.workspaceName as String, 'cyan')} ${ColorUtil.colorize('[' + (workspaceDisplayInfo.workspaceFullName as String) + ']', 'dim')}"
        }
    }

    private void printStatusTable(List<List<String>> rows) {
        if( !rows ) return

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
        if( path.startsWith(userHome) ) {
            return '~' + path.substring(userHome.length())
        }
        return path
    }

    private Map getConfigValue(Map config, Map login, String configKey, String envVarName, String defaultValue = null) {
        //Checks where the config value came from
        def loginValue = login[configKey]
        def configValue = config[configKey]
        def envValue = envVarName ? System.getenv(envVarName) : null
        def effectiveValue = configValue ?: envValue ?: defaultValue

        def source = null
        if( loginValue ) {
            source = shortenPath(getLoginFile().toString())
        } else if( configValue ) {
            source = shortenPath(getConfigFile().toString())
        } else if( envValue ) {
            source = "env var \$${envVarName}"
        } else if( defaultValue ) {
            source = "default"
        }

        return [
            value     : effectiveValue,
            source    : source,
            fromConfig: configValue != null,
            fromEnv   : envValue != null,
            isDefault : !configValue && !envValue
        ]
    }

    private String getWorkspaceNameFromApi(String accessToken, String endpoint, String workspaceId) {
        try {
            def workspaceDetails = getWorkspaceDetailsFromApi(accessToken, endpoint, workspaceId)
            if( workspaceDetails ) {
                return "${workspaceDetails.orgName} / ${workspaceDetails.workspaceName} [${workspaceDetails.workspaceFullName}]"
            }
            return null
        } catch( Exception e ) {
            return null
        }
    }

    private boolean checkApiConnection(String endpoint) {
        try {
            def connection = new URL("${endpoint}/service-info").openConnection() as HttpURLConnection
            connection.requestMethod = 'GET'
            connection.connectTimeout = API_TIMEOUT_MS
            connection.readTimeout = API_TIMEOUT_MS
            return connection.responseCode == 200
        } catch( Exception e ) {
            return false
        }
    }

    private String promptForApiUrl() {
        System.out.print("Seqera Platform API endpoint [Default ${DEFAULT_API_ENDPOINT}]: ")
        System.out.flush()

        def input = readUserInput()
        return input?.isEmpty() ? DEFAULT_API_ENDPOINT : input
    }

    private static String readUserInput() {
        def reader = new BufferedReader(new InputStreamReader(System.in))
        return reader.readLine()?.trim()
    }

    private HttpURLConnection createHttpConnection(String url, String method, String authToken = null) {
        def connection = new URL(url).openConnection() as HttpURLConnection
        connection.requestMethod = method
        connection.connectTimeout = API_TIMEOUT_MS
        connection.readTimeout = API_TIMEOUT_MS
        if (authToken) {
            connection.setRequestProperty('Authorization', "Bearer ${authToken}")
        }
        return connection
    }

    private Map getWorkspaceDetailsFromApi(String accessToken, String endpoint, String workspaceId) {
        try {
            def userInfo = callUserInfoApi(accessToken, endpoint)
            def userId = userInfo.id as String

            def connection = createHttpConnection("${endpoint}/user/${userId}/workspaces", 'GET', accessToken)

            if (connection.responseCode != 200) {
                return null
            }

            def response = connection.inputStream.text
            def json = new JsonSlurper().parseText(response) as Map
            def orgsAndWorkspaces = json.orgsAndWorkspaces as List

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

    private Map getCloudEndpointInfo(String apiUrl) {
        for (env in SEQERA_API_TO_AUTH0) {
            def standardUrl = env.key as String
            def legacyUrl = standardUrl.replace('://api.', '://') + '/api'
            if (apiUrl == standardUrl || apiUrl == legacyUrl) {
                return [isCloud: true, endpoint: env.key, auth: env.value]
            }
        }
        return [isCloud: false, endpoint: apiUrl, auth: null]
    }

    private boolean isCloudEndpoint(String apiUrl) {
        return getCloudEndpointInfo(apiUrl).isCloud
    }

    private Map callUserInfoApi(String accessToken, String apiUrl) {
        def connection = createHttpConnection("${apiUrl}/user-info", 'GET', accessToken)

        if (connection.responseCode != 200) {
            def error = connection.errorStream?.text ?: "HTTP ${connection.responseCode}"
            throw new RuntimeException("Failed to get user info: ${error}")
        }

        def response = connection.inputStream.text
        def json = new JsonSlurper().parseText(response) as Map
        return json.user as Map
    }

    private Path getConfigFile() {
        return Const.APP_HOME_DIR.resolve('config')
    }

    private Path getLoginFile() {
        return Const.APP_HOME_DIR.resolve('.login')
    }

    private Map readLoginFile() {
        def configFile = getLoginFile()
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

    private Map readConfig() {
        def builder = new ConfigBuilder().setHomeDir(Const.APP_HOME_DIR).setCurrentDir(Const.APP_HOME_DIR)
        return builder.buildConfigObject().flatten()
    }

    private void writeConfig(Map config, Map workspaceMetadata = null) {
        def configFile = getConfigFile()
        def loginFile = getLoginFile()

        // Create directory if it doesn't exist
        if (!Files.exists(configFile.parent)) {
            Files.createDirectories(configFile.parent)
        }

        // Write tower config to .login file
        def towerConfig = config.findAll { key, value ->
            key.toString().startsWith('tower.') && !key.toString().endsWith('.comment')
        }

        def loginConfigText = new StringBuilder()
        loginConfigText.append("// Seqera Platform configuration\n")
        loginConfigText.append("tower {\n")
        towerConfig.each { key, value ->
            def configKey = key.toString().substring(6) // Remove "tower." prefix

            if (value instanceof String) {
                def line = "    ${configKey} = '${value}'"
                // Add workspace comment if this is workspaceId and we have metadata
                if (configKey == 'workspaceId' && workspaceMetadata) {
                    line += "  // ${workspaceMetadata.orgName} / ${workspaceMetadata.workspaceName} [${workspaceMetadata.workspaceFullName}]"
                }
                loginConfigText.append("${line}\n")
            } else {
                loginConfigText.append("    ${configKey} = ${value}\n")
            }
        }
        loginConfigText.append("}\n")

        // Write the .login file
        Files.writeString(loginFile, loginConfigText.toString(), StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)

        // Add includeConfig line to main config file if it doesn't exist
        addIncludeConfigToMainFile(configFile)
    }

    private void addIncludeConfigToMainFile(Path configFile) {
        def includeConfigLine = "includeConfig '.login'"

        def configContent = ""
        if (Files.exists(configFile)) {
            configContent = Files.readString(configFile)
            // Check if includeConfig line already exists
            if (configContent.contains(includeConfigLine)) {
                return // Already exists, nothing to do
            }
        }

        // Append the includeConfig line
        def updatedContent = configContent
        if (!updatedContent.isEmpty() && !updatedContent.endsWith("\n")) {
            updatedContent += "\n"
        }
        updatedContent += "\n${includeConfigLine}\n"

        Files.writeString(configFile, updatedContent, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)
    }

    private String removeIncludeConfigLine(String content) {
        // Remove the includeConfig '.login' line
        def lines = content.split('\n')
        def filteredLines = lines.findAll { line ->
            !line.trim().equals("includeConfig '.login'")
        }
        return filteredLines.join('\n')
    }
}
