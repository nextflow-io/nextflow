package io.seqera.tower.plugin.auth

import groovy.json.JsonBuilder
import groovy.json.JsonSlurper
import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.http.HxClient
import nextflow.Const
import nextflow.SysEnv
import nextflow.cli.CmdAuth
import nextflow.cli.ColorUtil
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException
import nextflow.platform.PlatformHelper

import java.awt.Desktop
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardOpenOption
import java.time.Duration

@Slf4j
@CompileStatic
class AuthCommandImpl implements CmdAuth.AuthCommand {

    static final int API_TIMEOUT_MS = 10000
    static final int AUTH_POLL_TIMEOUT_RETRIES = 60
    static final int AUTH_POLL_INTERVAL_SECONDS = 5
    static final int WORKSPACE_SELECTION_THRESHOLD = 8  // Max workspaces to show in single list; above this uses org-first selection

    /**
     * Creates an HxClient instance with optional authentication token
     */
    protected HxClient createHttpClient(String authToken = null) {
        final builder = HxClient.newBuilder()
            .connectTimeout(Duration.ofMillis(API_TIMEOUT_MS))

        if (authToken) {
            builder.bearerToken(authToken)
        }

        return builder.build()
    }

    @Override
    void login(String apiUrl) {

        // Check if TOWER_ACCESS_TOKEN environment variable is set
        final envToken = SysEnv.get('TOWER_ACCESS_TOKEN')
        if( envToken ) {
            println ""
            ColorUtil.printColored("WARNING: Authentication token is already configured via TOWER_ACCESS_TOKEN environment variable.", "yellow bold")
            ColorUtil.printColored("${ColorUtil.colorize('nextflow auth login', 'cyan')} sets credentials using Nextflow config files, which take precedence over the environment variable.", "dim")
            ColorUtil.printColored(" however, caution is advised to avoid confusing behaviour.", "dim")
            println ""
        }

        // Check if seqera-auth.config file already exists
        final authFile = getAuthFile()
        if( Files.exists(authFile) ) {
            ColorUtil.printColored("Error: Authentication token is already configured in Nextflow config.", "red")
            ColorUtil.printColored("Auth file: ${ColorUtil.colorize(authFile.toString(), 'magenta')}", "dim")
            println " Run ${ColorUtil.colorize('nextflow auth logout', 'cyan')} to remove the current authentication."
            return
        }

        println("Nextflow authentication with Seqera Platform")
        ColorUtil.printColored(" - Authentication will be saved to: ${ColorUtil.colorize(getAuthFile().toString(), 'magenta')}", "dim")

        apiUrl = normalizeApiUrl(apiUrl)
        ColorUtil.printColored(" - Seqera Platform API endpoint: ${ColorUtil.colorize(apiUrl, 'magenta')} (can be customised with ${ColorUtil.colorize('-url', 'cyan')})", "dim")

        // Check if this is a cloud endpoint or enterprise
        final endpointInfo = getCloudEndpointInfo(apiUrl)
        if( endpointInfo.isCloud ) {
            try {
                performAuth0Login(endpointInfo.endpoint as String, endpointInfo.auth as Map)
            } catch( Exception e ) {
                log.debug("Authentication failed", e)
                throw new AbortOperationException("${e.message}")
            }
        } else {
            // Enterprise endpoint - use PAT authentication
            handleEnterpriseAuth(apiUrl)
        }
    }

    protected void performAuth0Login(String apiUrl, Map auth0Config) {

        // Start device authorization flow
        final deviceAuth = requestDeviceAuthorization(auth0Config)

        println ""
        println "Confirmation code: ${ColorUtil.colorize(deviceAuth.user_code as String, 'yellow')}"
        final urlWithCode = "${deviceAuth.verification_uri}?user_code=${deviceAuth.user_code}"
        println "Authentication URL: ${ColorUtil.colorize(urlWithCode, 'magenta')}"
        ColorUtil.printColored("\n[ Press Enter to open in browser ]", "bold")

        // Wait for Enter key with proper interrupt handling
        try {
            final input = System.in.read()
            if( input == -1 ) {
                throw new AbortOperationException("Authentication cancelled")
            }
        } catch( InterruptedException e ) {
            Thread.currentThread().interrupt()
            throw new AbortOperationException("Authentication cancelled")
        } catch( IOException e ) {
            throw new AbortOperationException("Failed to read input: ${e.message}")
        }

        // Try to open browser automatically
        boolean browserOpened = openBrowser(urlWithCode)

        if( !browserOpened ) {
            ColorUtil.printColored("Could not open browser automatically. Please copy the URL above and open it manually in your browser.", "yellow")
        }
        print("${ColorUtil.colorize('Waiting for authentication...', 'dim', true)}")

        try {
            // Poll for device token
            final tokenData = pollForDeviceToken(deviceAuth.device_code as String, deviceAuth.interval as Integer ?: AUTH_POLL_INTERVAL_SECONDS, auth0Config)
            final accessToken = tokenData['access_token'] as String

            // Verify login by calling /user-info
            final userInfo = callUserInfoApi(accessToken, apiUrl)
            println "\n\n${ColorUtil.colorize('✔', 'green', true)} Authentication successful"

            // Generate PAT
            final pat = generatePAT(accessToken, apiUrl)

            // Save to config
            saveAuthToConfig(pat, apiUrl)

            // Automatically run configuration
            try {
                config(false)
            } catch( Exception e ) {
                ColorUtil.printColored("Configuration setup failed: ${e.message}", "red")
                ColorUtil.printColored("You can run 'nextflow auth config' later to set up your configuration.", "dim")
            }

        } catch( InterruptedException e ) {
            Thread.currentThread().interrupt()
            throw new AbortOperationException("Authentication cancelled")
        } catch( Exception e ) {
            throw new RuntimeException("Authentication failed: ${e.message}", e)
        }
    }

    private boolean openBrowser(GString urlWithCode) {
        def browserOpened = false
        try {
            // Method 1: Java Desktop API
            if( Desktop.isDesktopSupported() ) {
                final desktop = Desktop.getDesktop()
                if( desktop.isSupported(Desktop.Action.BROWSE) ) {
                    desktop.browse(new URI(urlWithCode))
                    browserOpened = true
                }
            }

            // Method 2: Platform-specific commands
            if( !browserOpened ) {
                browserOpened = runBrowserCommand(urlWithCode)
            }
        } catch( Exception e ) {
            log.debug("Exception opening browser", e)
        }
        return browserOpened
    }

    protected boolean runBrowserCommand(GString urlWithCode) {
        def command = []
        def browserOpened = false
        final os = System.getProperty("os.name").toLowerCase()
        if( os.contains("mac") || os.contains("darwin") ) {
            command = ["open", urlWithCode]
        } else if( os.contains("win") ) {
            command = ["cmd", "/c", "start", urlWithCode]
        } else {
            // Linux and other Unix-like systems
            final browsers = ["xdg-open", "firefox", "google-chrome", "chromium", "safari"]
            for( browser in browsers ) {
                try {
                    new ProcessBuilder(browser, urlWithCode).start()
                    browserOpened = true
                    break
                } catch( Exception e ) {
                    log.debug("Failed to open browser ${browser}: ${e.message}")
                    // Try next browser
                }
            }
        }

        if( !browserOpened && command ) {
            new ProcessBuilder(command as String[]).start()
            browserOpened = true
        }
        return browserOpened
    }

    private Map requestDeviceAuthorization(Map auth0Config) {
        final params = [
            'client_id': auth0Config.clientId,
            'scope'    : 'openid profile email offline_access',
            'audience' : 'platform'
        ]
        return performAuth0Request("https://${auth0Config.domain}/oauth/device/code", params)
    }

    private Map pollForDeviceToken(String deviceCode, int intervalSeconds, Map auth0Config) {
        final tokenUrl = "https://${auth0Config.domain}/oauth/token"
        def retryCount = 0

        while( retryCount < AUTH_POLL_TIMEOUT_RETRIES ) {
            final params = [
                'grant_type' : 'urn:ietf:params:oauth:grant-type:device_code',
                'device_code': deviceCode,
                'client_id'  : auth0Config.clientId
            ]

            try {
                final result = performAuth0Request(tokenUrl, params)
                return result
            } catch( RuntimeException e ) {
                final message = e.message
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


    protected void handleEnterpriseAuth(String apiUrl) {
        println ""
        ColorUtil.printColored("Please generate a Personal Access Token from your Seqera Platform instance.", "cyan bold")
        println "You can create one at: ${ColorUtil.colorize(getWebUrlFromApiEndpoint(apiUrl) + '/tokens', 'magenta')}"
        println ""

        final pat = promptPAT()

        if( !pat ) {
            throw new AbortOperationException("Personal Access Token is required for Seqera Platform Enterprise authentication")
        }

        // Save to config
        saveAuthToConfig(pat, apiUrl)
        ColorUtil.printColored("Personal Access Token saved to Nextflow auth config (${getAuthFile().toString()})", "green")
    }

    private String getWebUrlFromApiEndpoint(String apiEndpoint) {
        /* Convert API endpoint to web URL
         * e.g., https://api.cloud.seqera.io -> https://cloud.seqera.io
         *      https://cloud.seqera.io/api -> https://cloud.seqera.io
         */
        return apiEndpoint.replace('://api.', '://').replace('/api', '')
    }

    private String promptPAT(){
        System.out.print("Enter your Personal Access Token: ")
        System.out.flush()

        final console = System.console()
        final pat = console ?
            new String(console.readPassword()) :
            new BufferedReader(new InputStreamReader(System.in)).readLine()
        return pat.trim()
    }

    private String generatePAT(String accessToken, String apiUrl) {
        final tokensUrl = "${apiUrl}/tokens"
        final username = System.getProperty("user.name")
        final timestamp = new Date().format("yyyy-MM-dd-HH-mm")
        final tokenName = "nextflow-auth-${username}-${timestamp}"

        final requestBody = new JsonBuilder([name: tokenName]).toString()

        final client = createHttpClient(accessToken)
        final request = HttpRequest.newBuilder()
            .uri(URI.create(tokensUrl))
            .header('Content-Type', 'application/json')
            .POST(HttpRequest.BodyPublishers.ofString(requestBody))
            .build()

        final response = client.send(request, HttpResponse.BodyHandlers.ofString())

        if( response.statusCode() != 200 ) {
            final error = response.body() ?: "HTTP ${response.statusCode()}"
            throw new RuntimeException("Failed to generate PAT: ${error}")
        }

        final json = new JsonSlurper().parseText(response.body()) as Map
        return json.accessKey as String
    }

    private String normalizeApiUrl(String url) {
        if( !url ) {
            // Read config to get the actual resolved endpoint value
            final builder = new ConfigBuilder().setHomeDir(Const.APP_HOME_DIR).setCurrentDir(Const.APP_HOME_DIR)
            final configObject = builder.buildConfigObject()
            final towerConfig = configObject.navigate('tower') as Map ?: [:]
            return PlatformHelper.getEndpoint(towerConfig, SysEnv.get())
        }
        if( !url.startsWith('http://') && !url.startsWith('https://') ) {
            return 'https://' + url
        }
        return url
    }

    protected Map performAuth0Request(String url, Map params) {
        final postData = params.collect { k, v ->
            "${URLEncoder.encode(k.toString(), 'UTF-8')}=${URLEncoder.encode(v.toString(), 'UTF-8')}"
        }.join('&')

        final client = createHttpClient()
        final request = HttpRequest.newBuilder()
            .uri(URI.create(url))
            .header('Content-Type', 'application/x-www-form-urlencoded')
            .POST(HttpRequest.BodyPublishers.ofString(postData))
            .build()

        final response = client.send(request, HttpResponse.BodyHandlers.ofString())

        if( response.statusCode() == 200 ) {
            final json = new JsonSlurper().parseText(response.body())
            return json as Map
        } else {
            final errorBody = response.body()
            if( errorBody ) {
                final errorJson = new JsonSlurper().parseText(errorBody) as Map
                final error = errorJson.error
                throw new RuntimeException("${error}: ${errorJson.error_description ?: ''}")
            } else {
                throw new RuntimeException("Request failed: HTTP ${response.statusCode()}")
            }
        }
    }

    private void saveAuthToConfig(String accessToken, String apiUrl) {
        final config = [:]
        config['tower.accessToken'] = accessToken
        config['tower.endpoint'] = apiUrl
        config['tower.enabled'] = true

        writeConfig(config, null)
    }

    @Override
    void logout() {
        // Check if seqera-auth.config file exists
        final authFile = getAuthFile()
        if( !Files.exists(authFile) ) {
            println "No previous login found.\n"
            return
        }
        // Read token from seqera-auth.config file
        final authConfig = readAuthFile()
        final existingToken = authConfig['tower.accessToken']
        // Extract tower config for PlatformHelper (strip 'tower.' prefix)
        final towerConfig = authConfig.findAll { it.key.toString().startsWith('tower.') }
            .collectEntries { k, v -> [(k.toString().substring(6)): v] }
        final apiUrl = PlatformHelper.getEndpoint(towerConfig, SysEnv.get())

        if( !existingToken ) {
            ColorUtil.printColored("WARN: No authentication token found in auth file.", "yellow bold")
            println "Removing file: ${ColorUtil.colorize(authFile.toString(), 'magenta')}"
            removeAuthFromConfig()
            return
        }

        // Check if TOWER_ACCESS_TOKEN environment variable is set
        final envToken = SysEnv.get('TOWER_ACCESS_TOKEN')
        if( envToken ) {
            println ""
            ColorUtil.printColored("WARNING: TOWER_ACCESS_TOKEN environment variable is set.", "yellow bold")
            println " ${ColorUtil.colorize('nextflow auth logout', 'dim cyan')}${ColorUtil.colorize(' only removes credentials from Nextflow config files.', 'dim')}"
            ColorUtil.printColored(" The environment variable will remain unaffected.", "dim")
            println ""
        }

        ColorUtil.printColored(" - Found authentication token in auth file: ${ColorUtil.colorize(authFile.toString(), 'magenta')}", "dim")
        ColorUtil.printColored(" - Using Seqera Platform endpoint: ${ColorUtil.colorize(apiUrl, 'magenta')}", "dim")

        // Validate token by calling /user-info API
        try {
            final userInfo = callUserInfoApi(existingToken as String, apiUrl)
            ColorUtil.printColored(" - Token is valid for user: ${ColorUtil.colorize(userInfo.userName as String, 'bold')}", "dim")
        } catch( Exception e ) {
            ColorUtil.printColored("Failed to validate token: ${e.message}", "red")
        }

        // Check if we need to delete from platform
        final shouldDeleteFromPlatform = isCloudEndpoint(apiUrl)

        ColorUtil.printColored("\nRunning this command will:", "yellow bold")
        ColorUtil.printColored("  • Remove local Nextflow configuration: ${ColorUtil.colorize(authFile.toString(), 'magenta')}", "yellow")
        if( shouldDeleteFromPlatform ) {
            ColorUtil.printColored("  • Delete the corresponding access token from Seqera Platform: ${ColorUtil.colorize(apiUrl, 'magenta')}", "yellow")
        } else {
            println ""
            ColorUtil.printColored("Warning: Access token not deleted, as using enterprise installation: ${ColorUtil.colorize(apiUrl, 'magenta')}", "yellow")
        }

        final confirmed = promptForYesNo("\n${ColorUtil.colorize('Continue with logout?', 'bold')} (${ColorUtil.colorize('Y', 'green')}/n): ", true)

        if( !confirmed ) {
            println("Logout cancelled.")
            return
        }

        // Only delete PAT from platform if this is a cloud endpoint
        if( shouldDeleteFromPlatform ) {
            try {
                final tokenId = decodeTokenId(existingToken as String)
                deleteTokenViaApi(existingToken as String, apiUrl, tokenId)
            } catch( Exception e ) {
                ColorUtil.printColored("Error removing token: ${e.message}", "red")
            }
        }

        removeAuthFromConfig()
    }

    private String decodeTokenId(String token) {
        try {
            // Decode base64 token
            final decoded = new String(Base64.decoder.decode(token), "UTF-8")

            // Parse JSON to extract token ID
            final json = new JsonSlurper().parseText(decoded) as Map
            final tokenId = json.tid

            if( !tokenId ) {
                throw new RuntimeException("No token ID found in decoded token")
            }

            return tokenId.toString()
        } catch( Exception e ) {
            throw new RuntimeException("Failed to decode token ID: ${e.message}")
        }
    }

    private void deleteTokenViaApi(String token, String apiUrl, String tokenId) {
        final client = createHttpClient(token)
        final request = HttpRequest.newBuilder()
            .uri(URI.create("${apiUrl}/tokens/${tokenId}"))
            .DELETE()
            .build()

        final response = client.send(request, HttpResponse.BodyHandlers.ofString())

        if( response.statusCode() != 200 && response.statusCode() != 204 ) {
            final error = response.body() ?: "HTTP ${response.statusCode()}"
            throw new RuntimeException("Failed to delete token: ${error}")
        }

        println "\n${ColorUtil.colorize('✔', 'green', true)} Token successfully deleted from Seqera Platform."
    }

    private void removeAuthFromConfig() {
        final configFile = getConfigFile()
        final authFile = getAuthFile()

        // Remove includeConfig line from main config file
        if( Files.exists(configFile) ) {
            final existingContent = Files.readString(configFile)
            final updatedContent = removeIncludeConfigLine(existingContent)
            Files.writeString(configFile, updatedContent.toString(), StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)
        }

        // Delete seqera-auth.config file
        if( Files.exists(authFile) ) {
            Files.delete(authFile)
        }

        println "${ColorUtil.colorize('✔', 'green', true)} Authentication removed from Nextflow config."
    }

    @Override
    void config(Boolean showHeader = true) {
        // Read from both main config and seqera-auth.config file
        final builder = new ConfigBuilder().setHomeDir(Const.APP_HOME_DIR).setCurrentDir(Const.APP_HOME_DIR)
        final configObject = builder.buildConfigObject()
        final config = configObject.flatten()

        // Navigate to tower config section (returns map without 'tower.' prefix)
        final towerConfig = configObject.navigate('tower') as Map ?: [:]
        final existingToken = PlatformHelper.getAccessToken(towerConfig, SysEnv.get())
        final endpoint = PlatformHelper.getEndpoint(towerConfig, SysEnv.get())

        if( !existingToken ) {
            println "No authentication found. Please run ${ColorUtil.colorize('nextflow auth login', 'cyan')} first."
            return
        }

        if (showHeader) {
            println "Nextflow Seqera Platform configuration"
            ColorUtil.printColored(" - Config file: ${ColorUtil.colorize(getAuthFile().toString(), 'magenta')}", "dim")

            // Check if token is from environment variable
            if( !towerConfig['accessToken'] && SysEnv.get('TOWER_ACCESS_TOKEN') ) {
                ColorUtil.printColored(" - Using access token from TOWER_ACCESS_TOKEN environment variable", "dim")
            }
        }

        try {
            // Get user info to validate token and get user ID
            final userInfo = callUserInfoApi(existingToken as String, endpoint as String)
            ColorUtil.printColored(" - Authenticated as: ${ColorUtil.colorize(userInfo.userName as String, 'cyan bold')}", "dim")
            println ""

            // Track if any changes are made
            def configChanged = false

            // Configure tower.enabled
            configChanged |= configureEnabled(config)

            // Configure workspace
            final workspaceResult = configureWorkspace(config, existingToken as String, endpoint as String, userInfo.id as String)
            configChanged = configChanged || (workspaceResult.changed as boolean)

            // Save updated config only if changes were made
            if( configChanged ) {
                writeConfig(config, workspaceResult.metadata as Map)
                println "\n${ColorUtil.colorize('✔', 'green', true)} Configuration saved to ${ColorUtil.colorize(getAuthFile().toString(), 'magenta')}"
            }

            // Configure compute environment for the workspace (always run after workspace selection)
            final currentWorkspaceId = config.get('tower.workspaceId') as String
            if( currentWorkspaceId ) {
                // Get workspace metadata if not already available (e.g., when user kept existing workspace)
                def workspaceMetadata = workspaceResult.metadata as Map
                if( !workspaceMetadata ) {
                    workspaceMetadata = getWorkspaceDetailsFromApi(existingToken as String, endpoint as String, currentWorkspaceId)
                }
                configureComputeEnvironment(existingToken as String, endpoint as String, currentWorkspaceId, workspaceMetadata)
            }

        } catch( Exception e ) {
            throw new AbortOperationException("Failed to configure settings: ${e.message}")
        }
    }

    private boolean configureEnabled(Map config) {
        final currentEnabled = config.get('tower.enabled', false)

        println "Workflow monitoring settings. Current setting: ${currentEnabled ? ColorUtil.colorize('enabled', 'green') : ColorUtil.colorize('disabled', 'red')}"
        ColorUtil.printColored("  When enabled, all workflow runs are automatically monitored by Seqera Platform", "dim")
        ColorUtil.printColored("  When disabled, you can enable per-run with the ${ColorUtil.colorize('-with-tower', 'cyan')} flag", "dim")
        println ""

        final promptText = "${ColorUtil.colorize('Enable workflow monitoring for all runs?', 'bold', true)} (${currentEnabled ? ColorUtil.colorize('Y', 'green') + '/n' : 'y/' + ColorUtil.colorize('N', 'red')}): "
        final input = promptForYesNo(promptText, currentEnabled)

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
        final envWorkspaceId = SysEnv.get('TOWER_WORKFLOW_ID')
        if( envWorkspaceId ) {
            println "\nDefault workspace: ${ColorUtil.colorize('TOWER_WORKFLOW_ID environment variable is set', 'yellow')}"
            ColorUtil.printColored("  Not prompting for default workspace configuration as environment variable takes precedence", "dim")
            return [changed: false, metadata: null]
        }

        // Get all workspaces for the user
        final workspaces = getUserWorkspaces(accessToken, endpoint, userId)

        if( !workspaces ) {
            println "\nNo workspaces found for your account."
            return [changed: false, metadata: null]
        }

        // Show current workspace setting
        final currentWorkspaceId = config.get('tower.workspaceId')

        String currentSetting = getCurrentWorkspaceName(workspaces, config.get('tower.workspaceId'))

        println "\nDefault workspace. Current setting: ${ColorUtil.colorize(currentSetting, 'cyan', true)}"
        ColorUtil.printColored("  Workflow runs use this workspace by default", "dim")
        // Group by organization
        final orgWorkspaces = workspaces.groupBy { ((Map) it).orgName ?: 'Personal' }

        // If threshold or fewer total options, show all at once
        if( workspaces.size() <= WORKSPACE_SELECTION_THRESHOLD ) {
            return selectWorkspaceFromAll(config, workspaces, currentWorkspaceId)
        } else {
            // Two-stage selection: org first, then workspace
            return selectWorkspaceByOrg(config, orgWorkspaces, currentWorkspaceId)
        }
    }

    private String getCurrentWorkspaceName(List<Object> workspaces, currentWorkspaceId) {
        final currentWorkspace = workspaces.find { ((Map) it).workspaceId.toString() == currentWorkspaceId?.toString() } as Map
        return currentWorkspace ? "${currentWorkspace.orgName} / ${currentWorkspace.workspaceName}" : "None (Personal workspace)"
    }

    private Map selectWorkspaceFromAll(Map config, List workspaces, final currentWorkspaceId) {
        println "\nAvailable workspaces:"
        final isPersonalWorkspace = !currentWorkspaceId
        final currentIndicator = isPersonalWorkspace ? ColorUtil.colorize(' (current)', 'yellow bold') : ''
        println "  0. ${ColorUtil.colorize('None (Personal workspace)', 'cyan', true)} ${ColorUtil.colorize('[no organization]', 'dim', true)}${currentIndicator}"

        // Sort workspaces by org name, then workspace name
        final sortedWorkspaces = workspaces.sort { a, b ->
            final aMap = a as Map
            final bMap = b as Map
            final orgCompare = (aMap.orgName as String ?: '').compareToIgnoreCase(bMap.orgName as String ?: '')
            orgCompare != 0 ? orgCompare : (aMap.workspaceName as String ?: '').compareToIgnoreCase(bMap.workspaceName as String ?: '')
        }

        sortedWorkspaces.eachWithIndex { workspace, index ->
            final ws = workspace as Map
            final isCurrent = ws.workspaceId.toString() == currentWorkspaceId?.toString()
            final currentInd = isCurrent ? ColorUtil.colorize(' (current)', 'yellow bold') : ''
            final prefix = ws.orgName ? "${ColorUtil.colorize(ws.orgName as String, 'cyan', true)} / " : ""
            println "  ${index + 1}. ${prefix}${ColorUtil.colorize(ws.workspaceName as String, 'magenta', true)} ${ColorUtil.colorize('[' + (ws.workspaceFullName as String) + ']', 'dim', true)}${currentInd}"
        }

        // Show current workspace and prepare prompt
        final currentWorkspaceName = getCurrentWorkspaceName(sortedWorkspaces, currentWorkspaceId)

        println("\n${ColorUtil.colorize('Leave blank to keep current setting', 'bold')} (${ColorUtil.colorize(currentWorkspaceName, 'cyan')}),")
        final selection = promptForNumber(ColorUtil.colorize("or select workspace (0-${sortedWorkspaces.size()}): ", 'bold', true), 0, sortedWorkspaces.size(), true)

        if( selection == null ) {
            return [changed: false, metadata: null]
        }

        if( selection == 0 ) {
            final hadWorkspaceId = config.containsKey('tower.workspaceId')
            config.remove('tower.workspaceId')
            return [changed: hadWorkspaceId, metadata: null]
        } else {
            final selectedWorkspace = sortedWorkspaces[selection - 1] as Map
            final selectedId = selectedWorkspace.workspaceId.toString()
            final currentId = config.get('tower.workspaceId')
            config['tower.workspaceId'] = selectedId
            final metadata = [
                orgName          : selectedWorkspace.orgName,
                workspaceName    : selectedWorkspace.workspaceName,
                workspaceFullName: selectedWorkspace.workspaceFullName
            ]
            return [changed: currentId != selectedId, metadata: metadata]
        }
    }

    private Map selectWorkspaceByOrg(Map config, Map orgWorkspaces, final currentWorkspaceId) {
        // Get current workspace info for prompts
        final allWorkspaces = []
        orgWorkspaces.values().each { workspaceList ->
            allWorkspaces.addAll(workspaceList as List)
        }
        final currentWorkspaceDisplay = getCurrentWorkspaceName(allWorkspaces, currentWorkspaceId)

        // First, select organization
        final orgs = orgWorkspaces.keySet().toList().sort { (it as String).toLowerCase() }

        // Always add Personal as first option (it's never returned by the API but should always be available)
        orgs.add(0, 'Personal')

        println "\nAvailable organizations:"
        orgs.eachWithIndex { orgName, index ->
            final displayName = orgName == 'Personal' ? 'None [Personal workspace]' : orgName
            println "  ${index + 1}. ${ColorUtil.colorize(displayName as String, 'cyan', true)}"
        }
        println("\n${ColorUtil.colorize('Leave blank to keep current setting', 'bold')} (${ColorUtil.colorize(currentWorkspaceDisplay, 'cyan')}),")
        final orgSelection = promptForNumber(ColorUtil.colorize("or select organization (1-${orgs.size()}): ", 'bold', true), 1, orgs.size(),true)
        if (!orgSelection)
            return [changed: false, metadata: null]

        final selectedOrgName = orgs[orgSelection - 1]

        // If Personal was selected, remove workspace ID (use personal workspace)
        if( selectedOrgName == 'Personal' ) {
            final hadWorkspaceId = config.containsKey('tower.workspaceId')
            config.remove('tower.workspaceId')
            return [changed: hadWorkspaceId, metadata: null]
        }

        final orgWorkspaceList = (orgWorkspaces[selectedOrgName] as List).sort { a, b ->
            final aMap = a as Map
            final bMap = b as Map
            (aMap.workspaceName as String ?: '').compareToIgnoreCase(bMap.workspaceName as String ?: '')
        }

        println ""
        println "Select workspace in ${selectedOrgName}:"

        orgWorkspaceList.eachWithIndex { workspace, index ->
            final ws = workspace as Map
            final isCurrent = ws.workspaceId.toString() == currentWorkspaceId?.toString()
            final currentInd = isCurrent ? ColorUtil.colorize(' (current)', 'yellow bold') : ''
            println "  ${index + 1}. ${ColorUtil.colorize(ws.workspaceName as String, 'magenta', true)} ${ColorUtil.colorize('[' + (ws.workspaceFullName as String) + ']', 'dim', true)}${currentInd}"
        }

        final maxSelection = orgWorkspaceList.size()
        final wsSelection = promptForNumber(ColorUtil.colorize("\nSelect workspace (1-${maxSelection}): ", 'bold', true), 1, maxSelection,false)

        final selectedWorkspace = orgWorkspaceList[wsSelection - 1] as Map
        final selectedId = selectedWorkspace.workspaceId.toString()
        final currentId = config.get('tower.workspaceId')
        config['tower.workspaceId'] = selectedId
        final metadata = [
            orgName          : selectedWorkspace.orgName,
            workspaceName    : selectedWorkspace.workspaceName,
            workspaceFullName: selectedWorkspace.workspaceFullName
        ]
        return [changed: currentId != selectedId, metadata: metadata]
    }

    private List getUserWorkspaces(String accessToken, String endpoint, String userId) {
        final client = createHttpClient(accessToken)
        final request = HttpRequest.newBuilder()
            .uri(URI.create("${endpoint}/user/${userId}/workspaces"))
            .GET()
            .build()

        final response = client.send(request, HttpResponse.BodyHandlers.ofString())

        if( response.statusCode() != 200 ) {
            final error = response.body() ?: "HTTP ${response.statusCode()}"
            throw new RuntimeException("Failed to get workspaces: ${error}")
        }

        final json = new JsonSlurper().parseText(response.body()) as Map
        final orgsAndWorkspaces = json.orgsAndWorkspaces as List

        return orgsAndWorkspaces.findAll { ((Map) it).workspaceId != null }
    }

    private List getComputeEnvironments(String accessToken, String endpoint, String workspaceId) {
        final client = createHttpClient(accessToken)
        final uri = workspaceId ?
            "${endpoint}/compute-envs?workspaceId=${workspaceId}" :
            "${endpoint}/compute-envs"

        final request = HttpRequest.newBuilder()
            .uri(URI.create(uri))
            .GET()
            .build()

        final response = client.send(request, HttpResponse.BodyHandlers.ofString())

        if( response.statusCode() != 200 ) {
            final error = response.body() ?: "HTTP ${response.statusCode()}"
            throw new RuntimeException("Failed to get compute environments: ${error}")
        }

        final json = new JsonSlurper().parseText(response.body()) as Map
        return json.computeEnvs as List ?: []
    }

    private void setPrimaryComputeEnvironment(String accessToken, String endpoint, String computeEnvId, String workspaceId) {
        final client = createHttpClient(accessToken)
        final uri = "${endpoint}/compute-envs/${computeEnvId}/primary?workspaceId=${workspaceId}"

        final requestBody = new JsonBuilder([:]).toString()

        final request = HttpRequest.newBuilder()
            .uri(URI.create(uri))
            .header('Content-Type', 'application/json')
            .POST(HttpRequest.BodyPublishers.ofString(requestBody))
            .build()

        final response = client.send(request, HttpResponse.BodyHandlers.ofString())

        if( response.statusCode() != 200 && response.statusCode() != 204 ) {
            final error = response.body() ?: "HTTP ${response.statusCode()}"
            throw new RuntimeException("Failed to set primary compute environment: ${error}")
        }
    }

    private void configureComputeEnvironment(String accessToken, String endpoint, String workspaceId, Map workspaceMetadata) {
        try {
            // Get compute environments for the workspace
            final computeEnvs = getComputeEnvironments(accessToken, endpoint, workspaceId)

            // If there are zero compute environments, log a warning and provide a link
            if( computeEnvs.isEmpty() ) {
                println ""
                ColorUtil.printColored("Warning: No compute environments found in this workspace.", "red")
                ColorUtil.printColored("  You must create a compute environment to run pipelines.", "yellow")

                // Generate the web URL for creating a compute environment
                if( workspaceMetadata ) {
                    final orgName = workspaceMetadata.orgName as String
                    final workspaceName = workspaceMetadata.workspaceName as String
                    final webUrl = getWebUrlFromApiEndpoint(endpoint)
                    final createUrl = "${webUrl}/orgs/${orgName}/workspaces/${workspaceName}/compute-envs/new"
                    println "  Create one at: ${ColorUtil.colorize(createUrl, 'magenta')}"
                }
                return
            }

            // Find which compute environment is already set as primary
            final primaryEnv = computeEnvs.find { ((Map) it).primary == true }
            final currentPrimaryName = primaryEnv ? (primaryEnv as Map).name as String : 'None'

            // Show current setting
            println ""
            println "Primary compute environment. Current setting: ${ColorUtil.colorize(currentPrimaryName, 'cyan', true)}"
            ColorUtil.printColored("  Setting a primary compute environment allows you to run pipelines without", "dim")
            ColorUtil.printColored("  explicitly specifying a compute environment every time.", "dim")
            println ""
            println "Available compute environments:"

            // Sort compute environments by name
            final sortedComputeEnvs = computeEnvs.sort { a, b ->
                final aMap = a as Map
                final bMap = b as Map
                (aMap.name as String ?: '').compareToIgnoreCase(bMap.name as String ?: '')
            }

            sortedComputeEnvs.eachWithIndex { ce, index ->
                final env = ce as Map
                final name = env.name as String
                final platform = env.platform as String
                final isPrimary = env.primary == true
                final currentIndicator = isPrimary ? ColorUtil.colorize(' (current)', 'yellow bold') : ''
                println "  ${index + 1}. ${ColorUtil.colorize(name, 'cyan', true)} ${ColorUtil.colorize('[' + platform + ']', 'dim yellow', true)}${currentIndicator}"
            }

            println ""
            println "${ColorUtil.colorize('Leave blank to keep current setting', 'bold')} (${ColorUtil.colorize(currentPrimaryName, 'cyan')}),".toString()
            final selection = promptForNumber(ColorUtil.colorize("or select compute environment (1-${sortedComputeEnvs.size()}): ", 'bold', true), 1, sortedComputeEnvs.size(), true)

            if( selection == null ) {
                return
            }

            final selectedEnv = sortedComputeEnvs[selection - 1] as Map
            final computeEnvId = selectedEnv.id as String

            // Only set if different from current primary
            if( primaryEnv && (primaryEnv as Map).id == computeEnvId ) {
                // User selected the same one that's already primary
                return
            }

            // Set the selected compute environment as primary
            setPrimaryComputeEnvironment(accessToken, endpoint, computeEnvId, workspaceId)
            ColorUtil.printColored("Primary compute environment set to: ${ColorUtil.colorize(selectedEnv.name as String, 'cyan')}", "green")

        } catch( Exception e ) {
            ColorUtil.printColored("Warning: Failed to configure compute environment: ${e.message}", "yellow")
        }
    }

    private Boolean promptForYesNo(String prompt, Boolean defaultValue) {
        while( true ) {
            final input = readUserInput(prompt)?.toLowerCase()

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
            final input = readUserInput(prompt)

            if( input?.isEmpty() && allowEmpty ) {
                return null
            }

            try {
                final number = Integer.parseInt(input)
                if( number >= min && number <= max ) {
                    return number
                }
            } catch( NumberFormatException ignored ) {
                // Fall through to error message
            }
            ColorUtil.printColored("Invalid input. Please enter a number between ${min} and ${max}.", "red")
        }
    }

    @Override
    void status() {
        final config = readConfig()
        printStatus(collectStatus(config))
    }

    private ConfigStatus collectStatus(Map config) {
        // Collect all status information
        final status = new ConfigStatus([], null, null)

        // Extract tower config and strip prefix for PlatformHelper
        final towerConfig = config.findAll { it.key.toString().startsWith('tower.') }
            .collectEntries { k, v -> [(k.toString().substring(6)): v] }

        // API endpoint - use PlatformHelper
        final String endpoint = PlatformHelper.getEndpoint(towerConfig, SysEnv.get())
        final endpointInfo = getConfigValue(config, 'tower.endpoint', 'TOWER_API_ENDPOINT')
        status.table.add(['API endpoint', ColorUtil.colorize(endpoint, 'magenta'), (endpointInfo.source ?: 'default') as String])

        // API connection check
        final apiConnectionOk = checkApiConnection(endpoint)
        final connectionColor = apiConnectionOk ? 'green' : 'red'
        status.table.add(['API connection', ColorUtil.colorize(apiConnectionOk ? '✔ OK' : 'ERROR', connectionColor), ''])

        // Authentication check - use PlatformHelper
        final String accessToken = PlatformHelper.getAccessToken(towerConfig, SysEnv.get())

        // Determine source for display
        final tokenInfo = getConfigValue(config, 'tower.accessToken', 'TOWER_ACCESS_TOKEN')
        final String tokenSource = tokenInfo.source ?: 'not set'

        if( accessToken ) {
            try {
                final userInfo = callUserInfoApi(accessToken, endpoint)
                final currentUser = userInfo.userName as String
                status.table.add(['Authentication', "${ColorUtil.colorize('✔ OK', 'green')} (user: ${ColorUtil.colorize(currentUser, 'cyan')})".toString(), tokenSource])
            } catch( Exception e ) {
                status.table.add(['Authentication', ColorUtil.colorize('✘ Connection check failed', 'red'), tokenSource])
            }
        } else {
            status.table.add(['Authentication', ColorUtil.colorize('Not set', 'red'), 'not set'])
        }

        // Monitoring enabled
        final enabledInfo = getConfigValue(config, 'tower.enabled', null)
        final enabledValue = enabledInfo.value?.toString()?.toLowerCase() in ['true', '1', 'yes'] ? 'Yes' : 'No'
        final enabledColor = enabledValue == 'Yes' ? 'green' : 'yellow'
        status.table.add(['Workflow monitoring', ColorUtil.colorize(enabledValue, enabledColor), (enabledInfo.source ?: 'default') as String])

        // Default workspace - use PlatformHelper
        final String workspaceId = PlatformHelper.getWorkspaceId(towerConfig, SysEnv.get())
        final workspaceInfo = getConfigValue(config, 'tower.workspaceId', 'TOWER_WORKSPACE_ID')
        if( workspaceId ) {
            // Try to get workspace name from API if we have a token
            def workspaceDetails = null
            if( accessToken ) {
                workspaceDetails = getWorkspaceDetailsFromApi(accessToken, endpoint, workspaceId)
            }

            if( workspaceDetails ) {
                // Add workspace ID row and remember its index
                status.workspaceRowIndex = status.table.size()
                status.table.add(['Default workspace', ColorUtil.colorize(workspaceId, 'cyan'), workspaceInfo.source as String])
                // Store workspace details for display after this row (outside table structure)
                status.workspaceInfo = workspaceDetails
            } else {
                status.table.add(['Default workspace', ColorUtil.colorize(workspaceId, 'cyan', true), workspaceInfo.source as String])
            }
        } else {
            if( accessToken ) {
                status.table.add(['Default workspace', ColorUtil.colorize('None (Personal workspace)', 'cyan', true), 'default'])
            } else {
                status.table.add(['Default workspace', ColorUtil.colorize('N/A', 'dim', true), 'default'])
            }
        }

        // Primary compute environment and work directory
        def primaryEnv = null
        if( accessToken && workspaceId ) {
            try {
                final computeEnvs = getComputeEnvironments(accessToken, endpoint, workspaceId)
                primaryEnv = computeEnvs.find { ((Map) it).primary == true } as Map
            } catch( Exception e ) {
                status.table.add(['Primary compute env', ColorUtil.colorize('Error fetching', 'red'), ''])
                status.table.add(['Default work dir', ColorUtil.colorize('N/A', 'dim', true), ''])
                return status
            }
        }

        if( primaryEnv ) {
            final envName = primaryEnv.name as String
            final envPlatform = primaryEnv.platform as String
            final displayValue = "${ColorUtil.colorize(envName, 'cyan')} ${ColorUtil.colorize('[' + envPlatform + ']', 'dim yellow', true)}".toString()
            status.table.add(['Primary compute env', displayValue, 'workspace'])
            status.table.add(['Default work dir', ColorUtil.colorize(primaryEnv.workDir as String, 'magenta'), 'compute env'])
        } else {
            final ceValue = (!accessToken || !workspaceId) ? 'N/A' : 'None'
            final ceColor = ceValue == 'None' ? 'yellow' : 'dim'
            status.table.add(['Primary compute env', ColorUtil.colorize(ceValue, ceColor, true), 'workspace'])
            status.table.add(['Default work dir', ColorUtil.colorize(ceValue, ceColor, true), 'compute env'])
        }

        return status
    }
    private void printStatus(ConfigStatus status){
        println ""

        // Generate all table lines with consistent column widths
        final tableLines = generateStatusTableLines(status.table)

        if( status.workspaceInfo && status.workspaceRowIndex != null ) {
            /* Insert workspace details after the workspace row
             * workspaceRowIndex is the row in the data array (0-indexed)
             * In tableLines: [0]=header, [1]=separator, [2]=first data row, etc.
             * So workspace row is at index: workspaceRowIndex + 2
             * We want to insert AFTER it, so: workspaceRowIndex + 2 + 1 = workspaceRowIndex + 3
             */
            final insertAfterLine = status.workspaceRowIndex + 3

            final workspaceDetails = [
                ColorUtil.colorize("${" " * 22}$status.workspaceInfo.orgName / $status.workspaceInfo.workspaceName", "cyan", true),
                ColorUtil.colorize("${" " * 22}$status.workspaceInfo.workspaceFullName", "dim", true)
            ]

            // Insert workspace details into the output
            tableLines.addAll(insertAfterLine, workspaceDetails)
        }

        // Print all lines
        tableLines.each { println it }
    }

    private List<String> generateStatusTableLines(List<List<String>> rows) {
        if( !rows ) return []

        final List<String> lines = []

        // Calculate column widths (accounting for ANSI codes)
        def col1Width = rows.collect { stripAnsiCodes(it[0]).length() }.max()
        def col2Width = rows.collect { stripAnsiCodes(it[1]).length() }.max()
        def col3Width = rows.collect { stripAnsiCodes(it[2]).length() }.max()

        // Add some padding
        col1Width = Math.max(col1Width, 15) + 2
        col2Width = Math.max(col2Width, 15) + 2
        col3Width = Math.max(col3Width, 10) + 2

        // Add table header
        lines.add(ColorUtil.colorize("${'Setting'.padRight(col1Width)} ${'Value'.padRight(col2Width)} Source", "bold"))
        lines.add("${'-' * col1Width} ${'-' * col2Width} ${'-' * col3Width}".toString())

        // Add rows
        rows.each { row ->
            final paddedCol1 = padStringWithAnsi(row[0], col1Width)
            final paddedCol2 = padStringWithAnsi(row[1], col2Width)
            final paddedCol3 = ColorUtil.colorize(row[2], 'dim', true)
            lines.add("${paddedCol1} ${paddedCol2} ${paddedCol3}".toString())
        }

        return lines
    }

    private String stripAnsiCodes(String text) {
        return text?.replaceAll(/\u001B\[[0-9;]*m/, '') ?: ''
    }

    private String padStringWithAnsi(String text, int width) {
        final plainText = stripAnsiCodes(text)
        final padding = width - plainText.length()
        return padding > 0 ? text + (' ' * padding) : text
    }

    private Map getConfigValue(Map config, String configKey, String envVarName) {
        //Checks where the config value came from
        final configValue = config[configKey]
        final envValue = envVarName ? SysEnv.get(envVarName) : null
        final effectiveValue = configValue ?: envValue

        def source = null
        if( configValue ) {
            source = "nextflow config"
        } else if( envValue ) {
            source = "env var \$${envVarName}"
        }

        return [
            value     : effectiveValue,
            source    : source,
            fromConfig: configValue != null,
            fromEnv   : envValue != null
        ]
    }

    protected boolean checkApiConnection(String endpoint) {
        try {
            final client = createHttpClient()
            final request = HttpRequest.newBuilder()
                .uri(URI.create("${endpoint}/service-info"))
                .GET()
                .build()

            final response = client.send(request, HttpResponse.BodyHandlers.ofString())
            return response.statusCode() == 200
        } catch( Exception e ) {
            log.debug("Failed to connect to API endpoint ${endpoint}: ${e.message}", e)
            return false
        }
    }

    protected static String readUserInput(String message = null) {
        if (message) {
            System.out.print(message)
            System.out.flush()
        }
        final console = System.console()
        final line = console != null
            ? console.readLine()
            : new BufferedReader(new InputStreamReader(System.in)).readLine()
        return line?.trim()
    }

    protected Map getWorkspaceDetailsFromApi(String accessToken, String endpoint, String workspaceId) {
        try {
            final userInfo = callUserInfoApi(accessToken, endpoint)
            final userId = userInfo.id as String

            final client = createHttpClient(accessToken)
            final request = HttpRequest.newBuilder()
                .uri(URI.create("${endpoint}/user/${userId}/workspaces"))
                .GET()
                .build()

            final response = client.send(request, HttpResponse.BodyHandlers.ofString())

            if (response.statusCode() != 200) {
                return null
            }

            final json = new JsonSlurper().parseText(response.body()) as Map
            final orgsAndWorkspaces = json.orgsAndWorkspaces as List

            final workspace = orgsAndWorkspaces.find { ((Map)it).workspaceId?.toString() == workspaceId }
            if (workspace) {
                final ws = workspace as Map
                return [
                    orgName: ws.orgName,
                    workspaceName: ws.workspaceName,
                    workspaceFullName: ws.workspaceFullName
                ]
            }

            return null
        } catch (Exception e) {
            log.debug("Failed to get workspace details for workspace ${workspaceId}: ${e.message}", e)
            return null
        }
    }

    private Map getCloudEndpointInfo(String apiUrl) {
        // Check if this is a standard cloud endpoint
        final authDomain = PlatformHelper.getAuthDomain(apiUrl)
        if (authDomain) {
            final clientId = PlatformHelper.getAuthClientId(apiUrl)
            return [
                isCloud: true,
                endpoint: apiUrl,
                auth: [domain: authDomain, clientId: clientId]
            ]
        }

        // Check for legacy URL format (e.g., https://cloud.seqera.io/api)
        if (apiUrl.contains('://cloud.') && apiUrl.endsWith('/api')) {
            final legacyToStandard = apiUrl.replace('://cloud.', '://api.cloud.').replaceAll('/api$', '')
            final legacyAuthDomain = PlatformHelper.getAuthDomain(legacyToStandard)
            if (legacyAuthDomain) {
                final clientId = PlatformHelper.getAuthClientId(legacyToStandard)
                return [
                    isCloud: true,
                    endpoint: legacyToStandard,
                    auth: [domain: legacyAuthDomain, clientId: clientId]
                ]
            }
        }

        return [isCloud: false, endpoint: apiUrl, auth: null]
    }

    private boolean isCloudEndpoint(String apiUrl) {
        return getCloudEndpointInfo(apiUrl).isCloud
    }

    protected Map callUserInfoApi(String accessToken, String apiUrl) {
        final client = createHttpClient(accessToken)
        final request = HttpRequest.newBuilder()
            .uri(URI.create("${apiUrl}/user-info"))
            .GET()
            .build()

        final response = client.send(request, HttpResponse.BodyHandlers.ofString())

        if (response.statusCode() != 200) {
            final error = response.body() ?: "HTTP ${response.statusCode()}"
            throw new RuntimeException("Failed to get user info: ${error}")
        }

        final json = new JsonSlurper().parseText(response.body()) as Map
        return json.user as Map
    }

    protected Path getConfigFile() {
        return Const.APP_HOME_DIR.resolve('config')
    }

    protected Path getAuthFile() {
        return Const.APP_HOME_DIR.resolve('seqera-auth.config')
    }

    private Map readAuthFile() {
        final configFile = getAuthFile()
        if (!Files.exists(configFile)) {
            return [:]
        }

        try {
            final configText = Files.readString(configFile)
            final config = new ConfigSlurper().parse(configText)
            return config.flatten()
        } catch (Exception e) {
            throw new RuntimeException("Failed to read config file ${configFile}: ${e.message}")
        }
    }

    private Map readConfig() {
        final builder = new ConfigBuilder().setHomeDir(Const.APP_HOME_DIR).setCurrentDir(Const.APP_HOME_DIR)
        return builder.buildConfigObject().flatten()
    }

    private void writeConfig(Map config, Map workspaceMetadata = null) {
        log.debug("Getting config ")
        final configFile = getConfigFile()
        final authFile = getAuthFile()

        // Create directory if it doesn't exist
        if (!Files.exists(configFile.parent)) {
            Files.createDirectories(configFile.parent)
        }

        // Write tower config to seqera-auth.config file
        final towerConfig = config.findAll { key, value ->
            key.toString().startsWith('tower.')
        }

        final authConfigText = new StringBuilder()
        authConfigText.append("// Seqera Platform configuration\n")
        authConfigText.append("tower {\n")
        towerConfig.each { key, value ->
            final configKey = key.toString().substring(6) // Remove "tower." prefix

            if (value instanceof String) {
                def line = "    ${configKey} = '${value}'"
                // Add workspace comment if this is workspaceId and we have metadata
                if (configKey == 'workspaceId' && workspaceMetadata) {
                    line += "  // ${workspaceMetadata.orgName} / ${workspaceMetadata.workspaceName} [${workspaceMetadata.workspaceFullName}]"
                }
                authConfigText.append("${line}\n")
            } else {
                authConfigText.append("    ${configKey} = ${value}\n")
            }
        }
        authConfigText.append("}\n")

        // Write the seqera-auth.config file
        Files.writeString(authFile, authConfigText.toString(), StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)

        // Add includeConfig line to main config file if it doesn't exist
        addIncludeConfigToMainFile(configFile)
    }

    private void addIncludeConfigToMainFile(Path configFile) {
        final includeConfigLine = "includeConfig 'seqera-auth.config'"

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
        // Remove the includeConfig 'seqera-auth.config' line
        final lines = content.split('\n')
        final filteredLines = lines.findAll { line ->
            !line.trim().equals("includeConfig 'seqera-auth.config'")
        }
        return filteredLines.join('\n')
    }
    @Canonical
    static class ConfigStatus {
        List<List<String>> table
        Map workspaceInfo
        Integer workspaceRowIndex  // Track which row has the workspace so we can insert details after it
    }
}


