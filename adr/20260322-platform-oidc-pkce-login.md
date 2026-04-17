# Replace Auth0 with Platform OIDC PKCE for Nextflow CLI login

- Authors: Paolo Di Tommaso
- Status: draft
- Date: 2026-03-22
- Tags: auth, oidc, pkce, nextflow

## Summary

Replace the Auth0 Device Authorization Grant used by `nextflow auth login` with an OAuth2 Authorization Code + PKCE flow against Seqera Platform itself acting as the OIDC identity provider. Implement the reusable OIDC PKCE client as `lib-platform-oidc` in libseqera.

## Problem Statement

The `nextflow auth login` command authenticates via Auth0, an external identity provider, using the OAuth2 Device Authorization Grant. This requires hardcoded Auth0 domain and client ID mappings per environment (dev, stage, prod) and creates an artificial distinction between "cloud" and "enterprise" endpoints — cloud uses Auth0 while enterprise falls back to manual PAT entry.

Seqera Platform now exposes a standards-compliant OIDC provider at `/.well-known/openid-configuration`, making it possible for CLI clients to authenticate directly against Platform using Authorization Code + PKCE — eliminating the Auth0 dependency and unifying the login flow for all Platform instances.

## Goals

- Single login flow for all Platform endpoints (cloud and enterprise)
- Eliminate Auth0 dependency and hardcoded domain/clientId mappings
- Deprecate `TOWER_AUTH_DOMAIN` and `TOWER_AUTH_CLIENT_ID` env vars
- Reusable OIDC PKCE library in libseqera for other CLI tools
- No changes to PAT generation, storage, runtime usage, refresh, or logout

## Non-goals

- Changing how PATs are generated or stored after login
- Modifying runtime token refresh (`TowerXAuth`, `TowerClient`, `WaveClient`)
- Removing backward compatibility with existing PATs in config

Note: The logout flow was simplified as a consequence of removing the cloud-vs-enterprise distinction — PAT deletion via Platform API is now always attempted for all endpoints.

## Login Flow

```
  Nextflow CLI                          Browser                         Platform
       │                                                                    │
       │  GET /.well-known/openid-configuration                             │
       │───────────────────────────────────────────────────────────────────>│
       │<───────────── { authorization_endpoint, token_endpoint } ──────────│
       │                                                                    │
       │  [generate code_verifier, code_challenge, state]                   │
       │  [start local HTTP server on 127.0.0.1:PORT]                       │
       │                                                                    │
       │  open browser ──>│                                                 │
       │                  │  GET /authorize?                                │
       │                  │    client_id=nextflow_cli                       │
       │                  │    &response_type=code                          │
       │                  │    &scope=openid+profile+email+offline_access   │
       │                  │    &redirect_uri=http://127.0.0.1:PORT/callback │
       │                  │    &state=<random>                              │
       │                  │    &code_challenge=<S256(verifier)>             │
       │                  │    &code_challenge_method=S256                  │
       │                  │────────────────────────────────────────────────>│
       │                  │                                                 │
       │                  │              (user authenticates on Platform)   │
       │                  │                                                 │
       │                  │<── redirect to 127.0.0.1:PORT/callback ─────────│
       │                  │      ?code=<auth_code>&state=<state>            │
       │                  │                                                 │
       │  [callback server receives code, validates state]                  │
       │  [returns HTML "Login successful" to browser]                      │
       │                                                                    │
       │  POST /token                                                       │
       │    grant_type=authorization_code                                   │
       │    &client_id=nextflow_cli                                         │
       │    &code=<auth_code>                                               │
       │    &code_verifier=<verifier>                                       │
       │    &redirect_uri=http://127.0.0.1:PORT/callback                    │
       │───────────────────────────────────────────────────────────────────>│
       │<──────────── { access_token, refresh_token } ──────────────────────│
       │                                                                    │
       │  GET /user-info  (Authorization: Bearer <access_token>)            │
       │───────────────────────────────────────────────────────────────────>│
       │<──────────────────── { user info } ────────────────────────────────│
       │                                                                    │
       │  POST /tokens  (generate PAT — same as current flow)               │
       │───────────────────────────────────────────────────────────────────>│
       │<──────────────────── { accessKey: <PAT> } ─────────────────────────│
       │                                                                    │
       │  [save PAT to ~/.nextflow/seqera-auth.config]                      │
```

## Solution

### 1. Register `nextflow_cli` client in Platform

**File:** `platform/tower-config/src/main/resources/application-oauth-client.yml`

```yaml
- client-id: "nextflow_cli"
  client-name: "Nextflow CLI"
  client-type: NATIVE
  client-secret: null
  application-type: "native"
  token-endpoint-auth-method: "none"
  id-token-signed-response-alg: "RS256"
  redirect-uris:
    - "http://127.0.0.1"
  allowed-flows:
    - "authorization_code_pkce"
  allowed-scopes:
    - "openid"
    - "profile"
    - "email"
    - "offline_access"
  third-party: false
```

### 2. New module `lib-platform-oidc` in libseqera

Plain Java production code, Groovy/Spock tests. No external dependencies — uses only JDK classes (`java.net.http.HttpClient`, `com.sun.net.httpserver.HttpServer`, `java.security.*`).

```
lib-platform-oidc/
  build.gradle
  VERSION
  changelog.txt
  src/
    main/java/io/seqera/platform/auth/oidc/
      OidcConfig.java         # authorization_endpoint + token_endpoint
      PkceChallenge.java      # code_verifier + code_challenge
      PkceUtil.java           # PKCE generation helpers
      OidcDiscovery.java      # GET /.well-known/openid-configuration
      OidcCallbackServer.java # Local HTTP server on 127.0.0.1:0
      OidcTokenExchange.java  # POST token endpoint for code exchange
      OidcLoginFlow.java      # Orchestrator
    test/groovy/io/seqera/platform/auth/oidc/
      PkceUtilTest.groovy
      OidcDiscoveryTest.groovy
      OidcCallbackServerTest.groovy
      OidcTokenExchangeTest.groovy
      OidcLoginFlowTest.groovy
```

**`OidcLoginFlow`** — public API:
```java
public class OidcLoginFlow {
    public OidcLoginFlow(String endpoint, String clientId) { ... }
    public String login(Consumer<String> browserLauncher) throws Exception { ... }
}
```

Flow: OIDC discovery → generate PKCE → start callback server (ephemeral port) → invoke `browserLauncher` with authorization URL → wait for callback → exchange code for tokens → return `access_token`.

The `browserLauncher` callback delegates browser-opening to the caller since it's platform-specific.

### 3. Modify `AuthCommandImpl.login()` in Nextflow

**File:** `nextflow/plugins/nf-tower/src/main/io/seqera/tower/plugin/auth/AuthCommandImpl.groovy`

Replace the cloud-vs-enterprise branching with a single OIDC flow:

```groovy
// Was: getCloudEndpointInfo → performAuth0Login / handleEnterpriseAuth
// Now:
performOidcLogin(apiUrl)
```

`performOidcLogin()` creates an `OidcLoginFlow`, gets the OAuth access token, then follows the existing post-auth steps unchanged: `getUserInfo()` → `generatePAT()` → `saveAuthToConfig()` → `config()`.

Remove: `performAuth0Login()`, `requestDeviceAuthorization()`, `pollForDeviceToken()`, `performAuth0Request()`, `handleEnterpriseAuth()`, `promptPAT()`, `getCloudEndpointInfo()`.

### 4. Deprecate Auth0 mappings in PlatformHelper

**File:** `nextflow/modules/nextflow/src/main/groovy/nextflow/platform/PlatformHelper.groovy`

Add `@Deprecated` to `getAuthDomain()` and `getAuthClientId()`.

## Verification

1. `./gradlew :lib-platform-oidc:test` in libseqera
2. `./gradlew :plugins:nf-tower:test` in Nextflow
3. Manual: `./launch.sh auth login -url https://api.cloud.dev-seqera.io` — browser opens, PKCE flow completes, PAT stored in config
4. Backward compat: existing PAT in `seqera-auth.config` continues to work

## Links

- [Platform OIDC provider PR](https://github.com/seqeralabs/platform/pull/10336)
- [OAuth client examples PR](https://github.com/seqeralabs/platform/pull/10473)
