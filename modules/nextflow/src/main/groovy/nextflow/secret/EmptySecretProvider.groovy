/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.secret

import groovy.transform.CompileStatic

/**
 * A mock secrets provider that returns empty strings for all secret requests.
 *
 * <p>This provider is a critical component of Nextflow's 2-phase configuration loading
 * strategy, which solves the chicken-and-egg problem between configuration parsing
 * and plugin loading.
 *
 * <h2>The Problem</h2>
 * Configuration files may reference secrets provided by plugins (e.g., AWS secrets),
 * but plugins are loaded AFTER configuration parsing. This creates a dependency cycle:
 * <ul>
 *   <li>Config parsing needs secret values to complete</li>
 *   <li>Plugin loading needs config to determine which plugins to load</li>
 *   <li>Secret providers are registered by plugins</li>
 * </ul>
 *
 * <h2>The Solution: 2-Phase Configuration Loading</h2>
 * <pre>
 * ┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
 * │   PHASE 1       │    │   BETWEEN       │    │   PHASE 2       │
 * │   (Mock Mode)   │    │   PHASES        │    │   (Real Mode)   │
 * ├─────────────────┤    ├─────────────────┤    ├─────────────────┤
 * │ EmptySecretProv │    │ Load Plugins    │    │ Real Secrets    │
 * │ secrets → ""    │    │ Load Secret     │    │ secrets → value │
 * │ Config: fallback│ -> │ Providers       │ -> │ Config: actual  │
 * │ accessed=true   │    │                 │    │ values          │
 * └─────────────────┘    └─────────────────┘    └─────────────────┘
 * </pre>
 *
 * <h3>Phase 1: Mock Configuration Loading</h3>
 * <ol>
 *   <li>EmptySecretProvider is used as the secrets provider</li>
 *   <li>Configuration is parsed normally</li>
 *   <li>When {@code secrets.SOME_NAME} is referenced, this provider returns empty string</li>
 *   <li>The {@code accessed} flag is set to {@code true}</li>
 *   <li>Config must use defensive patterns: {@code secrets.FOO ? "value-${secrets.FOO}" : "fallback"}</li>
 *   <li>Result: Configuration parses successfully with fallback values</li>
 * </ol>
 *
 * <h3>Between Phases: Plugin and Secrets Loading</h3>
 * <ol>
 *   <li>Plugins are loaded using the successfully parsed configuration</li>
 *   <li>Plugin-provided secrets providers are registered</li>
 *   <li>The real secrets loading system is initialized</li>
 * </ol>
 *
 * <h3>Phase 2: Real Configuration Loading (Conditional)</h3>
 * <ol>
 *   <li>Check if {@code usedSecrets()} returns {@code true}</li>
 *   <li>If secrets were accessed in Phase 1, reload the entire configuration</li>
 *   <li>This time, real secret providers are available</li>
 *   <li>Same config expressions now resolve with actual secret values</li>
 *   <li>Result: Configuration contains real secret values instead of fallbacks</li>
 * </ol>
 *
 * <h2>Example Usage</h2>
 * Given a config file:
 * <pre>
 * outputDir = secrets.MY_SECRET ? "results-${secrets.MY_SECRET}" : "results"
 * </pre>
 *
 * <b>Phase 1:</b> {@code secrets.MY_SECRET} returns {@code ""} → {@code outputDir = "results"}
 * <b>Phase 2:</b> {@code secrets.MY_SECRET} returns {@code "hello-world"} → {@code outputDir = "results-hello-world"}
 *
 * <h2>Key Benefits</h2>
 * <ul>
 *   <li><b>No parsing failures:</b> Configuration always parses successfully</li>
 *   <li><b>Plugin compatibility:</b> Supports plugin-provided secret providers</li>
 *   <li><b>Performance:</b> Only reloads config when secrets are actually used</li>
 *   <li><b>Defensive configs:</b> Forces robust configuration patterns</li>
 * </ul>
 *
 * <h2>Important Notes</h2>
 * <ul>
 *   <li>This provider does NOT remember which specific secrets were accessed</li>
 *   <li>It only tracks WHETHER any secrets were accessed at all</li>
 *   <li>Configuration files MUST handle empty secret values gracefully</li>
 *   <li>The 2-phase loading is transparent to the user</li>
 * </ul>
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 * @see nextflow.cli.CmdRun#run()
 * @see nextflow.secret.SecretsLoader
 */
@CompileStatic
class EmptySecretProvider extends NullProvider {

    private boolean accessed

    @Override
    Secret getSecret(String name) {
        accessed = true
        return new SecretImpl(name, '')
    }

    boolean usedSecrets() {
        return accessed
    }
}
