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

package nextflow.lineage.cli

import spock.lang.Narrative
import spock.lang.PendingFeature
import spock.lang.See
import spock.lang.Specification
import spock.lang.Title

/**
 * ADR conformance contract — validator core + difference report shape.
 *
 * Audience: tooling integrators reading or rendering the diff report.
 * Pins D11 (shared LineageValidator core), D13 (difference categories),
 * D15 (one-hop causality enrichment), and D17 informational runtime
 * reporting (the `--strict-runtime` flag itself lives in
 * LineageValidateCliFlagsSpec).
 */
@Title("Lineage Validate — validator core and difference report contract")
@Narrative('''
The shared LineageValidator kernel, the difference-category split, and
the one-hop causality enrichment on diverging outputs. Every method
here is a pin against the corresponding decision in
adr/20260521-lineage-validate.md.
''')
@See([
    "https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md",
    "https://github.com/nextflow-io/nextflow/pull/7167",
    "https://spockframework.org/spock/docs/2.4/all_in_one.html#_pendingfeature"
])
class LineageValidateReportingSpec extends Specification {

    // ───── D11 — Shared LineageValidator core ────────────────────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d11--cli-vs-lineagesnapshotter-extract-shared-lineagevalidator")
    def 'CLI and LineageSnapshotter delegate to one LineageValidator core (no drift)'() {
        expect: false
    }

    // ───── D13 — Difference categories ───────────────────────────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d13--difference-categories-outputs--params--workflow-identity--resources")
    def 'differences are bucketed into outputs / params / workflow-identity / resources / runtime / status'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d13--difference-categories-outputs--params--workflow-identity--resources")
    def 'outputs, params, workflow-identity, status fail validation; resources, runtime are informational'() {
        expect: false
    }

    // ───── D15 — One-hop causality enrichment ────────────────────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d15--causality-in-failure-reports-one-hop-upstream")
    def 'diverging outputs carry a one-hop cause block: producing TaskRun + differing inputs'() {
        expect: false
    }

    // ───── D17 — Runtime identity reporting (informational tier) ─────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d17--runtime-identity-informational-opt-in-strictness")
    def 'container / conda / executor drift is reported but informational by default'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d17--runtime-identity-informational-opt-in-strictness")
    def 'runtime fields not yet recorded by the model report as "unknown" rather than equivalent'() {
        expect: false
    }
}
