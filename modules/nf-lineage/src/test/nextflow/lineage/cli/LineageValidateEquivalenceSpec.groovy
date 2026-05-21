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
 * ADR conformance contract — semantic equivalence rules.
 *
 * Audience: reviewers debating "what does equivalent mean?". Pins the
 * equivalence model: D2 (outputs + key inputs), D3 (checksum-only file
 * equivalence), D5 (subworkflow flattening), D9 (workflow identity),
 * D14 (output join key), D16 (failed runs), D18 (non-file outputs gap).
 */
@Title("Lineage Validate — semantic equivalence rules contract")
@Narrative('''
The equivalence model: which fields participate, which are ignored, and
what makes two recorded runs "the same". Every method here is a pin
against the corresponding decision in adr/20260521-lineage-validate.md;
changing equivalence semantics requires updating both the ADR and the
spec replacing the stub here.
''')
@See([
    "https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md",
    "https://github.com/nextflow-io/nextflow/pull/7167",
    "https://spockframework.org/spock/docs/2.4/all_in_one.html#_pendingfeature"
])
class LineageValidateEquivalenceSpec extends Specification {

    // ───── D2 — Equivalence unit: outputs + key inputs ───────────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d2--equivalence-unit-outputs--key-inputs")
    def 'two runs with identical key inputs and published outputs validate as equivalent'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d2--equivalence-unit-outputs--key-inputs")
    def 'refactors that change internal task graph but preserve outputs do not fail validation'() {
        expect: false
    }

    // ───── D3 — File equivalence: checksum-only ──────────────────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d3--file-equivalence-checksum-only")
    def 'FileOutput comparison uses recorded checksum only — no rehashing, no type-aware diff'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d3--file-equivalence-checksum-only")
    def 'volatile checksums can be excluded via --ignore-fields outputs.<name>.checksum'() {
        expect: false
    }

    // ───── D5 — Subworkflow handling: flatten to terminal outputs ────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d5--subworkflow-handling-flatten-to-terminal-outputs")
    def 'sub-workflow restructuring with identical terminal outputs validates as equivalent'() {
        expect: false
    }

    // ───── D9 — Workflow identity: commitId + repository ─────────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d9--workflow-identity-commitid--repository")
    def 'workflow identity is repository plus commitId; scriptFiles are skipped when commitId is present'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d9--workflow-identity-commitid--repository")
    def 'dirty-tree fallback compares scriptFiles checksums when commitId is absent'() {
        expect: false
    }

    // ───── D14 — Output join key ─────────────────────────────────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d14--output-join-key-relative-path-under--outputdir")
    def 'two FileOutput records join on relative path under -outputDir'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d14--output-join-key-relative-path-under--outputdir")
    def 'longest-common-suffix fallback (with warning) applies when -outputDir is not recorded'() {
        expect: false
    }

    // ───── D16 — Failed runs: compare what exists ────────────────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d16--failed-runs-compare-what-exists-surface-status-divergence")
    def 'a FAILED run can be compared against a SUCCESS run; status divergence is reported'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d16--failed-runs-compare-what-exists-surface-status-divergence")
    def 'missing task runs or outputs on one side are reported in their categories'() {
        expect: false
    }

    // ───── D18 — Non-file outputs (model gap) ────────────────────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d18--non-file-outputs-model-gap-scope-honestly")
    def 'value emits are out of scope until the lineage model adds ValueOutput'() {
        expect: false
    }
}
