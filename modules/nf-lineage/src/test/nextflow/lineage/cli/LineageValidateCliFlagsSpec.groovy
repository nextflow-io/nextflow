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
import spock.lang.Unroll

/**
 * ADR conformance contract — CLI flags, env vars, exit codes.
 *
 * Audience: pipeline authors wiring `nextflow lineage validate` into CI.
 * Pins the user-facing surface for D1 (CI-first defaults), D6 (output
 * formats + exit codes), D7 (CI/agent env auto-detect), D8 (ignore
 * mechanism), and the `--strict-runtime` flag from D17.
 *
 * See LineageValidateEquivalenceSpec, LineageValidateBaselineSpec, and
 * LineageValidateReportingSpec for the other ADR decision groups.
 */
@Title("Lineage Validate — CLI flags, env vars, exit codes contract")
@Narrative('''
The flags, environment variables, and exit codes a pipeline author sees.
Every method here is a @PendingFeature pin against the corresponding
decision in adr/20260521-lineage-validate.md; replacing a stub with a
real spec is how an ADR decision graduates from "documented" to "shipped".
''')
@See([
    "https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md",
    "https://github.com/nextflow-io/nextflow/pull/7167",
    "https://spockframework.org/spock/docs/2.4/all_in_one.html#_pendingfeature"
])
class LineageValidateCliFlagsSpec extends Specification {

    // ───── D1 — CI-first defaults ────────────────────────────────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d1--primary-use-case-ci-led-three-use-case-design")
    def 'defaults optimise for CI: fast pass/fail, scriptable exit codes, structured output without flags'() {
        expect: false
    }

    // ───── D6 — Output formats: human / json / summary ───────────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d6--failure-output-human-diff-default---json-for-ci")
    def 'default failure output is a JGit-style unified diff grouped by category'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d6--failure-output-human-diff-default---json-for-ci")
    def '--json emits a stable structured difference schema'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d6--failure-output-human-diff-default---json-for-ci")
    def '--summary emits one-line counts (e.g. "5 outputs differ, 1 param differs")'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d7--auto-detect-ci-environments")
    def '--format human|json|summary is the explicit override for the env-driven default'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d6--failure-output-human-diff-default---json-for-ci")
    def '--trace switches to debug mode with deeper upstream walk on diverging outputs'() {
        expect: false
    }

    // ───── D6 — Exit codes ───────────────────────────────────────────────

    @Unroll
    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d6--failure-output-human-diff-default---json-for-ci")
    def 'exit code is #code when #scenario'() {
        expect: false
        where:
        code | scenario
        0    | 'runs are semantically equivalent'
        1    | 'runs differ in any failing category'
        2    | 'load failure, schema mismatch, or resolver error'
    }

    // ───── D7 — Auto-detect CI / agent environments ──────────────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d7--auto-detect-ci-environments")
    def 'CI=true switches default output mode to --json'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d7--auto-detect-ci-environments")
    def 'GITHUB_ACTIONS=true emits ::error:: annotations and writes to GITHUB_STEP_SUMMARY'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d7--auto-detect-ci-environments")
    def 'AGENT=1 or CLAUDECODE=1 switches default to one-line summary plus JSON tail'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d7--auto-detect-ci-environments")
    def 'precedence is explicit flag > environment variable > nextflow.config > built-in default'() {
        expect: false
    }

    // ───── D8 — Ignore mechanism: flag + config, dotted paths ────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d8--ignore-mechanism-flag--config-jsonpath-style")
    def '--ignore-fields accepts dotted paths (params.outdir, outputs.*.modifiedAt)'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d8--ignore-mechanism-flag--config-jsonpath-style")
    def 'wildcards * (single segment) and ** (recursive) are supported in ignore paths'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d8--ignore-mechanism-flag--config-jsonpath-style")
    def 'lineage.validate.ignore list in nextflow.config is honoured'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d8--ignore-mechanism-flag--config-jsonpath-style")
    def 'built-in EPHEMERAL_FIELDS always apply, unaffected by user config'() {
        expect: false
    }

    // ───── D17 — Strict runtime flag ─────────────────────────────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d17--runtime-identity-informational-opt-in-strictness")
    def '--strict-runtime promotes the runtime category to fail-if-different'() {
        expect: false
    }
}
