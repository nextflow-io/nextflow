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
 * ADR conformance contract — baseline resolution and schema.
 *
 * Audience: plugin authors implementing a LineageResolver SPI (e.g.
 * `nf-tower`'s `tower://` resolver). Pins D4 (peer LID + snapshot file +
 * --save-snapshot), D10 (schema evolution), D12 (LineageResolver SPI
 * dispatch).
 */
@Title("Lineage Validate — baseline resolution and schema contract")
@Narrative('''
How `--against` references are resolved into trees, how schema versions
are reconciled, and which resolvers the SPI must support. Every method
here is a pin against the corresponding decision in
adr/20260521-lineage-validate.md.
''')
@See([
    "https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md",
    "https://github.com/nextflow-io/nextflow/pull/7167",
    "https://spockframework.org/spock/docs/2.4/all_in_one.html#_pendingfeature"
])
class LineageValidateBaselineSpec extends Specification {

    // ───── D4 — Baseline model: peer LID + snapshot file ─────────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d4--baseline-model-peer-lid-and-snapshot-file-from-day-one")
    def '--against accepts a peer LID baseline (lid://...)'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d4--baseline-model-peer-lid-and-snapshot-file-from-day-one")
    def '--against accepts a local snapshot file baseline'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d4--baseline-model-peer-lid-and-snapshot-file-from-day-one")
    def '--save-snapshot writes the normalized tree of the new run to disk'() {
        expect: false
    }

    // ───── D10 — Schema evolution ────────────────────────────────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d10--schema-evolution-refuse-on-major-mismatch-normalize-on-minor")
    def 'identical schemaVersion: direct comparison'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d10--schema-evolution-refuse-on-major-mismatch-normalize-on-minor")
    def 'minor schemaVersion mismatch: normalise both sides, drop newer-only fields'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d10--schema-evolution-refuse-on-major-mismatch-normalize-on-minor")
    def 'major schemaVersion mismatch: hard error with exit code 2 and migration message'() {
        expect: false
    }

    // ───── D12 — Pluggable LineageResolver SPI ───────────────────────────

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d12--platform-integration-pluggable-lineageresolver-spi")
    def 'lid:// resolves via the built-in LinStore resolver'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d12--platform-integration-pluggable-lineageresolver-spi")
    def 'file paths resolve via the built-in SnapshotFileResolver'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d12--platform-integration-pluggable-lineageresolver-spi")
    def 'tower:// fails with an install hint when nf-tower plugin is absent'() {
        expect: false
    }

    @PendingFeature
    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d12--platform-integration-pluggable-lineageresolver-spi")
    def 'nf-tower contributes a tower:// resolver that joins Platform API on task hash'() {
        expect: false
    }
}
