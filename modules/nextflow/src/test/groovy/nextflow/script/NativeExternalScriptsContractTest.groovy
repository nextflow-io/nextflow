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

package nextflow.script

import spock.lang.Narrative
import spock.lang.PendingFeature
import spock.lang.See
import spock.lang.Specification
import spock.lang.Title
import spock.lang.Unroll

/**
 * Executable contract for the native external scripts ADR.
 *
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
@Title('Native external scripts run unchanged and receive typed Nextflow task context')
@Narrative('''
The native external scripts ADR proposes a replacement for template-time source
rewriting: external Bash, Python, R, Julia, and related scripts should be staged
as real native-language files, while Nextflow injects task context through
sidecars and language-specific runtime namespaces.

The contract is agent-facing as well as human-facing: module-local scripts must
be easy to inspect, edit, lint, and execute through a fast Nextflow feedback
loop instead of forcing agents back into heredocs or preprocessed templates.

This contract is intentionally pending. Each feature method names one design
promise from the ADR so implementation work can replace the stub with a real
assertion and remove @PendingFeature when the behavior lands.
''')
@See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
class NativeExternalScriptsContractTest extends Specification {

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'external script source is staged into the task work directory unchanged'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'native script execution does not perform Groovy or GString template expansion on the source file'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'task context is serialized to a canonical sidecar before the native script runs'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'Bash scripts receive process context through associative arrays instead of rewritten variables'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'Python scripts can import a nextflow namespace object backed by the task sidecar'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'generated Python type stubs specialize process inputs outputs params and resources'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'R scripts receive an idiomatic nextflow context object backed by the task sidecar'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'native external scripts are staged through module and workflow resource bundles'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'native external scripts are staged as code assets not declared input files'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'task inputFiles remains limited to declared process data inputs'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'module bin helpers and native scripts are staged as one coherent runtime bundle'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'nextflow module run provides a fast feedback loop for module-local script changes'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'script failures report the original external script path and line number'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'script-relative assets remain resolvable from the external script body'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'legacy template expansion remains available for existing template users'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'Rust scripts can be supported later through rust-script and generated typed context structs'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'R Markdown documents can be supported later as report-oriented R script assets'() {
        expect:
        false
    }

    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def 'Jupyter notebooks can be supported later through the same sidecar context contract'() {
        expect:
        false
    }

    @Unroll
    @PendingFeature
    @See('https://github.com/nextflow-io/nextflow/blob/master/adr/20260522-native-external-scripts.md')
    def '#language native scripts expose the common process context fields'() {
        expect:
        false

        where:
        language | fields
        'Bash'   | ['input', 'output', 'params', 'resources', 'cpus', 'attempt']
        'Python' | ['input', 'output', 'params', 'resources', 'cpus', 'attempt']
        'R'      | ['input', 'output', 'params', 'resources', 'cpus', 'attempt']
        'Julia'  | ['input', 'output', 'params', 'resources', 'cpus', 'attempt']
        'Rust'   | ['input', 'output', 'params', 'resources', 'cpus', 'attempt']
    }
}
