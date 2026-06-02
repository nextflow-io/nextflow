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

package io.seqera.executor

import nextflow.util.MemoryUnit
import spock.lang.Specification

/**
 * Validates capture of resource-derived command values from the live command GString
 * plus the raw script source, with no re-execution of the user script.
 *
 * <p>Each test builds a real {@code command} GString (what {@code TaskRun.scriptGString}
 * would hold, with the requested values already interpolated) and the matching raw
 * {@code source} (a non-interpolating triple-single-quoted string, so {@code ${task.x}}
 * stays verbatim) — exactly the two artifacts Nextflow presents at materialization.
 *
 * @author Paolo Di Tommaso
 */
class ResourceInterpolatorTest extends Specification {

    private static final Set<String> TASK_REF = ['task'] as Set
    private static final ResourceInterpolator.Resource MEM = ResourceInterpolator.Resource.MEMORY
    private static final ResourceInterpolator.Resource CPU = ResourceInterpolator.Resource.CPUS

    /** Substitute each binding's original value back into the placeholdered script. */
    private static String restore(ResourceInterpolator.Result res) {
        String s = res.script
        for( ResourceInterpolator.Binding b : res.bindings )
            s = s.replace('${' + b.name + '}', b.value.toPlainString())
        return s
    }

    // ---- direct self-sizing values -----------------------------------------

    def 'placeholders a direct self-sizing memory expression (the sched#492 bcftools case)'() {
        given:
        def reqMem = MemoryUnit.of('36864 MB')
        GString command = "bcftools sort --max-mem ${(reqMem.toMega() * 0.9) as int}M out.vcf"
        def source = '''"""bcftools sort --max-mem ${(task.memory.toMega()*0.9) as int}M out.vcf"""'''

        expect: 'the request is baked into the rendered command'
        command.toString() == 'bcftools sort --max-mem 33177M out.vcf'

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, 6, reqMem)

        then:
        res.changed()
        res.script == 'bcftools sort --max-mem ${SEQERA_TASK_MEM_0}M out.vcf'
        res.bindings.size() == 1
        res.bindings[0].resource == MEM
        res.bindings[0].value == 33177.0
        restore(res) == command.toString()
    }

    def 'placeholders a realistic multi-line command (continuations + surrounding interpolations)'() {
        given:
        def reqMem = MemoryUnit.of('36864 MB')
        def prefix = 'SAMPLE_21'
        def vcf = 'SAMPLE_21.vcf'
        GString command = """\
            bcftools sort \\
                --output ${prefix}.vcf.gz \\
                --max-mem ${(reqMem.toMega() * 0.9) as int}M \\
                --temp-dir . \\
                ${vcf}
            """
        def source = '''\
            """
            bcftools sort \\
                --output ${prefix}.vcf.gz \\
                --max-mem ${(task.memory.toMega() * 0.9) as int}M \\
                --temp-dir . \\
                ${vcf}
            """
            '''

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, 6, reqMem)

        then: 'only the memory value is placeholdered; prefix/vcf are left untouched'
        res.changed()
        res.bindings.size() == 1
        res.bindings[0].resource == MEM
        res.bindings[0].value == 33177.0
        res.script.contains('--max-mem ${SEQERA_TASK_MEM_0}M')
        res.script.contains('--output SAMPLE_21.vcf.gz')
        !res.script.contains('33177')
        restore(res) == command.toString()
    }

    def 'placeholders a JVM -Xmx flag sized from giga (unit-agnostic)'() {
        given:
        def reqMem = MemoryUnit.of('32 GB')
        GString command = "java -Xmx${reqMem.toGiga()}g -jar picard.jar"
        def source = '''"""java -Xmx${task.memory.toGiga()}g -jar picard.jar"""'''

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, 4, reqMem)

        then:
        res.changed()
        res.script == 'java -Xmx${SEQERA_TASK_MEM_0}g -jar picard.jar'
        res.bindings[0].resource == MEM
        res.bindings[0].value == 32.0
        restore(res) == command.toString()
    }

    def 'placeholders a memory value rendered in bytes (large number)'() {
        given:
        def reqMem = MemoryUnit.of('32 GB')
        GString command = "STAR --limitBAMsortRAM ${reqMem.toBytes()}"
        def source = '''"""STAR --limitBAMsortRAM ${task.memory.toBytes()}"""'''

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, 4, reqMem)

        then:
        res.changed()
        res.script == 'STAR --limitBAMsortRAM ${SEQERA_TASK_MEM_0}'
        res.bindings[0].value == new BigDecimal(reqMem.toBytes())
        restore(res) == command.toString()
    }

    def 'placeholders a cpus value with arithmetic'() {
        given:
        def reqCpus = 4
        GString command = "bwa mem -t ${reqCpus * 2} ref.fa reads.fq"
        def source = '''"""bwa mem -t ${task.cpus * 2} ref.fa reads.fq"""'''

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, reqCpus, MemoryUnit.of('8 GB'))

        then: 'the rendered value (8) is captured; the launch-time ratio preserves the x2'
        res.changed()
        res.script == 'bwa mem -t ${SEQERA_TASK_CPUS_0} ref.fa reads.fq'
        res.bindings[0].resource == CPU
        res.bindings[0].value == 8.0
        restore(res) == command.toString()
    }

    def 'placeholders both cpus and memory in a multi-line command'() {
        given:
        def reqMem = MemoryUnit.of('16 GB')
        def reqCpus = 8
        GString command = """\
            mytool \\
                --threads ${reqCpus} \\
                --memory ${reqMem.toMega()}
            """
        def source = '''\
            """
            mytool \\
                --threads ${task.cpus} \\
                --memory ${task.memory.toMega()}
            """
            '''

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, reqCpus, reqMem)

        then:
        res.changed()
        res.bindings.size() == 2
        res.bindings.find { it.resource == CPU }.value == 8.0
        res.bindings.find { it.resource == MEM }.value == 16384.0
        res.script.contains('--threads ${SEQERA_TASK_CPUS_0}')
        res.script.contains('--memory ${SEQERA_TASK_MEM_0}')
        restore(res) == command.toString()
    }

    def 'placeholders two distinct memory references in one command'() {
        given:
        def reqMem = MemoryUnit.of('8 GB')
        GString command = "tool --heap ${reqMem.toMega()} --buffer ${reqMem.toMega()}"
        def source = '''"""tool --heap ${task.memory.toMega()} --buffer ${task.memory.toMega()}"""'''

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, 4, reqMem)

        then: 'each reference gets its own placeholder'
        res.changed()
        res.script == 'tool --heap ${SEQERA_TASK_MEM_0} --buffer ${SEQERA_TASK_MEM_1}'
        res.bindings.size() == 2
        res.bindings.every { it.resource == MEM && it.value == 8192.0 }
        restore(res) == command.toString()
    }

    def 'placeholders a value inside double-quotes (which still expands in bash)'() {
        given:
        def reqMem = MemoryUnit.of('4 GB')
        GString command = "sh -c \"java -Xmx${reqMem.toGiga()}g\" && echo done"
        def source = '''"""sh -c "java -Xmx${task.memory.toGiga()}g" && echo done"""'''

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, 2, reqMem)

        then: 'inside double-quotes ${VAR} still expands, so the rewrite is applied'
        res.changed()
        res.script == 'sh -c "java -Xmx${SEQERA_TASK_MEM_0}g" && echo done'
        restore(res) == command.toString()
    }

    // ---- indirection through local variables --------------------------------

    def 'follows one level of local-variable indirection (skesa: ternary + nested GString)'() {
        given:
        def reqMem = MemoryUnit.of('200 MB')
        def reqCpus = 8
        def memory = reqMem ? "--memory ${reqMem.toMega()} MB" : ''
        GString command = "skesa --cores ${reqCpus} ${memory}"
        def source = '''\
            def memory = task.memory ? "--memory ${task.memory.toMega()} MB" : ''
            """skesa --cores ${task.cpus} ${memory}"""
            '''

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, reqCpus, reqMem)

        then: 'both the direct cpus value and the indirected memory value are placeholdered'
        res.changed()
        res.script == 'skesa --cores ${SEQERA_TASK_CPUS_0} --memory ${SEQERA_TASK_MEM_0} MB'
        res.bindings.find { it.resource == CPU }.value == 8.0
        res.bindings.find { it.resource == MEM }.value == 200.0
        restore(res) == command.toString()
    }

    def 'limitation: a scalar (non-GString) local derived from a resource is left unchanged (safe)'() {
        given: 'def avail = task.memory.toMega() is a number, not a nested GString'
        def reqMem = MemoryUnit.of('36864 MB')
        def avail = (reqMem.toMega() * 0.9) as int
        GString command = "bcftools sort --max-mem ${avail}M out.vcf"
        def source = '''\
            def avail = (task.memory.toMega() * 0.9) as int
            """bcftools sort --max-mem ${avail}M out.vcf"""
            '''

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, 6, reqMem)

        then: 'not optimized, but never corrupted — identical to today'
        !res.changed()
        res.script == command.toString()
    }

    // ---- joint expressions & no-ops -----------------------------------------

    def 'leaves a joint memory+cpus expression untouched, placeholders the pure cpus value'() {
        given:
        def reqMem = MemoryUnit.of('16 GB')
        def reqCpus = 4
        GString command = """\
            samtools sort \\
                -@ ${reqCpus} \\
                -m ${(reqMem.toMega() / reqCpus) as int}M \\
                -o out.bam in.bam
            """
        def source = '''\
            """
            samtools sort \\
                -@ ${task.cpus} \\
                -m ${(task.memory.toMega() / task.cpus) as int}M \\
                -o out.bam in.bam
            """
            '''

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, reqCpus, reqMem)

        then: 'per-thread memory depends on both resources -> left literal; only -@ is placeholdered'
        res.changed()
        res.bindings.size() == 1
        res.bindings[0].resource == CPU
        res.script.contains('-@ ${SEQERA_TASK_CPUS_0}')
        res.script.contains('-m 4096M')
        restore(res) == command.toString()
    }

    def 'no-op when no command value derives from resources'() {
        given:
        GString command = "echo ${'hello'} attempt=${1}"
        def source = '''"""echo ${greeting} attempt=${task.attempt}"""'''

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, 4, MemoryUnit.of('4 GB'))

        then:
        !res.changed()
        res.script == 'echo hello attempt=1'
    }

    def 'gate: skips when the script does not reference task'() {
        given:
        def reqMem = MemoryUnit.of('1 GB')
        GString command = "do --mem ${reqMem.toMega()}"
        def source = '''"""do --mem ${someVar}"""'''

        when:
        def res = ResourceInterpolator.interpolate(command, source, ['someVar'] as Set, 1, reqMem)

        then:
        !res.changed()
        res.script == command.toString()
    }

    def 'does not match longer property names like task.memoryUsage'() {
        given:
        GString command = "echo ${999}"
        def source = '''"""echo ${task.memoryUsage}"""'''

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, 2, MemoryUnit.of('1 GB'))

        then:
        !res.changed()
        res.script == 'echo 999'
    }

    // ---- safety fallbacks ---------------------------------------------------

    def 'falls back when the placeholder would land inside shell single-quotes (would not expand)'() {
        given:
        def reqMem = MemoryUnit.of('512 MB')
        GString command = "tool --opt 'mem=${reqMem.toMega()}'"
        def source = '''"""tool --opt 'mem=${task.memory.toMega()}'"""'''

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, 2, reqMem)

        then: 'a ${VAR} inside single-quotes is inert in bash -> leave unchanged'
        !res.changed()
        res.script == "tool --opt 'mem=512'"
    }

    def 'consistency guard: falls back when the source template does not match the command'() {
        given: 'the source template parses to fewer interpolations than the command has'
        def reqMem = MemoryUnit.of('1 GB')
        GString command = "x ${reqMem.toMega()} y ${reqMem.toMega()}"
        def source = '''"""x ${task.memory.toMega()} y"""'''

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, 1, reqMem)

        then:
        !res.changed()
        res.script == command.toString()
    }

    def 'falls back for unrecognized (single double-quoted) command forms'() {
        given: 'lastStringLiteral only handles """ and $/.../$'
        def reqMem = MemoryUnit.of('1 GB')
        GString command = "do --mem ${reqMem.toMega()}"
        def source = '''"do --mem ${task.memory.toMega()}"'''

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, 1, reqMem)

        then:
        !res.changed()
    }

    def 'no-op when there is no command GString'() {
        expect:
        !ResourceInterpolator.interpolate(null, 'whatever', TASK_REF, 4, MemoryUnit.of('4 GB')).changed()
    }

    // ---- alternate command forms --------------------------------------------

    def 'supports dollar-slashy command templates'() {
        given:
        def reqMem = MemoryUnit.of('512 MB')
        GString command = "tool --mem ${reqMem.toMega()}"
        def source = '''\
            def x = 1
            $/tool --mem ${task.memory.toMega()}/$
            '''

        when:
        def res = ResourceInterpolator.interpolate(command, source, TASK_REF, 2, reqMem)

        then:
        res.changed()
        res.script == 'tool --mem ${SEQERA_TASK_MEM_0}'
        restore(res) == command.toString()
    }

    // ---- low-level helpers --------------------------------------------------

    def 'parses interpolation expressions in order'() {
        expect:
        ResourceInterpolator.parseInterpolations('a ${task.cpus} b $task.memory c') == ['task.cpus', 'task.memory']
        ResourceInterpolator.parseInterpolations('--max-mem ${(task.memory.toMega()*0.9) as int}M') == ['(task.memory.toMega()*0.9) as int']
        ResourceInterpolator.parseInterpolations('no interpolation here') == []
    }

    def 'placeholderInSingleQuotes detects only single-quoted placeholders'() {
        expect:
        // single-quoted Groovy literals keep ${...} verbatim (no GString interpolation)
        ResourceInterpolator.placeholderInSingleQuotes('a \'${SEQERA_TASK_MEM_0}\' b')
        !ResourceInterpolator.placeholderInSingleQuotes('a "X${SEQERA_TASK_MEM_0}" b')
        !ResourceInterpolator.placeholderInSingleQuotes('a ${SEQERA_TASK_MEM_0} b')
    }
}
