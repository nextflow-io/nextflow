// One process per hash key, so a cache miss localises immediately.
// See ../runbook.md for how this is run.
//
// NOTE on P_BAG_INPUT: the plan called for a second H1-vs-H2 discriminator
// exercising a multi-file `path` input (delivered as an ArrayBag<FileHolder>)
// to cover `cacheFunnelFirst`, on the premise that ArrayBag implements both
// `Bag` and `CacheFunnel`. Verified against
// modules/nextflow/src/main/groovy/nextflow/util/ArrayBag.groovy: it does not
// implement CacheFunnel (only `Bag<E>, List<E>, KryoSerializable`), and no
// type in the codebase implements both a Map/Bag/Set and CacheFunnel at once
// (checked FileHolder, GroupKey, SecretImpl, CmdLineOptionMap, PluginRef,
// AzPoolOpts/AzFileShareOpts). `HashBuilder.with()` only lets `cacheFunnelFirst`
// change behaviour for such a dual object; against a Bag of plain FileHolders
// the funnel-vs-collection branch order is unobservable either way. Separately,
// `StdSpecs` never varies `cacheFunnelFirst` independently of
// `orderIndependentMaps` (EncodingRules.LEGACY sets both false,
// RECORD_TYPES sets both true), so no pair of shipped specs could isolate it
// even with a dual object. A P_BAG_INPUT process was therefore not added: it
// would not discriminate anything. See runbook.md section 3 for the
// consequence for std/v1-vs-std/v2 coverage.

include { P_MODULE_BUNDLE } from './modules/bundled/main.nf'

params.greeting = 'hello'

process P_BASIC {
    ext flavour: 'vanilla'

    input:
    val word
    path sample

    output:
    stdout

    script:
    """
    echo "${word} ${params.greeting} ${task.ext.flavour} \$(cat ${sample})"
    """
}

// H1 vs H2 discriminator: a genuine Map reaches the hasher, so the encoding
// rules (values() vs entrySet()) actually bite.
process P_MAP_INPUT {
    input:
    val settings

    output:
    stdout

    script:
    """
    echo "${settings.alpha} ${settings.beta}"
    """
}

// H3 vs H4 discriminator: without an eval output the two specs are identical.
process P_EVAL {
    output:
    stdout
    eval 'echo evaluated', emit: probe

    script:
    """
    echo eval-process
    """
}

process P_BIN {
    output:
    stdout

    script:
    """
    helper.sh one
    """
}

process P_CONTAINER {
    container 'quay.io/nextflow/bash@sha256:bea0e244b7c5367b2b0de687e7d28f692013aa18970941c7dd184450125163ac'

    output:
    stdout

    script:
    """
    echo containerised
    """
}

process P_CONDA {
    conda 'bioconda::fastqc=0.12.1'

    output:
    stdout

    script:
    """
    echo conda-backed
    """
}

process P_STUB {
    output:
    stdout

    script:
    """
    echo real
    """

    stub:
    """
    echo stubbed
    """
}

workflow {
    P_BASIC(Channel.value('word'), file("${projectDir}/data/input.txt"))
    P_MAP_INPUT(Channel.value([alpha: 'one', beta: 'two']))
    P_EVAL()
    P_BIN()
    P_MODULE_BUNDLE()
    P_CONTAINER()
    P_CONDA()
    P_STUB()
}
