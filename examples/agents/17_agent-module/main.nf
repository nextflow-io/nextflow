/*
 * Agents as modules.
 *
 * The whole agent — its model, instruction, prompt, two skills and its tool — lives
 * in `mods/reporter/`. This script only includes it, under an alias, and calls it.
 *
 * Note there is no `nextflow.enable.types = true` here: the flag is per FILE, and
 * this script declares no typed process of its own. The agent module and the tool
 * module each carry whatever they need.
 */

include { reporter as qc } from './mods/reporter'

workflow {
    qc(channel.of('Assembly for isolate SRR001: N50 = 45 kb, completeness = 96.4%, total length = 5.1 Mb. Should we proceed?'))
        .view { a -> "ANSWER=\n${a}" }
}
