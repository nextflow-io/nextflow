#!/usr/bin/env nextflow

nextflow.enable.types = true

record Vcf {
    id: String
    vcf: Path
}

/*
 * files nested in collections of collections and in tuples
 * should be staged implicitly
 */
process check_nested {
    debug true

    input:
    batches: List<List<Vcf>>
    bams: List<Tuple<String,Path>>
    bam_groups: List<List<Tuple<String,Path>>>
    tuple(id: String, reads: Tuple<Path,Path>)

    script:
    def files = batches.sum([])*.vcf + bams*.getAt(1) + bam_groups.sum([])*.getAt(1) + [reads[0], reads[1]]
    """
    for f in ${files.join(' ')}; do
        [[ "\$f" == "\$(basename \$f)" ]] || { echo "not a staged path: \$f"; exit 1; }
        [[ -e \$f ]] || { echo "file not staged: \$f"; exit 1; }
    done
    echo ${id}: ${files.join(' ')}
    """
}

def make(name: String) -> Path {
    def f = file("${workflow.workDir}/${name}")
    f.parent.mkdirs()
    f.text = name
    return f
}

workflow {
    batches = [
        [record(id: 'a', vcf: make('1/a.vcf')), record(id: 'b', vcf: make('2/b.vcf'))],
        [record(id: 'c', vcf: make('3/c.vcf'))]
    ]
    bams = [tuple('a', make('1/a.bam')), tuple('b', make('2/b.bam'))]
    bam_groups = [[tuple('c', make('3/c.bam'))], [tuple('d', make('4/d.bam'))]]
    reads = tuple(make('1/r_1.fq'), make('1/r_2.fq'))

    check_nested( batches, bams, bam_groups, tuple('s', reads) )
}
