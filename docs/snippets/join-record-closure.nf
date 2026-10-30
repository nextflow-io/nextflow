samples = channel.of( record(sample: 'X', genome: 'g1'), record(sample: 'Y', genome: 'g2'), record(sample: 'Z', genome: 'g3') )
genomes = channel.of( record(id: 'g3', ref: 6), record(id: 'g2', ref: 5), record(id: 'g1', ref: 4) )

samples.join(genomes, by: { sample, genome -> sample.genome == genome.id }).view()
