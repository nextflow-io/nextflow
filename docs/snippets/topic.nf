nextflow.preview.topic = true

Channel.of(1,2,3) | topic('foo')
Channel.topic('foo').view()