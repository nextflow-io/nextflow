Channel.of('alpha', 'beta', 'gamma')
    .collectFile(name: 'sample.txt', newLine: true)
    .subscribe { file ->
        println "Entries are saved to file: $file"
        println "File content is: ${file.text}"
    }