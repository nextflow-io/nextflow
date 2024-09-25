Channel.of('alpha', 'beta', 'gamma')
    .collectFile(name: 'sample.txt', newLine: true)
    .subscribe { txt ->
        println "Entries are saved to file: $txt"
        println "File content is: ${txt.text}"
    }