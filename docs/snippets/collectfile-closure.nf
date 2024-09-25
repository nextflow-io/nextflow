Channel.of('Hola', 'Ciao', 'Hello', 'Bonjour', 'Halo')
    .collectFile { item ->
        [ "${item[0]}.txt", item + '\n' ]
    }
    .subscribe { txt ->
        println "File '${txt.name}' contains:"
        println txt.text
    }