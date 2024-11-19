Channel.of('Hola', 'Ciao', 'Hello', 'Bonjour', 'Halo')
    .collectFile { item ->
        [ "${item[0]}.txt", item + '\n' ]
    }
    .subscribe { file ->
        println "File '${file.name}' contains:"
        println file.text
    }