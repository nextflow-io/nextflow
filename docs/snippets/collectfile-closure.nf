Channel.of('Hola', 'Ciao', 'Hello', 'Bonjour', 'Halo')
    .collectFile { item ->
        [ "${item[0]}.txt", item + '\n' ]
    }
    .subscribe {
        println "File '${it.name}' contains:"
        println it.text
    }