   source = Channel.from 'Hello world', 'Hola', 'Hello John'
    queue1 = Channel.create()
    queue2 = Channel.create()
    
    source.choice(queue1, queue2) { a -> a =~ /^Hello.*/ ? 0 : 1 }

    queue1.subscribe { println it }
