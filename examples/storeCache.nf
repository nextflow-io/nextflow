
cheers = Channel.from 'Hello', 'Ciao', 'Hola'


process storeCache  {
    storeDir 'cache'

    input:
    val cheers

    output:
    file "${cheers}.txt" into salut

    "printf $cheers > ${cheers}.txt"

}

salut.subscribe { println it }