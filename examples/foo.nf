process sayHello {
        input:
        val intValue
        val stringValue
        val listValue
        val mapValue
        // path 'foo.txt'
        // tuple val(intValue), val(stringValue)

        output:
        stdout
        val intValue
        val stringValue
        val listValue
        val mapValue
        path 'foo.txt'
        tuple val(intValue), val(stringValue), val(listValue), val(mapValue), path('foo.txt')

        """
        #!python3
        a = "$stringValue" + str($intValue)
        print(a)

        with open("foo.txt", "w") as f:
            f.write(str(a))
        """
}
process testBash {
        input:
        val x

        output:
        stdout

        """
        echo $x
        """
}

process testBash2 {
        input:
        val x

        output:
        stdout

        """
        echo $x
        """
}

workflow {

    Channel
            .of(1, 2, 3)
            .set { a }

    Channel
            .of(4, 5, 6)
            .set { b }

    (ch_stdout, ch_x, ch_y, ch_path, ch_tuple) = sayHello(testBash(a), testBash2(b), [1, 2, 3], ["john": "smith"])

    ch_x
    .map({it.toInteger() * it.toInteger()})
    .branch {
            small: it <= 20
            large: it > 20
    }
    .small.view({it})

    ch_stdout
    .map({it.toInteger() * it.toInteger()})
    .branch {
            small: it <= 20
            large: it > 20
    }
    .small.view({it})

   numbers = Channel.of(2, 2, 3)
   words = Channel.of('hello', 'ciao')
   numbers
           .combine(words)
           .view()

}
