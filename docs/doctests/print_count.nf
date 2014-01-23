    process printCount {

      input:
      val cheers from 'Bonjour', 'Ciao', 'Hello', 'Hola'

      share:
      val count from 1 into result

      script:
      count += count
      "echo $cheers world! "

    }


    result.subscribe  { println "Process result: $it" }
