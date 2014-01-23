    process saveHello {
      input:
      val cheers from 'Bonjour', 'Ciao', 'Hello', 'Hola'

      share:
      file greetings into result
      
      "echo '$cheers' >> $greetings "

    }


    result.subscribe { println it.text  }  


