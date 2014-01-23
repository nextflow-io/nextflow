    process randomNum {
    
       output:
       file 'random_out' into numbers

       '''
       echo $RANDOM > random_out
       '''

    }
    
   numbers.subscribe { println "Received: " + it.text }
