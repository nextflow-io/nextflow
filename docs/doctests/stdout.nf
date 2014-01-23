    process echoSomething {
        output: 
        stdout channel
    
        "echo Hello world!"
      
    }
    
    channel.subscribe { print "I say..  $it" }
