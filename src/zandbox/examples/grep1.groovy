import nextflow.Channel

Channel
        .from( 'a', 'b', 'aa', 'bc', 3, 4.5 )
        .grep( ~/a+/ )
        .subscribe{ println it  }

sleep 100