 methods = ['prot','dna', 'rna']

   process anyValue {
     input:
     val x from methods

     output:
     val x into receiver

     "echo $x > file"

   }

   receiver.subscribe { println "Received: $it" }
