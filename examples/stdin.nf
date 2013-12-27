str = Channel.from('hello', 'hola', 'bonjour', 'ciao').map { it+'\n' }

process printAll {
   input: stdin str

   """
   cat -
   """
}	