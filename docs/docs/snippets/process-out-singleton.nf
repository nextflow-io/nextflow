process echo {
  input:
  val greeting

  output:
  val greeting

  exec:
  true
}

process greet {
  input:
  val greeting
  val name

  output:
  val "$greeting, $name!"

  exec:
  true
}

workflow {
  names = channel.of( 'World', 'Mundo', 'Welt' )
  greeting = echo('Hello')
  result = greet(greeting, names)
  result.view()
}
