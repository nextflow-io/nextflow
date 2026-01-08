
def CMD =  """
    mkdir -p a/a b/b c/c
    touch a/1.txt
    touch b/1.txt
    touch c/1.txt
    touch a/a/2.txt
    touch b/b/2.txt
    touch c/c/2.txt
    """

process foo {
  output:
  file("a/*/*.txt")
  script:
  CMD
}

process bar {
  scratch true
  output:
  file("a/*/*.txt")
  script:
  CMD
}

workflow {
  foo()
  foo.out.view { "foo: $it" }
  bar()
  bar.out.view { "bar: $it" }
} 
