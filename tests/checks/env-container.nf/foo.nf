process foo {
  echo true
  input: 
  env BAR from Channel.value('Hello World!')
  '''
  env | sort 
  hello.sh
  '''
}
