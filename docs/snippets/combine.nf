numbers = Channel.of(1, 2, 3)
words = Channel.of('hello', 'ciao')

numbers
    .combine(words)
    .view()