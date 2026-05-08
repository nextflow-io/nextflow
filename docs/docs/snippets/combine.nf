numbers = channel.of(1, 2, 3)
words = channel.of('hello', 'ciao')

numbers
    .combine(words)
    .view()