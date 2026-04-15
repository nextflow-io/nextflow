left  = channel.of( record(id: 'X', a: 1), record(id: 'X', a: 3) )
right = channel.of( record(id: 'X', b: 2), record(id: 'X', b: 4) )

left.join(right, by: 'id').view()