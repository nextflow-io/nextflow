left  = channel.of( ['X', 1], ['Y', 2], ['Z', 3], ['P', 7] )
right = channel.of( ['Z', 6], ['Y', 5], ['X', 4] )

left.join(right, remainder: true).view()