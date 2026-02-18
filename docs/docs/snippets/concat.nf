a = channel.of( 'a', 'b', 'c' )
b = channel.of( 1, 2, 3 )
c = channel.of( 'p', 'q' )

c.concat( b, a ).view()