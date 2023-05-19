a = Channel.of( 'a', 'b', 'c' )
b = Channel.of( 1, 2, 3 )
c = Channel.of( 'p', 'q' )

c.concat( b, a ).view()