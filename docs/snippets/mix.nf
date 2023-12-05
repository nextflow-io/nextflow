c1 = Channel.of( 1, 2, 3 )
c2 = Channel.of( 'a', 'b' )
c3 = Channel.of( 'z' )

c1.mix(c2, c3).view()