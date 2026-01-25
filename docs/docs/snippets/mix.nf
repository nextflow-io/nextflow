c1 = channel.of( 1, 2, 3 )
c2 = channel.of( 'a', 'b' )
c3 = channel.of( 'z' )

c1.mix(c2, c3).view()