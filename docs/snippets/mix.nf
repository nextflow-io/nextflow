c1 = channel.of( '1', '2', '3' )
c2 = channel.of( 'a', 'b' )
v3 = channel.value( 'z' )

c1.mix(c2).mix(v3).view()