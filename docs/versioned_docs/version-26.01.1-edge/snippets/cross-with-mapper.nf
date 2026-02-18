source = channel.of( [1, 'alpha'], [2, 'beta'] )
target = channel.of( [1, 'a'], [1, 'b'], [2, 'a'], [2, 'b'] )

source .cross(target) { v -> v[1][0] } .view()