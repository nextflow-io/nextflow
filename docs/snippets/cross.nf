source = Channel.of( [1, 'alpha'], [2, 'beta'] )
target = Channel.of( [1, 'x'], [1, 'y'], [1, 'z'], [2, 'p'], [2, 'q'], [2, 't'] )

source.cross(target).view()