alpha: $bash_var
delta: $(bash expr)
gamma: ${bash var}
long \
     | line
groovy1: !foo
groovy2: !obj.beta
groovy3: !{obj.pico}
groovy4: !{(x + y) / z}
groovy5: '!{(x + z) / y}'
groovy6: "!{x + y / z}"
do not resolve escape: $1\t$2\t$3\n
ignore this: !&
.. and this: ! hola
.. and this: !\t