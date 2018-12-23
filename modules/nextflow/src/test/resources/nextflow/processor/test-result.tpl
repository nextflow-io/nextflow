alpha: $bash_var
delta: $(bash expr)
gamma: ${bash var}
long \
     | line
groovy1: !foo
groovy2: !obj.beta
groovy3: World
groovy4: 5
groovy5: '11'
groovy6: "9.5"
do not resolve escape: $1\t$2\t$3\n
ignore this: !&
.. and this: ! hola
.. and this: !\t