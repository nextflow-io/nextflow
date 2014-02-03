.. _operator-page:

*******************
Operators
*******************

Dataflow `operators` are methods that allow you to connect channels to each other or to transform values
emitted by a channel applying some user provided rules.

Operators can be separated in to five groups:

* `Filtering operators`_
* `Transforming operators`_
* `Combining operators`_
* `Forking operators`_
* `Maths operators`_



Filtering operators
===================

Given a channel, filtering operators allow you to select only the items that comply with a given rule.

The available filter operators are:

* `distinct`_
* `filter`_
* `first`_
* `grep`_
* `last`_
* `take`_
* `unique`_

filter
---------

The ``filter`` operator allows you to get only the items emitted by a channel, that satisfy a boolean `predicate`. For example::

    Channel
        .from( 1, 2, 3, 4, 5 )
        .filter { it % 2 == 1 }
        .subscribe { println it }

::

    1
    3
    5

.. note:: Using Nextflow a predicate is represented with a :ref:`closure <script-closure>` retuning a boolean value

grep
-------

The ``grep`` operator allows you to `filter` a channel by discarding all the items that do not meet a specified condition.
It works in a similar manner to the `filter`_ operator but conditions can be specified in a more flexible manner. A grep
condition can be either a :ref:`regular expression <script-regexp>`, a Java `class`, a literal value or a boolean predicate.

The following example shows how to filter a channel by using a regular expression::

    Channel
        .from( 'a', 'b', 'aa', 'bc', 3, 4.5 )
        .grep( ~/a+/ )
        .subscribe { println it }

::

    a
    aa


The following example shows how to filter a channel by specifying the class ``Number`` so that only numbers are returned::

    Channel
        .from( 'a', 'b', 'aa', 'bc', 3, 4.5 )
        .grep( Number )
        .subscribe { println it }

::

    3
    4.5


The following example shows how to filter a channel by using a `boolean` expression::

     Channel
         .from( 'a', 'b', 'aa', 'bc', 3, 4.5 )
         .grep { it.toString().size() == 1 }
         .subscribe { println it }

::

     a
     aa
     3

.. tip:: In the above example the filtering expression is wrapped in curly brackets instead of normal
  round brackets, since it specifies a :ref:`closure <script-closure>` as the operator's argument.
  This just is an shorter syntax for ``grep({ it.toString().size() == 1 })``



unique
---------

The ``unique`` operator allows you to remove duplicate items from a channel and only emit single items with no repetition.

For example::

    Channel
        .from( 1,1,1,5,7,7,7,3,3 )
        .unique()
        .subscribe { println it }

::

    1
    5
    7
    3


You can also specify an optional :ref:`closure <script-closure>` that customizes the way it distinguishes between unique items.
For example::

    Channel
        .from(1,3,4,5)
        .unique { it % 2 }
        .subscribe { println it }

::

    1
    4


distinct
-----------

The ``distinct`` operator allows you to remove `consecutive` duplicated items from a channel, so that each emitted item
is different from the preceding one. For example::


    Channel
        .from( 1,1,2,2,2,3,1,1,2,2,3 )
        .distinct()
        .subscribe onNext: { println it }, onComplete: { println 'Done' }

::

    1
    2
    3
    1
    2
    3
    Done



You can also specify an optional :ref:`closure <script-closure>` that customizes the way it distinguishes between distinct items.
For example::

    Channel
        .from( 1,1,2,2,2,3,1,1,2,4,6 )
        .distinct { it % 2 }
        .subscribe onNext: { println it }, onComplete: { println 'Done' }


::

    1
    2
    3
    2
    Done



first
--------

The ``first`` operator creates a channel that returns the first item emitted by the source channel, or eventually
the first item that matches an optional condition. The condition can be specified by using a :ref:`regular expression<script-regexp>`,
a Java `class` type or any boolean `predicate`. For example::


    // no condition is specified, emits the very first item: 1
    Channel
        .from( 1, 2, 3 )
        .first()
        .subscribe { println it }


    // emits the first String value: 'a'
    Channel
        .from( 1, 2, 'a', 'b', 3 )
        .first( String )
        .subscribe { println it }

    // emits the first item matching the regular expression: 'aa'
    Channel
        .from( 'a', 'aa', 'aaa' )
        .first( ~/aa.*/ )
        .subscribe { println it }

    // emits the first item for which the predicate evaluates to true: 4
    Channel
        .from( 1,2,3,4,5 )
        .first { it > 3 }
        .subscribe { println it }


take
-------

The ``take`` operator allows you to filter only the first `n` items emitted by a channel. For example::

    Channel
        .from( 1,2,3,4,5,6 )
        .take( 3 )
        .subscribe onNext: { println it }, onComplete: { println 'Done' }

::

    1
    2
    3
    Done


last
-------

The ``last`` operator creates a channel that only returns the last item emitted by the source channel. For example::

    Channel
        .from( 1,2,3,4,5,6 )
        .last()
        .subscribe { println it }

::

    6


Transforming operators
======================

This section explains operators with which you can transform items that are emitted by channels.

The available operators are:

* `map`_
* `mapMany`_
* `reduce`_
* `groupBy`_
* `buffer`_
* `collate`_
* `flatten`_
* `toList`_
* `toSortedList`_



map
------

Transform the items emitted by a channel by applying a function to each of them. For example::

    Channel
        .from( 1, 2, 3, 4, ,5 )
        .map { it * it  }
        .subscribe onNext: { println it }, onComplete: { println 'Done' }

::

    1
    4
    9
    16
    25
    Done



mapMany
----------

Transform the items emitted by a channel by applying a function to each of them and then flattens the results of that function.


For example::

    // create a channel of numbers
    numbers = Channel.from( 1, 2, 3 )

    // map each number to a tuple (array), which items are emitted separately
    results = numbers.mapMany { n -> [ n*2, n*3 ] }

    // print the final results
    results.subscribe onNext: { println it }, onComplete: { println 'Done' }

::

    2
    3
    4
    6
    6
    9
    Done


Associative arrays are manged in the same way. Fo example::

    Channel.from ( 1, 2, 3 )
           .mapMany { it -> [ number: it, square: it*it ] }
           .subscribe { println it.key + ': ' + it.value }

::

    number: 1
    square: 1
    number: 2
    square: 4
    number: 3
    square: 9


reduce
---------

The ``reduce`` operator applies a function of your choosing to the first item emitted by a source channel,
then feeds the result of that function along with the second item emitted by the source channel into the same function,
then feeds the result of that function along with the third item into the same function, and so on until all items have been emitted by the source channel.

Finally it emits the final result from the final call to your function as the sole output from the returned channel.

For example::

    Channel
        .from( 1, 2, 3, 4, 5 )
        .reduce { a, b ->  a+b }
        .subscribe onNext: { println it }, onComplete: { println 'Done' }


::

    15
    Done

There is also a version of ``reduce`` to which you can pass a `seed` item in addition to an accumulator function::

    myChannel.reduce( seedValue ) {  a, b -> ... }


The seed value is used to initialize the accumulator argument i.e. the ``a`` argument in the above example.

groupBy
----------

The ``groupBy`` operator collects all the values emitted by the source channel grouping them by a key defined by
a user defined mapping function. When finished emits an associative array that maps each key with 
the list of emitted values having that key.

For example::

    Channel
    	.from('hello','ciao','hola', 'hi', 'bonjour')
    	.groupBy { String str -> str[0] } 
    	.subscribe { println it }

:: 

    [ b:['bonjour'], c:['ciao'], h:['hello','hola','hi'] ]
    

When the mapping function is omitted, the grouping key is defined as following: 

* Any value of type ``Map`` is associated to the value of its first entry, or ``null`` when the map itself is empty.
* Any value of type ``Map.Entry`` is associated to the value of its ``key`` attribute.
* Any value of type ``Collection`` or ``Array`` is associated to its first entry.
* Any other value, the value itself is used as key.


buffer
---------

The ``buffer`` operator gathers the items emitted by the source channel into bundles and
and emits these bundles as its own emissions. 


There are a number of ways with which you can regulate how ``buffer`` gathers items from
the source channel into bundles:

* ``buffer( closingCondition )``: starts to collect the values emitted by the channel into 
  a bundle until the ``closingCondition`` condition is verified, after that the bundle is emitted 
  to the returned channel and new values are gathered into a new bundle. The process is repeated until 
  the last value in the source channel is sent. The ``closingCondition`` can be specified as a 
  `regular expression`, a type or a value that have to be matched or any closure evaluating a 
  boolean expression to be satisfied. For example:: 
  
    Channel
        .from( 1,2,3,1,2,3 ) 
        .buffer { it == 2 } 
        .subscribe {  println it }

    // emitted values
    [1,2]
    [3,1,2]
  
  

* ``buffer( openingCondition, closingCondition )``: starts to gather the value emitted by the channel 
  when the ``openingCondition`` is verified and on continue to buffer new emitted values until one of 
  them verify the ``closingCondition``. After that the bundle is emitted and it continues applying the 
  describe logic until the last value is emitted by the source channel. 
  Both condition can be defined as a `regular expression`, a value literal, a value type or any closure
  evaluating to a boolean expression. For example:: 
 
    Channel
        .from( 1,2,3,4,5,1,2,3,4,5,1,2 ) 
        .buffer( 2, 4 ) 
        .subscribe {  println it }

    // emits bundles starting with '2' and ending with'4'
    [2,3,4]
    [2,3,4]      
  

* ``buffer( size: n )``: emits tuple made up of exactly ``n`` values. Incomplete tuple
  are discarded. For example::

    Channel
        .from( 1,2,3,1,2,3,1 ) 
        .buffer( size: 2 )
        .subscribe {  println it }
        
    // emitted values 
    [1, 2]
    [3, 1]
    [2, 3]

If you want to emit also the last, eventually incomplete tuple, add the parameter ``remainder`` specifying ``true``, for example::

    Channel
        .from( 1,2,3,1,2,3,1 )
        .buffer( size: 2, remainder: true )
        .subscribe {  println it }

    // emitted values
    [1, 2]
    [3, 1]
    [2, 3]
    [1]



* ``buffer( size: n, skip: m )``: like the previous emits tuple containing exactly ``n`` values,
  but skip `m` values before to start collecting values for the next tuple (including the first).
  For example::

    Channel
        .from( 1,2,3,4,5,1,2,3,4,5,1,2 ) 
        .buffer( size:3, skip:2 )
        .subscribe {  println it }
        
    // emitted values 
    [3, 4, 5]
    [3, 4, 5]

As before, incomplete tuple are not emitted. If you need to emit them add the parameter ``remainder``
specifying ``true`` as showed in the previous example.

See also: `collate`_ operator.


collate
---------

The ``collate`` operator transforms a channel in a such way that the emitted values are grouped in tuple containing `n` items. For example::

    Channel
        .from(1,2,3,1,2,3,1)
        .collate( 3 )
        .subscribe { println it }

::

        [1, 2, 3]
        [1, 2, 3]
        [1]

As shown in the above example the last tuple may be incomplete e.g. contain less elements than the specified size.
If you need toi avoid this specify ``false`` as second parameter. For example::

    Channel
        .from(1,2,3,1,2,3,1)
        .collate( 3, false )
        .subscribe { println it }

::

        [1, 2, 3]
        [1, 2, 3]


A second version of the ``collate`` operator allows to specify, after the `size`, the `step` by which elements are collected in tuple.
For example::

    Channel
      .from(1,2,3,4)
      .collate( 3, 1 )
      .subscribe { println it }

::

    [1, 2, 3]
    [2, 3, 4]
    [3, 4]
    [4]

Like before, if you need to emit only complete tuple, specify ``false`` as third parameter.


See also: `buffer`_ operator.


flatten
----------

The ``flatten`` operator transforms the items emitted by a channel so that any ``Collection`` or ``Array`` object
value is converted to its single items and emitted separately. For example:: 

    Channel
    	.from( [1,[2,3]], 4, [5,[6]] )
    	.flatten()
    	.subscribe { println it }

:: 
    
    1
    2
    3
    4
    5
    6
    
    
See also: `mapMany`_ operator.



toList
---------

The ``toList`` operator collects all the values emitted by a channel to a ``List`` object
and emits the resulting collection to the channel returned by the operator itself. For example::

    Channel
    	.from( 1, 2, 3, 4 )
    	.toList() 
    	.subscribe onNext: {  println it }, onComplete: 'Done'
    	
::
 
    [1,2,3,4]
    Done


toSortedList
---------------

TODO


Combining operators
===================

The combining operators are:

* `cross`_
* `into`_
* `merge`_
* `mix`_
* `phase`_
* `spread`_
* `tap`_


into
-------

The ``into`` operator connects two channels so that values emitted by the source channel
are forwarded to target channel. For example:: 

    target = Channel.create()
    target.subscribe { println it }
	
    Channel
        .from( 'a', 'b', 'c', 'd' )     
        .into( target )
      
::
  
    a
    b
    c
    d
    

See also `tap`_ and `route`_ operators.


tap
------

The ``tap`` operator, combines the functions of `into`_ and `split`_ operators, is a such way that
it connects two channels copying the values from the source into the `tapped` channel, at the same time 
it splits the source channel into a newly create channel that is returned by the operator itself. 

It may be convenient is some scenarios where may be required to concatenate multiple operation. 

For example::


    log1 = Channel.create().subscribe { println "Log 1: $it" }  
    log2 = Channel.create().subscribe { println "Log 2: $it" }
  
    Channel
        .from ( 'a', 'b', 'c' ) 
  	    .tap( log1 ) 
  	    .map { it * 2 }
  	    .tap( log2 ) 
  	    .subscribe { println "Result: $it" }

:: 

    Log 1: a
    Log 1: b
    Log 1: c

    Log 2: aa
    Log 2: bb
    Log 2: cc

    Result: aa
    Result: bb
    Result: cc



See also `into`_ and `split`_ operators


merge
--------

The ``merge`` operator let to join the items emitted by two (or more) channels into a new one.

For example, the following code merges together two channels, one of which emits a series of odd integers and the other
of which emits a series of even integers::

    odds  = Channel.from([1, 3, 5, 7, 9]);
    evens = Channel.from([2, 4, 6]);

    odds
        .merge( evens ) { o, e -> [o, e] }
        .subscribe { println it }

::

    [1, 2]
    [3, 4]
    [5, 6]


When merging more than two channels specify them by wrapping them with a list, as showed in the following example::

   channel1.merge( [channel2, channel3] ) {  item1, item2, item3 -> item1 + item2 + item3 }


mix
------

The ``mix`` operator combines the items emitted by two (or more) channels into a single channel.


For example::

        c1 = Channel.from( 1,2,3 )
        c2 = Channel.from( 'a','b' )
        c3 = Channel.form( 'z' )

        c1 .mix(c2,c3)
           .subscribe onNext: { println it }, onComplete: { println 'Done' }

::

        1
        2
        3
        'a'
        'b'
        'z'

.. note:: The items emitted by the resulting mixed channel may appear in any order,
  regardless of which source channel they came from. So also the following example
  may be a valid result of the above example.

::

          'z'
          1
          'a'
          2
          'b'
          3


phase
--------

The ``phase`` operator creates a channel that synchronize the values emitted by two source channels,
so that emitted items have the same key.

The key is defined by default as the first entry in a array, list or map object,
or the object itself for any other data type.

For example::

        ch1 = Channel.from( 1,2,3 )
        ch2 = Channel.from( 1,0,0,2,7,8,9,3 )

        result = ch1.phase(ch2).subscribe { println it }

::

    [1,1]
    [2,2]
    [3,3]


An optional function can be provided in order to specify a custom key mapping strategy. For example::


    ch1 = Channel.from( [sequence: 'aaaaaa', key: 1], [sequence: 'bbbbbb', key: 2] )
    ch2 = Channel.from( [val: 'zzzz', id: 3], [val: 'xxxxx', id: 1], [val: 'yyyyy', id: 2])

    // provide a custom function in order to get the right id depending the map
    result = ch1.phase(ch2) { Map it ->

        if( it.containsKey('key') ) {
            return it.key
        }
        else if( it.containsKey('id') ) {
            return it.id
        }
        return null

    }


::

    [ [sequence: 'aaaaaa', key: 1], [val: 'xxxxx', id: 1] ]
    [ [sequence: 'bbbbbb', key: 2], [val: 'yyyyy', id: 2] ]



cross
-------

TODO


spread
---------

The ``spread`` operator combines the items emitted by the source channel with all the values defined by an array or an ``Iterable`` object specified as argument.

For example::

    Channel
        .from(1,2,3)
        .spread(['a','b'])
        .subscribe onNext: { println it }, onComplete: { println 'Done' }

::

    [1, 'a']
    [1, 'b']
    [2, 'a']
    [2, 'b']
    [3, 'a']
    [3, 'b']
    Done




Forking operators
=================

The forking operators are:

* `split`_
* `choice`_
* `separate`_
* `route`_


split
--------

The ``split`` operator copies the items emitted by a source channel into one or multiple channels specified as argument.
For example::

    source = Channel.from( 'a', 'b', 'c' )
    copy1 = Channel.create().subscribe { println "Channel 1 emit: " + it }
    copy2 = Channel.create().subscribe { println "Channel 2 emit: " + it }

    source.split( copy1, copy2 )

::

    Channel 1 emit: a
    Channel 1 emit: b
    Channel 1 emit: c

    Channel 2 emit: a
    Channel 2 emit: b
    Channel 2 emit: c


.. TODO split(n)


See also: `tap`_ operator.


choice
----------

TODO

::

    queue1.choice([queue2, queue3, queue4]) {a -> a % 3}



separate
------------

Separation is the opposite operation to merge. The supplied closure returns a list of values, each of which will be output into an output channel with the corresponding position index.

::

    queue1.separate([queue2, queue3, queue4]) {a -> [a-1, a, a+1]}

    def (queue2, queue3, queue4) = queue1.separate(3) {a -> [a-1, a, a+1]}




route
----------

The ``route`` operator forwards the items emitted by the source channel to a set of channels according
a the map specified as a parameter. For example::

    channel1 = Channel.create()

    Channel
        .from(1,3,2,1,1,2,1,4,2)


TODO

See also: `into`_ operator.


Maths operators
==============

This section explains operators that perform mathematical operations on the items emitted by channels.

The math operator are:

* `count`_
* `countBy`_
* `min`_
* `max`_
* `sum`_


count
--------

The ``count`` operator creates a channel that emits a single item: a number that represents the total number of
items emitted by the source channel.

An optional parameter can be specified representing the condition to be satisfied by the item to count. For example::

        Channel
            .from(9,1,7,5)
            .count()
            .subscribe { println it }
        // -> 4


        Channel
            .from(4,1,7,1,1)
            .count(1)
            .subscribe { println it }
         // -> 3

        Channel
            .from('a','c','c','q','b')
            .count ( ~/c/ )
            .subscribe { println it }
        // -> 2


countBy
----------

The ``countBy`` creates a channel which emits an associative array (map object) counting the occurence of the emitted
items in the source channel. For example::

    Channel
        .from( 'x', 'y', 'x', 'x', 'z', 'y' )
        .countBy()
        .subscribe { println it }

::

    [x:3, y:2, z:1]


A optional grouping criteria can be specified by a function (closure) mapping each item to a grouping key. For example::


    Channel
        .from( 'hola', 'hello', 'ciao', 'bonjour', 'halo' )
        .countBy { it[0] }
        .subscribe { println it }


::

    [h:3, c:1, b:1]



min
------

The ``min`` operator waits until the source channel completes, and then emits the value that had the lowest value.
For example::

    Channel
        .from( 8, 6, 2, 5 )
        .min()
        .subscribe { println "Min value is $it" }

::

  Min value is 2



max
------

The ``max`` operator waits until the source channel completes, and then emits the value that had the greatest value.
For example::

    Channel
        .from( 8, 6, 2, 5 )
        .min()
        .subscribe { println "Max value is $it" }

::

  Max value is 8



sum
------

The ``sum`` operators crates a channel that emits the sum of all values emitted by the source channel to which is applied.
For example::

    Channel
        .from( 8, 6, 2, 5 )
        .min()
        .subscribe { println "The sum is $it" }



::

    The sum is 21










