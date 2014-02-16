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

Transforming operators are used to transform the items emitted by a channel to new values.

These operators are:

* `buffer`_
* `collate`_
* `flatten`_
* `flatMap`_
* `groupBy`_
* `map`_
* `reduce`_
* `toList`_
* `toSortedList`_



map
------

The ``map`` operator applies a function of your choosing to every item emitted by a channel, and 
returns the items so obtained as a new channel. The function applied is called the `mapping` function 
and is expressed with a :ref:`closure <script-closure>` as shown in the example below::

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



flatMap
----------

The ``flatMap`` operator applies a function of your choosing to every item emitted by a channel, and
returns the items so obtained as a new channel. Whenever the `mapping` function returns a list of items,
this list is flattened so that each single item is emitted on its own.  

For example::

    // create a channel of numbers
    numbers = Channel.from( 1, 2, 3 )

    // map each number to a tuple (array), which items are emitted separately
    results = numbers.flatMap { n -> [ n*2, n*3 ] }

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


Associative arrays are handled in the same way, so that each array entry is emitted as a single `key-value` item. For example::

    Channel.from ( 1, 2, 3 )
           .flatMap { it -> [ number: it, square: it*it ] }
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

The ``reduce`` operator applies a function of your choosing to every item emitted by a channel.
Each time this function is invoked it takes two parameters: firstly the `i-th` emitted item
and secondly the result of the previous invocation of the function itself. The result is 
passed on to the next function call, along with the `i+1 th` item, until all the items are 
processed.

Finally, the ``reduce`` operator emits the result of the last invocation of your function 
as the sole output.

For example::

    Channel
        .from( 1, 2, 3, 4, 5 )
        .reduce { a, b -> println "a: $a b: $b"; return a+b }
        .subscribe { println "result = $it" }


It prints the following output::

	a: 1	b: 2
	a: 3	b: 3
	a: 6	b: 4
	a: 10	b: 5
	result = 15


.. note:: In a common usage scenario the first function parameter is used as an `accumulator` and
  the second parameter represents the `i-th` item to be processed.

Optionally you can specify a `seed` value in order to initialise the accumulator parameter
as shown below::

    myChannel.reduce( seedValue ) {  a, b -> ... }



groupBy
----------

The ``groupBy`` operator collects the values emitted by the source channel grouping them together using a `mapping`
function that associates each item with a key. When finished, it emits an associative
array that maps each key to the set of items identified by that key.  

For example::

    Channel
    	.from('hello','ciao','hola', 'hi', 'bonjour')
    	.groupBy { String str -> str[0] } 
    	.subscribe { println it }

:: 

    [ b:['bonjour'], c:['ciao'], h:['hello','hola','hi'] ]
    

The `mapping` function is an optional parameter. When omitted the values are grouped 
following these rules: 

* Any value of type ``Map`` is associated with the value of its first entry, or ``null`` when the map itself is empty.
* Any value of type ``Map.Entry`` is associated with the value of its ``key`` attribute.
* Any value of type ``Collection`` or ``Array`` is associated with its first entry.
* For any other value, the value itself is used as a key.


buffer
---------

The ``buffer`` operator gathers the items emitted by the source channel into subsets and emits these subsets separately.


There are a number of ways you can regulate how ``buffer`` gathers the items from
the source channel into subsets:

* ``buffer( closingCondition )``: starts to collect the items emitted by the channel into 
  a subset until the `closing condition` is verified. After that the subset is emitted 
  to the resulting channel and new items are gathered into a new subset. The process is repeated 
  until the last value in the source channel is sent. The ``closingCondition`` can be specified 
  either as a :ref:`regular expression <script-regexp>`, a Java class, a literal value, or a `boolean predicate`
  that has to be satisfied. For example::
  
    Channel
        .from( 1,2,3,1,2,3 ) 
        .buffer { it == 2 } 
        .subscribe {  println it }

    // emitted values
    [1,2]
    [3,1,2]
  
  

* ``buffer( openingCondition, closingCondition )``: starts to gather the items emitted by the channel 
  as soon as one of the them verify the `opening condition` and it continues until there is one item
  which verify the `closing condition`. After that the subset is emitted and it continues applying the 
  described logic until the last channel item is emitted.
  Both conditions can be defined either as a :ref:`regular expression <script-regexp>`, a literal value,
  a Java class, or a `boolean predicate` that need to be satisfied. For example:: 
 
    Channel
        .from( 1,2,3,4,5,1,2,3,4,5,1,2 ) 
        .buffer( 2, 4 ) 
        .subscribe {  println it }

    // emits bundles starting with '2' and ending with'4'
    [2,3,4]
    [2,3,4]      
  

* ``buffer( size: n )``: transform the source channel in such a way that it emits tuples 
  made up of ``n`` elements. An incomplete tuple is discarded. For example::

    Channel
        .from( 1,2,3,1,2,3,1 ) 
        .buffer( size: 2 )
        .subscribe {  println it }
        
    // emitted values 
    [1, 2]
    [3, 1]
    [2, 3]

If you want to emit the last items in a tuple containing less than ``n`` elements, simply 
add the parameter ``remainder`` specifying ``true``, for example::

    Channel
        .from( 1,2,3,1,2,3,1 )
        .buffer( size: 2, remainder: true )
        .subscribe {  println it }

    // emitted values
    [1, 2]
    [3, 1]
    [2, 3]
    [1]



* ``buffer( size: n, skip: m )``: as in the previous example, it emits tuples containing ``n`` elements, 
  but skips `m` values before starting to collect the values for the next tuple (including the first emission). For example::

    Channel
        .from( 1,2,3,4,5,1,2,3,4,5,1,2 ) 
        .buffer( size:3, skip:2 )
        .subscribe {  println it }
        
    // emitted values 
    [3, 4, 5]
    [3, 4, 5]

If you want to emit the remaining items in a tuple containing less than ``n`` elements, simply
add the parameter ``remainder`` specifying ``true``, as shown in the previous example.

See also: `collate`_ operator.


collate
---------

The ``collate`` operator transforms a channel in such a way that the emitted values are grouped in tuples containing `n` items. For example::

    Channel
        .from(1,2,3,1,2,3,1)
        .collate( 3 )
        .subscribe { println it }

::

        [1, 2, 3]
        [1, 2, 3]
        [1]

As shown in the above example the last tuple may be incomplete e.g. contain less elements than the specified size.
If you want to avoid this, specify ``false`` as the second parameter. For example::

    Channel
        .from(1,2,3,1,2,3,1)
        .collate( 3, false )
        .subscribe { println it }

::

        [1, 2, 3]
        [1, 2, 3]


A second version of the ``collate`` operator allows you to specify, after the `size`, the `step` by which elements
are collected in tuples. For example::

    Channel
      .from(1,2,3,4)
      .collate( 3, 1 )
      .subscribe { println it }

::

    [1, 2, 3]
    [2, 3, 4]
    [3, 4]
    [4]

As before, if you don't want to emit the last items which do not complete a tuple, specify ``false`` as the third parameter.


See also: `buffer`_ operator.


flatten
----------

The ``flatten`` operator transforms a channel in such a way that every item of type ``Collection`` or ``Array``
is flattened so that each single entry is emitted separately by the resulting channel. For example::

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
    
    
See also: `flatMap`_ operator.



toList
---------

The ``toList`` operator collects all the items emitted by a channel to a ``List`` object
and emits the resulting collection as a single item. For example::

    Channel
    	.from( 1, 2, 3, 4 )
    	.toList() 
    	.subscribe onNext: { println it }, onComplete: 'Done'
    	
::
 
    [1,2,3,4]
    Done


toSortedList
---------------


The ``toSortedList`` operator collects all the items emitted by a channel to a ``List`` object where they are sorted
and emits the resulting collection as a single item. For example::

    Channel
    	.from( 3, 2, 1, 4 )
    	.toSortedList()
    	.subscribe onNext: { println it }, onComplete: 'Done'

::

    [1,2,3,4]
    Done


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

The ``into`` operator connects two channels in such a way the values emitted by the source channel
are forwarded to the target channel. For example::

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

The ``tap`` operator combines the functions of `into`_ and `split`_ operators in such a way that
it connects two channels, copying the values from the source into the `tapped` channel. At the same
time it splits the source channel into a newly created channel that is returned by the operator itself.

The ``tap`` can be useful in certain scenarios where you may be required to concatenate multiple operations,
as in the following example::


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

The ``merge`` operator lets you join items emitted by two (or more) channels into a new channel.

For example the following code merges two channels together, one which emits a series of odd integers
and the other which emits a series of even integers::

    odds  = Channel.from([1, 3, 5, 7, 9]);
    evens = Channel.from([2, 4, 6]);

    odds
        .merge( evens ) { o, e -> [o, e] }
        .subscribe { println it }

::

    [1, 2]
    [3, 4]
    [5, 6]

If you want to merge more than to channels, they must be specified as a list i.e. by wrapping 
them in square brackets, as shown in the following example::

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
  regardless of which source channel they came from. Thus, the following example
  it could be a possibile result of the above example as well.

::

          'z'
          1
          'a'
          2
          'b'
          3


phase
--------

The ``phase`` operator creates a channel that synchronizes the values emitted by two other channels,
in such a way that it emits pairs of items that have a matching key.

The key is defined, by default, as the first entry in an array, a list or map object,
or the value itself for any other data type.

For example::

        ch1 = Channel.from( 1,2,3 )
        ch2 = Channel.from( 1,0,0,2,7,8,9,3 )

        result = ch1.phase(ch2).subscribe { println it }

::

    [1,1]
    [2,2]
    [3,3]


Optionally, a mapping function can be specified in order to provide a custom rule to associate an item to a key,
as shown in the following example::


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

The ``cross`` operators allows you to combine the items of two channes in such a way that 
the items of the source channel are emitted along with the items emitted by the target channel 
for which they have a matching key.  

The key is defined, by default, as the first entry in an array, a list or map object,
or the value itself for any other data type. For example:: 

	source = Channel.from( [1, 'alpha'], [2, 'beta'] )
	target = Channel.from( [1, 'x'], [1, 'y'], [1, 'z'], [2,'p'], [2,'q'], [2,'t'] )

	source.cross(target).subscribe { println it }

It will output:: 

	[ [1, alpha], [1, x] ]
	[ [1, alpha], [1, y] ]
	[ [1, alpha], [1, z] ]
	[ [2, beta],  [2, p] ]
	[ [2, beta],  [2, q] ]
	[ [2, beta],  [2, t] ]

The above example shows how the items emitted by the source channels are associated to the ones
emitted by the target channel (on the right) having the same key. 

There are two important cavets when using the ``cross`` operator: 

	#. The operator is not `reflexive`, i.e. the result of ``a.cross(b)`` is different from ``b.cross(a)`` 
	#. The source channel should emits items for which there's no key repetition i.e. the emitted 
	   items have an unique key identifier. 

Optionally, a mapping function can be specified in order to provide a custom rule to associate an item to a key,
in a similar manner as shown for the `phase`_ operator.

spread
---------

The ``spread`` operator combines the items emitted by the source channel with all the values in an array
or a ``Collection`` object specified as the operator argument. For example::

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

* `choice`_
* `separate`_
* `split`_
* `route`_


split
--------

The ``split`` operator copies the items emitted by a source channel into one or more channels specified as an argument.
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


A second version of the ``split`` operator takes an integer `n` as an argument and returns 
a list of `n` channels, each of which emits a copy of the items there were emitted by the 
source channel. For example::


    (split1,split2) = Channel.from( 'a','b','c').split(2)
    split1.subscribe { println "Channel 1 emit: " + it }
    split2.subscribe { println "Channel 2 emit: " + it }

::

    Channel 1 emit: a
    Channel 1 emit: b
    Channel 1 emit: c

    Channel 2 emit: a
    Channel 2 emit: b
    Channel 2 emit: c


.. note:: The above example takes advantage of the :ref:`multiple assignment <script-multiple-assignment>` syntax
  in order to assign two variables at once using the list of channels returned by the ``split`` operator.



See also `tap`_ operator.


choice
----------

The ``choice`` operator allows you to forward the items emitted by a source channel to two 
(or more) output channels, `choosing` one out of them at a time. 

The destination channel is selected by using a :ref:`closure <script-closure>` that must return the `index` number of the channel
where the item has to be sent. The first channel is identified by the index ``0``, the second as ``1`` and so on. 

The following example sends all string items beginning with ``Hello`` into ``queue1``, 
the others into ``queue2``  

::
  
    source = Channel.from 'Hello world', 'Hola', 'Hello John'
    queue1 = Channel.create()
    queue2 = Channel.create()

    source.choice( queue1, queue2 ) { a -> a =~ /^Hello.*/ ? 0 : 1 }

    queue1.subscribe { println it }



If you want to use more than two channels in the ``choice`` operator, they must be specified as a list 
i.e. by wrapping them in square brackets, as shown in the following example::

    source.choice( [queue1, queue2, queue3] ) { a -> a % 3 }
     

separate
------------

The ``separate`` operator lets you copy the items emitted by the source channel into multiple 
channels, which each of these can receive a `separate` version of the same item. 

The operator applies a function of your choosing to every item emitted by the source channel. 
This function must return a list of as many values as there are output channels. Each entry in the result 
list will be assigned to the output channel with the corresponding position index. For example:: 

    queue1 = Channel.create()
    queue2 = Channel.create()

    Channel
        .from ( 2,4,8 ) 
        .separate( queue1, queue2 ) { a -> [a+1, a*a] }

    queue1.subscribe { println "Channel 1: $it" }
    queue2.subscribe { println "Channel 2: $it" }
	
::

	Channel 1: 3
	Channel 2: 4
	Channel 1: 5
	Channel 2: 16
	Channel 2: 64
	Channel 1: 9

.. note:: If you want to use more than two channels in the ``separate`` operator, they must be specified 
	as a list  i.e. by wrapping them in square brackets.
	

A second version of the ``separate`` operator takes an integer `n` as an argument and returns a list of `n` channels, 
each of which gets a value from the corresponding element in the list returned by the closure as explained above.
For example::	

  source = Channel.from(1,2,3)
  (queue1, queue2, queue3) = source.separate(3) { a -> [a, a+, a*a] }

.. note:: The above example takes advantage of the :ref:`multiple assignment <script-multiple-assignment>` syntax
  in order to assign two variables at once using the list of channels returned by the ``separate`` operator.

See also: `choice`_, `split`_ and `map`_ operators.


route
----------

The ``route`` operator allows you to forward the items emitted by the source channel 
to a channel which is associated with the item's key. 

The channel's keys are specified by using a map parameter as the operator argument, 
that associates each channel with a key identifier. 

The item's key is defined, by default, as the first entry in an array, a list or map object,
or the value itself for any other data type.

Optionally, a mapping function can be specified as a parameter in order to provide a custom 
rule to associate an item with a key, as shown in the example below::

For example::

    r1 = Channel.create()
    r2 = Channel.create()
    r3 = Channel.create()
	
    Channel
        .from('hello','ciao','hola', 'hi', 'x', 'bonjour')
        .route ( b: r1, c: r2, h: r3 ) { it[0] }

	r3.subscribe { println it }        

:: 

	hello
	hola
	hi

In the above example all the string items starting with the letter ``b`` are copied to the 
channel ``r1``, the items that begin with ``c`` to the channel ``r2`` and the ones beginning
with ``h`` are copied to the channel ``r3``. Other items eventually existing are discarded. 
 

See also: `into`_, `choice`_  and `separate`_ operators.


Maths operators
================

This section talks about operators that performs maths operations on channels.

The maths operators are:

* `count`_
* `countBy`_
* `min`_
* `max`_
* `sum`_


count
--------

The ``count`` operator creates a channel that emits a single item: a number that represents the total number of
items emitted by the source channel. For example:: 

        Channel
            .from(9,1,7,5)
            .count()
            .subscribe { println it }
        // -> 4


An optional parameter can be provided in order to select which items are to be counted. 
The selection criteria can be specified either as a :ref:`regular expression <script-regexp>`, 
a literal value, a Java class, or a `boolean predicate` that needs to be satisfied. For example::


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
        
        Channel
            .from('a','c','c','q','b')
            .count { it <= 'c' }
            .subscribe { println it }
        // -> 4


countBy
----------

The ``countBy`` operator creates a channel which emits an associative array (i.e. ``Map`` object) 
that counts the occurrences of the emitted items in the source channel having the same key. 
For example::

    Channel
        .from( 'x', 'y', 'x', 'x', 'z', 'y' )
        .countBy()
        .subscribe { println it }

::

    [x:3, y:2, z:1]


An optional grouping criteria can be specified by using a :ref:`closure <script-closure>` 
that associates each item with the grouping key. For example::


    Channel
        .from( 'hola', 'hello', 'ciao', 'bonjour', 'halo' )
        .countBy { it[0] }
        .subscribe { println it }


::

    [h:3, c:1, b:1]



min
------

The ``min`` operator waits until the source channel completes, and then emits the item that has the lowest value.
For example::

    Channel
        .from( 8, 6, 2, 5 )
        .min()
        .subscribe { println "Min value is $it" }

::

  Min value is 2

An optional :ref:`closure <script-closure>` parameter can be specified in order to provide 
a function that returns the value to be compared. The example below shows how to find the string 
item that has the minimum length:: 

    Channel
    	.from("hello","hi","hey")
    	.min { it.size() } 
    	.subscribe {  println it }

::

	 "hi"

Alternatively it is possibile to specify a comparator function i.e. a :ref:`closure <script-closure>` 
taking two parameters that represent two emitted items to be compared. For example:: 


    Channel
    	.from("hello","hi","hey")
    	.min { a,b -> a.size() <=> b.size() } 
    	.subscribe {  println it }



max
------

The ``max`` operator waits until the source channel completes, and then emits the item that has the greatest value.
For example::

    Channel
        .from( 8, 6, 2, 5 )
        .min()
        .subscribe { println "Max value is $it" }

::

  Max value is 8


An optional :ref:`closure <script-closure>` parameter can be specified in order to provide 
a function that returns the value to be compared. The example below shows how to find the string 
item that has the maximum length:: 

    Channel
    	.from("hello","hi","hey")
    	.max { it.size() } 
    	.subscribe {  println it }

::

	 "hello"

Alternatively it is possibile to specify a comparator function i.e. a :ref:`closure <script-closure>` 
taking two parameters that represent two emitted items to be compared. For example:: 


    Channel
    	.from("hello","hi","hey")
    	.max { a,b -> a.size() <=> b.size() } 
    	.subscribe {  println it }



sum
------

The ``sum`` operator creates a channel that emits the sum of all the items emitted by the channel itself.
For example::

    Channel
        .from( 8, 6, 2, 5 )
        .sum()
        .subscribe { println "The sum is $it" }

::

    The sum is 21


An optional :ref:`closure <script-closure>` parameter can be specified in order to provide 
a function that, given an item, returns the value to be summed. For example:: 

	Channel
		.from( 4, 1, 7, 5 )
		.sum { it * it } 
		.subscribe {  println "Square: $it" } 

::

	Square: 91




