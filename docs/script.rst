.. _pipeline-page:

*****************
Pipeline script
*****************


The Nextflow scripting language is an extension of the Groovy programming language whose syntax has been
specialized to ease the writing of computational pipelines in a declarative manner.

This means that Nextflow can execute any Groovy piece of code or use any library for the Java JVM platform.

For a detailed description of the Groovy programming language, references these links:

* `Groovy User Guide <http://groovy-lang.org/documentation.html>`_
* `Groovy Cheat sheet <http://refcardz.dzone.com/refcardz/groovy>`_
* `Groovy in Action <http://www.manning.com/koenig2/>`_


Below you can find a crash course in the most important language constructs used in the Nextflow scripting language.

.. warning:: Nextflow uses ``UTF-8`` as default file character encoding for source and application files. Make sure
  to use the ``UTF-8`` encoding when editing Nextflow scripts with your favourite text editor.

Language basics
==================


Hello world
------------

To print something is as easy as using the ``print`` or the ``println`` methods.
::

    println "Hello, World!"

The only difference between the two is that the ``println`` method appends implicitly a `new line` character
to the printed string.


Variables
----------

To define a variable, simply assign a value to it. For example::

    x = 1
    println x

    x = new java.util.Date()
    println x

    x = -3.1499392
    println x

    x = false
    println x

    x = "Hi"
    println x


Lists
------

A List object can be defined by placing the list items in square brackets. For example::

    myList = [1776, -1, 33, 99, 0, 928734928763]

You can access a given item in the list with square bracket notation (indexes start at 0)::

    println myList[0]

In order to get the length of the list use the ``size`` method::

    println myList.size()


Learn more about lists:

* `Groovy Lists tutorial <http://groovy-lang.org/groovy-dev-kit.html#Collections-Lists>`_
* `Groovy List SDK <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/List.html>`_
* `Java List SDK <http://docs.oracle.com/javase/7/docs/api/java/util/List.html>`_


Maps
-----

Maps are used to store `associative arrays` or `dictionaries`. They are unordered collections of heterogeneous, named data.
For example::

    scores = [ "Brett":100, "Pete":"Did not finish", "Andrew":86.87934 ]


Note that each of the values stored in the map can be of a different type. Brett's is an integer, Pete's is a string,
and Andrew's is a floating point number.

We can access the values in a map in two main ways, as shown below::

    println scores["Pete"]
    println scores.Pete


To modify or add data to a map, the syntax is similar to adding values to list.
::

    scores["Pete"] = 3
    scores["Cedric"] = 120


Learn more about maps at the following links:

* `Groovy Maps tutorial <http://http://groovy-lang.org/groovy-dev-kit.html#Collections-Maps>`_
* `Groovy Map SDK <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html>`_
* `Java Map SDK <http://docs.oracle.com/javase/7/docs/api/java/util/Map.html>`_


.. _script-multiple-assignment:

Multiple assignment
----------------------

An array or a list object can used to assign multiple variables at once. For example::

    (a, b, c) = [10, 20, 'foo']
    assert a == 10 && b == 20 && c == 'foo'

The three variables on the left of the assignment operator are initialized by the correspondent item in the list.

Read more about `Multiple assignment <http://www.groovy-lang.org/semantics.html#_multiple_assignment>`_ on the Groovy documentation.


Conditional Execution
----------------------

One of the most important features of any programming language is the ability to execute different code under
different conditions. The simplest way to do this is to use the ``if`` construct. For example::

    x = Math.random()
    if( x < 0.5 ) {
        println "You lost."
    }
    else {
        println "You won!"
    }



Strings
--------

Strings can be defined by delimiting them with a single-quote ``'`` or a double-quote ``"``.
Using either type of string allows you to use strings with quotations easily, as shown below::

    println "he said 'cheese' once"
    println 'he said "cheese!" again'


Strings may be concatenated with ``+``. For example::

    a = "world"
    print "hello " + a + "\n"


.. _string-interpolation:

String interpolation
----------------------

There is an important difference between single ``'`` and double ``"`` quoted strings.
Double quoted strings supports variable interpolations, while single quoted strings do not.

In practice, strings that are declared inside double-quotes can contain arbitrary variables prefixing them with the ``$`` character
or any expressions by using the ``${expression}`` syntax in a very similar way to Bash/shell scripts.
::

    foxtype = 'quick'
    foxcolor = ['b', 'r', 'o', 'w', 'n']
    println "The $foxtype ${foxcolor.join()} fox"

    x = 'Hello'
    println '$x + $y'

It prints::

    The quick brown fox
    $x + $y


Multi-line strings
-------------------

A block of text that span multiple lines can be defined by delimiting it with triple single or double quotes, as shown below::

    text = """
        hello there James
        how are you today?
        """

.. note:: Like before, multi-line strings delimited by double-quotes characters supports variable interpolation, while
   single-quote string do not.


As in Bash/shell scripts, when terminating a multi-line text block with a ``\`` character, the resulting string is
not broken up by `new line` character(s)::

    myLongCmdline = """ blastp \
                    -in $input_query \
                    -out $output_file \
                    -db $blast_database \
                    -html
                    """

    result = myLongCmdline.execute().text




.. _script-closure:

Closures
=========

In very few words a closure is a block of code that can be passed as an argument to a function.
Thus you can define a chunk of code and then pass it around as if it were a string or an integer.

More formally, you can create functions that are defined `first class objects`.

::

    square = { it * it }


The curly brackets around the expression ``it * it`` tells the script interpreter to treat this expression as code.
In this case, the designator ``it`` refers to whatever value is given to the function. Then this compiled function is
assigned to the variable `'square`` much like those above. So now we can do something like this::

    println square(9)

and get the value 81.


This is not very interesting until we find that we can pass this function ``square`` around as a method argument.
There are some built in functions that take a function like this as an argument. One example is the ``collect`` method on lists.
For example::

    [ 1, 2, 3, 4 ].collect(square)


This expression says, create an array with the values 1,2,3 and 4, then call the `collect` method, passing in the
closure we defined above. The collect method runs through each item in the array, calls the closure on the item,
then puts the result in a new array, resulting in::

    [ 1, 4, 9, 16 ]


For more methods that you can call with closures as arguments, see the `Groovy GDK documentation <http://groovy.codehaus.org/groovy-jdk/>`_.


By default closures take 1 parameter called ``it``, you can also create closures with named parameters.
For example the method ``Map.each()`` can take a closure with two variables, to which it binds the `key` and the associated `value`::


    printMapClosure = { key, value ->
        println "$key = $value"
    }

    [ "Yue" : "Wu", "Mark" : "Williams", "sudha" : "Kumari" ].each(printMapClosure)


Prints::


    Yue=Wu
    Mark=Williams
    Sudha=Kumari




A closure has another two important features. First it can access variables in the scope where it is defined and
so it can `interact` with them.

The second thing is that a closure can be defined in an `anonymous` manner, meaning that it is not given a name,
and is defined in the place where it needs to be used.

As an example showing both these features see the following code fragment::

    myMap = ["China": 1 , "India" : 2, "USA" : 3]

    result = 0
    myMap.keySet().each( { result+= myMap[it] } )

    println result



.. _script-regexp:

Regular expressions
====================

Regular expressions are the Swiss Army knife of text processing. They provide the programmer with the ability to match
and extract patterns from strings.

Regular expressions are supported by using the ``~/pattern/`` syntax and the ``=~`` and the ``==~`` operators.

Use the ``=~`` to check if there's any occurrence of a given pattern into a given string, thus::

    assert 'foo' =~ /foo/       // return TRUE
    assert 'foobar' =~ /foo/    // return TRUE


Use the ``==~`` to check whenever a string matches a given regular expression pattern.
::

    assert 'foo' ==~ /foo/       // return TRUE
    assert 'foobar' ==~ /foo/    // return FALSE


It is worth noting that the ``~`` operator creates a Java ``Pattern`` object from the given string,
while the ``=~`` creates a Java ``Matcher`` object.
::

    x = ~/abc/
    println x.class
    // prints java.util.regex.Pattern

    y = 'some string' =~ /abc/
    println y.class
    // prints java.util.regex.Matcher


Regular expression support is imported from Java. Java's regular expression language and API is documented in the
`Pattern Java documentation <http://download.oracle.com/javase/7/docs/api/java/util/regex/Pattern.html>`_.

You may also be interested in this post `Groovy: Don't Fear the RegExp <http://naleid.com/blog/2008/05/19/dont-fear-the-regexp>`_.


String replacements
--------------------

To replace pattern occurrences into a given string use the ``replaceFirst`` and ``replaceAll`` methods. For example::

     x = "colour".replaceFirst(/ou/, "o")
     println x
     // prints: color

     y = "cheesecheese".replaceAll(/cheese/,"nice")
     println y
     // prints: nicenice



Capturing groups
----------------

You can match a pattern that includes groups.  First create a matcher object with the ``=~`` operator.
Then, you can index the matcher object to find the matches, ``matcher[0]`` returns a list representing the first match
of the regular expression in the string. The first element is the string that matches the entire regular expression, and
the remaining elements are the strings that match each group.

Here's how it works::

    programVersion = '2.7.3-beta'
    m = programVersion =~ /(\d+)\.(\d+)\.(\d+)-?(.+)/

    assert m[0] ==  ['2.7.3-beta', '2', '7', '3', 'beta']
    assert m[0][1] == '2'
    assert m[0][2] == '7'
    assert m[0][3] == '3'
    assert m[0][4] == 'beta'



Applying some syntax sugar, you can do the same in just one line of code::

    programVersion = '2.7.3-beta'
    (full, major, minor, patch, flavor) = (programVersion =~ /(\d+)\.(\d+)\.(\d+)-?(.+)/)[0]

    println full    // 2.7.3-beta
    println major   // 2
    println minor   // 7
    println patch   // 3
    println flavor  // beta


Removing part of a string
-------------------------

You can subtract a part of a String value using a regular expression pattern. The first match found is
replaced with an empty String. For example::

    // define the regexp pattern
    wordStartsWithGr = ~/(?i)\s+Gr\w+/

    // apply and verify the result
    ('Hello Groovy world!' - wordStartsWithGr) == 'Hello world!'
    ('Hi Grails users' - wordStartsWithGr) == 'Hi users'



Remove first match of a word with 5 characters::

    assert ('Remove first match of 5 letter word' - ~/\b\w{5}\b/) == 'Remove  match of 5 letter word'


Remove first found numbers followed by a whitespace character::

    assert ('Line contains 20 characters' - ~/\d+\s+/) == 'Line contains characters'



.. _script-file-io:

Files and I/O
==============

In order to access and work with files, you need to use the ``file`` method which returns a file system object
given a file path string. For example::

  myFile = file('some/path/to/my_file.file')


The ``file`` method can reference either `files` or `directories` depending on what the string path is locating in the
file system.

When using wildcard characters i.e. ``*``, ``?``, ``[]`` and ``{}`` the argument is interpreted as a `glob`_ path matcher
and the ``file`` method returns a list object holding the path of files
whose name matches the specified pattern, or an empty list if no match is found. For example::

  listOfFiles = file('some/path/*.fa')

.. note:: Two asterisks, i.e. ``**``, works like ``*`` but crosses directory boundaries. Also note that,
  by default wildcard characters do not match against hidden and directory files.

For example, if you want to include hidden files in the result list, add the optional parameter ``hidden`` as shown below::

  listWithHidden = file('some/path/*.fa', hidden: true)

The list of available options is shown below:

=============== ===================
Name            Description
=============== ===================
glob            When ``true`` interprets characters ``*``, ``?``, ``[]`` and ``{}`` as glob wildcards, otherwise handles them as normal characters (default: ``true``)
type            Type of paths returned, either ``file``, ``dir`` or ``any`` (default: ``file``)
hidden          When ``true`` includes hidden files in the resulting paths (default: ``false``)
maxDepth        Maximum number of directory levels to visit (default: `no limit`)
followLinks     When ``true`` it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: ``true``)
=============== ===================


.. tip:: If you are a Java geek you will be interested to know that the ``file`` method returns a
  `Path <http://docs.oracle.com/javase/7/docs/api/java/nio/file/Path.html>`_ object, which allows
  you to use the usual methods as you would in a Java program.

See also: :ref:`Channel.fromPath <channel-path>` .

Basic read/write
------------------

Given a file variable, declared by using the ``file`` method as shown in the previous example, reading a file
is as easy as getting the value of the file's ``text`` property, which returns the file content
as a string value. For example::

  print myFile.text


In the same way you can save a string value to a file by simply assigning it to the file's ``text`` property,
as shown below::

  myFile.text = 'Hello world!'


.. note:: The actual file content is overwritten by the assignment operation. If the file does not exist, it is created
  implicitly by the assignment operation.

In order to append a string value to a file without erasing the actual content, you can use the ``append`` method::

  myFile.append('Add this line\n')

or by using the `left shift` operator, which is just a more idiomatic way to append text content to a file::

  myFile << 'Add a line more\n'


Binary data can managed in the same way, just using the file property ``bytes`` instead of ``text``. Thus, the following
example read the file and returns its content as a byte array::

  binaryContent = myFile.bytes

Or you can save a byte array data buffer to a file, by simply writing::

  myFile.bytes = binaryBuffer


.. warning:: The above methods read and write ALL the file content at once, in a single variable or buffer. For this
  reason they are not suggested when dealing with big files, which require a more memory efficient approach, for example
  reading a file line by line or by using a fixed size buffer.


Read a file line by line
--------------------------

In order to read a text file line by line you can use the method ``readLines()`` provided by the file object which
returns the file content as a list of strings. For example::

    myFile = file('some/my_file.txt')
    allLines  = myFile.readLines()
    for( line : allLines ) {
        println line
    }


The same example can be written in a more idiomatic syntax, as shown below::

    file('some/my_file.txt')
        .readLines()
        .each { println it }


.. note:: The method ``readLines()`` reads all the file content at once and returns a list containing all the lines. For
  this reason do not use it to read big files.


When you need to process a big file line by line, use the method ``eachLine`` which allows you to read a file
processing each line one by one, thus avoiding the loading of all the file content in the memory. For example::

    count = 0
    myFile.eachLine {  str ->
            println "line ${count++}: $str"
        }



Advanced file reading operations
-----------------------------------

The ``Reader`` and the ``InputStream`` classes allow you to gain fine control on read operations for
text and binary files respectively. 


The file method ``newReader`` creates a `Reader <http://docs.oracle.com/javase/7/docs/api/java/io/Reader.html>`_ object
for the given file and allows you to read the content in single characters, lines or arrays of characters. For example::

    myReader = myFile.newReader()
    String line
    while( line = myReader.readLine() ) {
        println line
    }
    myReader.close()


The method ``withReader`` works in a similar manner but saves you from calling the method ``close`` when you have finished
processing the file since it is managed automatically by the method itself. The same example can be rewritten as shown below::

    myFile.withReader {
        String line
        while( line = myReader.readLine() ) {
            println line
        }
    }

The methods ``newInputStream`` and ``withInputStream`` work in a similar manner. The main difference is that they create an
`InputStream <http://docs.oracle.com/javase/7/docs/api/java/io/InputStream.html>`_ object that is useful to write binary
data.

Table of the most important methods to read file content:

=============== ==============
Name            Description
=============== ==============
getText         Returns the file content as a string value
getBytes        Returns the file content as byte array
readLines       Reads the file line by line and returns the content as a list of strings
eachLine        Iterates over the file line by line applying the specified :ref:`closure <script-closure>`
eachByte        Iterates over the file by each single byte applying the specified :ref:`closure <script-closure>`
withReader      Opens a file for reading and lets you access it with a `Reader <http://docs.oracle.com/javase/7/docs/api/java/io/Reader.html>`_ object
withInputStream Opens a file for reading and lets you access it with an `InputStream <http://docs.oracle.com/javase/7/docs/api/java/io/InputStream.html>`_ object
newReader       Returns a `Reader <http://docs.oracle.com/javase/7/docs/api/java/io/Reader.html>`_ object to read a text file
newInputStream  Returns an `InputStream <http://docs.oracle.com/javase/7/docs/api/java/io/InputStream.html>`_ object to read a binary file
=============== ==============


Read the Java documentation for `Reader <http://docs.oracle.com/javase/7/docs/api/java/io/Reader.html>`_ and
`InputStream <http://docs.oracle.com/javase/7/docs/api/java/io/InputStream.html>`_ classes to learn more about
the methods concerning these classes.


Advanced file writing operations
----------------------------------

When you need to access low-level write operations to handle single bytes or characters, or if you are working with big files you will
need to use the ``Writer`` and ``OutputStream`` classes which provide fine control on write operations.

For example, given two file objects ``sourceFile`` and ``targetFile``, the following code snippet shows how to copy the
file content from the first file into the second one replacing all the ``U`` characters with ``X``::

    sourceFile.withReader { source ->
        targetFile.withWriter { target ->
            String line
            while( line=source.readLine() ) {
                target << line.replaceAll('U','X')
            }
        }
    }


Table of the most important methods to write data into a file:

=================== ==============
Name                Description
=================== ==============
setText             Saves a string value to a file
setBytes            Saves a bytes array to a file
write               Saves a string to a file truncating the actual content
append              Appends a string value to a file without truncating the actual content
newWriter           Creates a `Writer <http://docs.oracle.com/javase/7/docs/api/java/io/Writer.html>`_ object that allows you to save text data to a file
newPrintWriter      Creates a `PrintWriter <http://docs.oracle.com/javase/7/docs/api/java/io/PrintWriter.html>`_ object that allows you to write formatted text to a file
newOutputStream     Creates an `OutputStream <http://docs.oracle.com/javase/7/docs/api/java/io/OutputStream.html>`_ object that allows you to write binary data to a file
withWriter          Applies the specified closure to a `Writer <http://docs.oracle.com/javase/7/docs/api/java/io/Writer.html>`_ object, closing it when finished.
withPrintWriter     Applies the specified closure to a `PrintWriter <http://docs.oracle.com/javase/7/docs/api/java/io/PrintWriter.html>`_ object, closing it when finished.
withOutputStream    Applies the specified closure to an `OutputStream <http://docs.oracle.com/javase/7/docs/api/java/io/OutputStream.html>`_ object, closing it when finished.
=================== ==============

Read the Java documentation for the `Writer <http://docs.oracle.com/javase/7/docs/api/java/io/Writer.html>`_,
`PrintWriter <http://docs.oracle.com/javase/7/docs/api/java/io/PrintWriter.html>`_ and
`OutputStream <http://docs.oracle.com/javase/7/docs/api/java/io/OutputStream.html>`_ classes to learn more about
methods available for writing data.


List directory content
----------------------

Let's assume that you need to walk through a directory of your choice. You can define the ``myDir`` variable
that points to it::

    myDir = file('any/path')

The simplest way to get a directory list is by using the methods ``list`` or ``listFiles``,
that return a collection of first-level elements (files and directories) in a directory. For example::

    allFiles = myDir.list()
    for( def file : allFiles ) {
        println file
    }

.. note:: The only difference between ``list`` and ``listFiles`` is that the first returns a list of strings, while the latter a
   list of file objects, that allow you to access file dependent data e.g. size, last modified time, etc.


The ``eachFile`` method allows you to iterate through the first-level elements only
(just like ``listFiles``). As with other `each-` methods, they take a closure as an input parameter. For example::

    myDir.eachFile { item ->
        if( item.isFile() ) {
            println "${item.getName()} - size: ${item.size()}"
        }
        else if( item.isDirectory() ) {
            println "${item.getName()} - DIR"
        }
    }


Several variants of the above method are available. See the table below for a complete list.

=================== ==================
Name                Description
=================== ==================
eachFile            Iterates through first-level elements (files and directories). `Read more <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachFile(groovy.io.FileType,%20groovy.lang.Closure)>`_
eachDir             Iterates through first-level directories only. `Read more <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachDir(groovy.lang.Closure)>`_
eachFileMatch       Iterates through files and dirs whose name match the given filter. `Read more <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachFileMatch(java.lang.Object,%20groovy.lang.Closure)>`_
eachDirMatch        Iterates through directories whose name match the given filter. `Read more <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachDirMatch(java.lang.Object,%20groovy.lang.Closure)>`_
eachFileRecurse     Iterates through directory elements in a depth-first fashion. `Read more <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachFileRecurse(groovy.lang.Closure)>`_
eachDirRecurse      Iterates through directories in a depth-first fashion (regular files are ignored). `Read more <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachDirRecurse(groovy.lang.Closure)>`_
=================== ==================


See also: Channel :ref:`channel-path` method.


Create directories
------------------

Given a file variable representing a nonexistent directory, like the following::

  myDir = file('any/path')

the method ``mkdir`` allows you to create a folder at the given path. It returns a ``true`` value if the folder is created
successfully, of ``false`` otherwise. ::

   result = myDir.mkdir()
   println result ? "OK" : "Cannot create folder: $myDir"

.. note:: If the parent directories do not exist, the above method will fail returning a ``false`` value.

The method ``mkdirs`` allows you to create the directory named by the file object, including any necessary but
nonexistent parent directories. For example::

    myDir.mkdirs()


Create links
------------

Given a file, the method ``mklink`` creates a *file system link* for that file using the path specified as a parameter.
For example::

  myFile = file('/some/path/file.txt')
  myFile.mklink('/user/name/link-to-file.txt')


Table of optional parameters:

==================  ================
Method              Description
==================  ================
hard                When ``true`` it creates an *hard* link, otherwise it creates a *soft* (aka *symbolic*) link. (default: ``false``).
overwrite           When ``true`` it overwrite any exiting file with the same name, otherwise it will throws a `FileAlreadyExistsException <http://docs.oracle.com/javase/8/docs/api/java/nio/file/FileAlreadyExistsException.html>`_ (default: ``false``).
==================  ================


Copy files
----------

The method ``copyTo`` allows you to copy a file into a new file or into a directory, or copy a directory to a new
directory. Having a file variable ``myFile``, the following example shows how to copy a file into a new file
with a different name::

  myFile.copyTo('new_name.txt')


.. note:: If the target file already exists, it will be replaced by the new one. Note also that if the target is
  a directory, the source file will be copied into that folder maintaining the original name.


When the source file is a directory, all the directory content is copied to the target folder. For example::


  myDir = file('/some/path')
  myDir.copyTo('/some/new/path')


If the target path does not exist, it will be created automatically.

.. tip:: The ``copyTo`` method mimics the semantic of the Linux command ``cp -r <source> <target>``,
  with the following caveats: Unix BASH distinguish between paths having or not having and ending slash, for example:
  ``/some/path/name`` and ``/some/path/name/``. The first locates a regular file while the latter identifies a directory
  location. With Nextflow, due to Java files API implementation, this is not possibles and both strings represents the same path.
  If that path exists on the file systems it is handled accordingly (as a regular file or as a directory). If the path does not
  exist, it is supposed to locate a regular file (and any parent directory will be created automatically).



Move files
----------

You can move a file by using the method ``moveTo`` as shown in the example below::

  myFile = file('/some/path/file.txt')
  myFile.moveTo('/another/path/new_file.txt')


.. note:: When a file with the same name as the target already exists, it will be replaced by the new one. Note also that
   when the target specifies a folder name instead of a file, the source file is moved in that folder maintaining the original name.

When the source file is a directory, all the directory content is moved to the new the destination folder::

  myDir = file('/any/dir_a')
  myDir.moveTo('/any/dir_b')


Please note that the result of the above example depends on the existence of the destination folder. If the destination
folder exists, the source is moved into the destination folder, thus the resulting path will be::

  /any/dir_b/dir_a

If the destination folder does not exist, the source is just renamed to the target name and so the resulting
path of the above move operation would be::

    /any/dir_b


.. tip:: The ``moveTo`` uses the same semantic as the Linux command ``mv <source> <target>``. The same caveats as for
  the method ``copyTo`` are applied.


Rename files
------------

You can rename a file or a directory by simply using the ``renameTo`` file method a shown below::

  myFile = file('my_file.txt')
  myFile.renameTo('new_file_name.txt')


Delete files
------------

The file method ``delete`` allows you to delete a file or a directory with a given path. It returns the value ``true``
when the operation is completed  successfully or ``false`` if it fails to delete it. For example::

  myFile = file('some/file.txt')
  result = myFile.delete
  println result ? "OK" : "Can delete: $myFile"


.. note:: This method deletes a directory ONLY if it does not contain any file or sub-directory. In order to delete a
  directory and ALL its content, i.e. removing all the files and sub-directories it may contain, use the method ``deleteDir``
  instead.


Check file attributes
---------------------

The following methods can be used on a file variable created by using the ``file`` method.

==================  ================
Method              Description
==================  ================
getName             Gets the file name e.g. ``/some/path/file.txt`` -> ``file.txt``
getBaseName         Gets the file name e.g. ``/some/path/file.txt`` -> ``file``
getExtension        Gets the file extension e.g. ``/some/path/file.txt`` -> ``txt``
getParent           Gets the file parent path e.g. ``/some/path/file.txt`` -> ``/some/path``
size                Gets the file size in bytes
exists              Returns ``true`` if the file exists, or ``false`` otherwise
isEmpty             Returns ``true`` if the file is zero length or does not exist, ``false`` otherwise
isFile              Returns ``true`` if it is a regular file e.g. not a directory
isDirectory         Returns ``true`` if the file is a directory
isHidden            Returns ``true`` if the file is hidden
lastModified        Returns the file last modified timestamp i.e. a long as Linux epoch time
==================  ================


For example, the following line prints a file name and size::

  println "File ${myFile.getName() size: ${myFile.size()}"



Get and modify file permissions
-------------------------------

Given a file variable representing any file or a directory, the method
``getPermissions`` returns a string of nine characters representing the file permission using the
`Linux symbolic notation <http://en.wikipedia.org/wiki/File_system_permissions#Symbolic_notation>`_
e.g. ``rw-rw-r--``. For example::

    permissions = myFile.getPermissions()


In the same way the method ``setPermissions`` allows you to set the file access permissions using the same string
notation. For example::

    myFile.setPermissions('rwxr-xr-x')


A second version of the ``setPermissions`` method allows you to set file permissions specifying three digits, representing
respectively the `owner`, `group` and `other` permissions. For example::

    myFile.setPermissions(7,5,5)


Learn more about `File permissions numeric notation <http://en.wikipedia.org/wiki/File_system_permissions#Numeric_notation>`_.

HTTP/FTP files
--------------

Nextflow provides a transparent integration for HTTP/S and FTP protocols that allows you to handle remote resources
as local file system objects. Simply specify the resource URL as the argument of the `file` object. For example::

    pdb = file('http://files.rcsb.org/header/5FID.pdb')

Then, you will be able to access it as a local file as described in the previous sections. For example::

    println pdb.text

The above one-liner prints the content of the PDB file. See in the previous sections how to stream or copy the content
of files.

.. note:: Write and list operations are not supported by HTTP/S and FTP files.

.. _glob: http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
