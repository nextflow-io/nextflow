.. _pipeline-page:

*****************
Pipeline script
*****************


The Nextflow scripting language is an extension of the Groovy programming language whose syntax has been
specialized to ease the writing of computational pipelines in a declarative manner.

This means that Nextflow can execute any piece of Groovy code or use any library for the JVM platform.

For a detailed description of the Groovy programming language, reference these links:

* `Groovy User Guide <http://groovy-lang.org/documentation.html>`_
* `Groovy Cheat sheet <http://refcardz.dzone.com/refcardz/groovy>`_
* `Groovy in Action <http://www.manning.com/koenig2/>`_


Below you can find a crash course in the most important language constructs used in the Nextflow scripting language.

.. warning:: Nextflow uses ``UTF-8`` as the default file character encoding for source and application files. Make sure
  to use the ``UTF-8`` encoding when editing Nextflow scripts with your favourite text editor.

Language basics
==================


Hello world
------------

To print something is as easy as using one of the ``print`` or ``println`` methods.
::

    println "Hello, World!"

The only difference between the two is that the ``println`` method implicitly appends a `new line` character
to the printed string.


Variables
----------

To define a variable, simply assign a value to it::

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

A List object can be defined by placing the list items in square brackets::

    myList = [1776, -1, 33, 99, 0, 928734928763]

You can access a given item in the list with square-bracket notation (indexes start at 0)::

    println myList[0]

In order to get the length of the list use the ``size`` method::

    println myList.size()


Learn more about lists:

* `Groovy Lists tutorial <http://groovy-lang.org/groovy-dev-kit.html#Collections-Lists>`_
* `Groovy List SDK <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/List.html>`_
* `Java List SDK <http://docs.oracle.com/javase/7/docs/api/java/util/List.html>`_


Maps
-----

Maps are used to store `associative arrays` or `dictionaries`. They are unordered collections of heterogeneous, named data::

    scores = [ "Brett":100, "Pete":"Did not finish", "Andrew":86.87934 ]


Note that each of the values stored in the map can be of a different type. ``Brett`` is an integer, ``Pete`` is a string,
and ``Andrew`` is a floating-point number.

We can access the values in a map in two main ways::

    println scores["Pete"]
    println scores.Pete


To add data to or modify a map, the syntax is similar to adding values to list::

    scores["Pete"] = 3
    scores["Cedric"] = 120


Learn more about maps:

* `Groovy Maps tutorial <http://groovy-lang.org/groovy-dev-kit.html#Collections-Maps>`_
* `Groovy Map SDK <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html>`_
* `Java Map SDK <http://docs.oracle.com/javase/7/docs/api/java/util/Map.html>`_


.. _script-multiple-assignment:

Multiple assignment
----------------------

An array or a list object can used to assign to multiple variables at once::

    (a, b, c) = [10, 20, 'foo']
    assert a == 10 && b == 20 && c == 'foo'

The three variables on the left of the assignment operator are initialized by the corresponding item in the list.

Read more about `Multiple assignment <http://www.groovy-lang.org/semantics.html#_multiple_assignment>`_ in the Groovy documentation.


Conditional Execution
----------------------

One of the most important features of any programming language is the ability to execute different code under
different conditions. The simplest way to do this is to use the ``if`` construct::

    x = Math.random()
    if( x < 0.5 ) {
        println "You lost."
    }
    else {
        println "You won!"
    }



Strings
-------

Strings can be defined by enclosing text in single or double quotes (``'`` or ``"`` characters)::

    println "he said 'cheese' once"
    println 'he said "cheese!" again'


Strings can be concatenated with ``+``::

    a = "world"
    print "hello " + a + "\n"


.. _string-interpolation:

String interpolation
--------------------

There is an important difference between single- and double-quoted strings:
Double-quoted strings support variable interpolations, while single-quoted strings do not.

In practice, double-quoted strings can contain the value of an arbitrary variable by prefixing its name with the ``$`` character,
or the value of any expression by using the ``${expression}`` syntax, similar to Bash/shell scripts::

    foxtype = 'quick'
    foxcolor = ['b', 'r', 'o', 'w', 'n']
    println "The $foxtype ${foxcolor.join()} fox"

    x = 'Hello'
    println '$x + $y'

This code prints::

    The quick brown fox
    $x + $y


Multi-line strings
-------------------

A block of text that span multiple lines can be defined by delimiting it with triple single or double quotes::

    text = """
        hello there James
        how are you today?
        """

.. note:: Like before, multi-line strings inside double quotes support variable interpolation, while
   single-quoted multi-line strings do not.


As in Bash/shell scripts, terminating a line in a multi-line string with a ``\`` character prevents a
a `new line` character from separating that line from the one that follows::

    myLongCmdline = """ blastp \
                    -in $input_query \
                    -out $output_file \
                    -db $blast_database \
                    -html
                    """

    result = myLongCmdline.execute().text

In the preceding example, ``blastp`` and its ``-in``, ``-out``, ``-db`` and ``-html`` switches and
their arguments are effectively a single line.


.. _script-closure:

Closures
=========

Briefly, a closure is a block of code that can be passed as an argument to a function.
Thus, you can define a chunk of code and then pass it around as if it were a string or an integer.

More formally, you can create functions that are defined as `first class objects`.

::

    square = { it * it }


The curly brackets around the expression ``it * it`` tells the script interpreter to treat this expression as code.
In this case, the designator ``it`` refers to whatever value is passed to the function when it is called. This compiled function is
assigned to the variable ``square``, much like variable assignments shown previously. Now we can do something like this::

    println square(9)

and get the value 81.


This is not very interesting until we find that we can pass the function ``square`` as an argument to other functions or methods.
Some built-in functions take a function like this as an argument. One example is the ``collect`` method on lists::

    [ 1, 2, 3, 4 ].collect(square)


This expression says: Create an array with the values 1, 2, 3 and 4, then call its ``collect`` method, passing in the
closure we defined above. The ``collect`` method runs through each item in the array, calls the closure on the item,
then puts the result in a new array, resulting in::

    [ 1, 4, 9, 16 ]


For more methods that you can call with closures as arguments, see the `Groovy GDK documentation <http://docs.groovy-lang.org/latest/html/groovy-jdk/>`_.


By default, closures take a single parameter called ``it``, but you can also create closures with multiple, custom-named parameters.
For example, the method ``Map.each()`` can take a closure with two arguments, to which it binds the `key` and the associated `value`
for each key-value pair in the ``Map``. Here, we use the obvious variable names ``key`` and ``value`` in our closure::


    printMapClosure = { key, value ->
        println "$key = $value"
    }

    [ "Yue" : "Wu", "Mark" : "Williams", "sudha" : "Kumari" ].each(printMapClosure)


Prints::


    Yue=Wu
    Mark=Williams
    Sudha=Kumari

A closure has two other important features. First, it can access variables in the scope where it is defined,
so that it can `interact` with them.

Second, a closure can be defined in an `anonymous` manner, meaning that it is not given a name,
and is defined in the place where it needs to be used.

As an example showing both these features, see the following code fragment::

    myMap = ["China": 1 , "India" : 2, "USA" : 3]

    result = 0
    myMap.keySet().each( { result+= myMap[it] } )

    println result


Learn more about closures in the `Groovy documentation <http://groovy-lang.org/closures.html>`_

.. _script-regexp:

Regular expressions
====================

Regular expressions are the Swiss Army knife of text processing. They provide the programmer with the ability to match
and extract patterns from strings.

Regular expressions are available via the ``~/pattern/`` syntax and the ``=~`` and ``==~`` operators.

Use ``=~`` to check whether a given pattern occurs anywhere in a string::

    assert 'foo' =~ /foo/       // return TRUE
    assert 'foobar' =~ /foo/    // return TRUE


Use ``==~`` to check whether a string matches a given regular expression pattern exactly.
::

    assert 'foo' ==~ /foo/       // return TRUE
    assert 'foobar' ==~ /foo/    // return FALSE


It is worth noting that the ``~`` operator creates a Java ``Pattern`` object from the given string,
while the ``=~`` operator creates a Java ``Matcher`` object.
::

    x = ~/abc/
    println x.class
    // prints java.util.regex.Pattern

    y = 'some string' =~ /abc/
    println y.class
    // prints java.util.regex.Matcher


Regular expression support is imported from Java. Java's regular expression language and API is documented in the
`Pattern Java documentation <http://download.oracle.com/javase/7/docs/api/java/util/regex/Pattern.html>`_.

You may also be interested in this post: `Groovy: Don't Fear the RegExp <https://web.archive.org/web/20170621185113/http://www.naleid.com/blog/2008/05/19/dont-fear-the-regexp>`_.


String replacement
--------------------

To replace pattern occurrences in a given string, use the ``replaceFirst`` and ``replaceAll`` methods::

     x = "colour".replaceFirst(/ou/, "o")
     println x
     // prints: color

     y = "cheesecheese".replaceAll(/cheese/, "nice")
     println y
     // prints: nicenice



Capturing groups
----------------

You can match a pattern that includes groups.  First create a matcher object with the ``=~`` operator.
Then, you can index the matcher object to find the matches: ``matcher[0]`` returns a list representing the first match
of the regular expression in the string. The first list element is the string that matches the entire regular expression, and
the remaining elements are the strings that match each group.

Here's how it works::

    programVersion = '2.7.3-beta'
    m = programVersion =~ /(\d+)\.(\d+)\.(\d+)-?(.+)/

    assert m[0] ==  ['2.7.3-beta', '2', '7', '3', 'beta']
    assert m[0][1] == '2'
    assert m[0][2] == '7'
    assert m[0][3] == '3'
    assert m[0][4] == 'beta'



Applying some syntactic sugar, you can do the same in just one line of code::

    programVersion = '2.7.3-beta'
    (full, major, minor, patch, flavor) = (programVersion =~ /(\d+)\.(\d+)\.(\d+)-?(.+)/)[0]

    println full    // 2.7.3-beta
    println major   // 2
    println minor   // 7
    println patch   // 3
    println flavor  // beta


Removing part of a string
-------------------------

You can remove part of a ``String`` value using a regular expression pattern. The first match found is
replaced with an empty String::

    // define the regexp pattern
    wordStartsWithGr = ~/(?i)\s+Gr\w+/

    // apply and verify the result
    ('Hello Groovy world!' - wordStartsWithGr) == 'Hello world!'
    ('Hi Grails users' - wordStartsWithGr) == 'Hi users'


Remove the first 5-character word from a string::

    assert ('Remove first match of 5 letter word' - ~/\b\w{5}\b/) == 'Remove  match of 5 letter word'


Remove the first number with its trailing whitespace from a string::

    assert ('Line contains 20 characters' - ~/\d+\s+/) == 'Line contains characters'



.. _script-file-io:

Files and I/O
==============

To access and work with files, use the ``file`` method, which returns a file system object
given a file path string::

  myFile = file('some/path/to/my_file.file')


The ``file`` method can reference either `files` or `directories`, depending on what the string path refers to in the
file system.

When using the wildcard characters ``*``, ``?``, ``[]`` and ``{}``, the argument is interpreted as a `glob`_ path matcher
and the ``file`` method returns a list object holding the paths of files whose names match the specified pattern, or an
empty list if no match is found::

  listOfFiles = file('some/path/*.fa')

.. note:: Two asterisks (``**``) in a glob pattern works like ``*`` but matches any number of directory components in a
          file system path.

By default, wildcard characters do not match directories or hidden files. For example, if you want to include hidden
files in the result list, add the optional parameter ``hidden``::

  listWithHidden = file('some/path/*.fa', hidden: true)

Here are ``file``'s available options:

=============== ===================
Name            Description
=============== ===================
glob            When ``true`` interprets characters ``*``, ``?``, ``[]`` and ``{}`` as glob wildcards, otherwise handles them as normal characters (default: ``true``)
type            Type of paths returned, either ``file``, ``dir`` or ``any`` (default: ``file``)
hidden          When ``true`` includes hidden files in the resulting paths (default: ``false``)
maxDepth        Maximum number of directory levels to visit (default: `no limit`)
followLinks     When ``true`` follows symbolic links during directory tree traversal, otherwise treats them as files (default: ``true``)
=============== ===================


.. tip:: If you are a Java geek you will be interested to know that the ``file`` method returns a
  `Path <http://docs.oracle.com/javase/7/docs/api/java/nio/file/Path.html>`_ object, which allows
  you to use the usual methods you would in a Java program.

See also: :ref:`Channel.fromPath <channel-path>` .

Basic read/write
------------------

Given a file variable, declared using the ``file`` method as shown in the previous example, reading a file
is as easy as getting the value of the file's ``text`` property, which returns the file content
as a string value::

  print myFile.text


Similarly, you can save a string value to a file by simply assigning it to the file's ``text`` property::

  myFile.text = 'Hello world!'


.. note:: Existing file content is overwritten by the assignment operation, which also implicitly creates
          files that do not exist.

In order to append a string value to a file without erasing existing content, you can use the ``append`` method::

  myFile.append('Add this line\n')

Or use the `left shift` operator, a more idiomatic way to append text content to a file::

  myFile << 'Add a line more\n'


Binary data can managed in the same way, just using the file property ``bytes`` instead of ``text``. Thus, the following
example reads the file and returns its content as a byte array::

  binaryContent = myFile.bytes

Or you can save a byte array data buffer to a file, by simply writing::

  myFile.bytes = binaryBuffer


.. warning:: The above methods read and write ALL the file content at once, in a single variable or buffer. For this
  reason they are not suggested when dealing with big files, which require a more memory efficient approach, for example
  reading a file line by line or by using a fixed size buffer.


Read a file line by line
--------------------------

In order to read a text file line by line you can use the method ``readLines()`` provided by the file object, which
returns the file content as a list of strings::

    myFile = file('some/my_file.txt')
    allLines  = myFile.readLines()
    for( line : allLines ) {
        println line
    }


This can also be written in a more idiomatic syntax::

    file('some/my_file.txt')
        .readLines()
        .each { println it }


.. note:: The method ``readLines()`` reads all the file content at once and returns a list containing all the lines. For
  this reason, do not use it to read big files.


To process a big file, use the method ``eachLine``, which reads only a single line at a time into memory::

    count = 0
    myFile.eachLine {  str ->
            println "line ${count++}: $str"
        }



Advanced file reading operations
-----------------------------------

The classes ``Reader`` and ``InputStream`` provide fine control for reading text and binary files, respectively._


The method ``newReader`` creates a `Reader <http://docs.oracle.com/javase/7/docs/api/java/io/Reader.html>`_ object
for the given file that allows you to read the content as single characters, lines or arrays of characters::

    myReader = myFile.newReader()
    String line
    while( line = myReader.readLine() ) {
        println line
    }
    myReader.close()


The method ``withReader`` works similarly, but automatically calls the ``close`` method for you when you have finished
processing the file. So, the previous example can be written more simply as::

    myFile.withReader {
        String line
        while( line = myReader.readLine() ) {
            println line
        }
    }

The methods ``newInputStream`` and ``withInputStream`` work similarly. The main difference is that they create an
`InputStream <http://docs.oracle.com/javase/7/docs/api/java/io/InputStream.html>`_ object useful for writing binary
data.

Here are the most important methods for reading from files:

=============== ==============
Name            Description
=============== ==============
getText         Returns the file content as a string value
getBytes        Returns the file content as byte array
readLines       Reads the file line by line and returns the content as a list of strings
eachLine        Iterates over the file line by line, applying the specified :ref:`closure <script-closure>`
eachByte        Iterates over the file byte by byte, applying the specified :ref:`closure <script-closure>`
withReader      Opens a file for reading and lets you access it with a `Reader <http://docs.oracle.com/javase/7/docs/api/java/io/Reader.html>`_ object
withInputStream Opens a file for reading and lets you access it with an `InputStream <http://docs.oracle.com/javase/7/docs/api/java/io/InputStream.html>`_ object
newReader       Returns a `Reader <http://docs.oracle.com/javase/7/docs/api/java/io/Reader.html>`_ object to read a text file
newInputStream  Returns an `InputStream <http://docs.oracle.com/javase/7/docs/api/java/io/InputStream.html>`_ object to read a binary file
=============== ==============


Read the Java documentation for `Reader <http://docs.oracle.com/javase/7/docs/api/java/io/Reader.html>`_ and
`InputStream <http://docs.oracle.com/javase/7/docs/api/java/io/InputStream.html>`_ classes to learn more about
methods available for reading data from files.


Advanced file writing operations
----------------------------------

The ``Writer`` and ``OutputStream`` classes provide fine control for writing text and binary files,
respectively, including low-level operations for single characters or bytes, and support for big files.

For example, given two file objects ``sourceFile`` and ``targetFile``, the following code copies the
first file's content into the second file, replacing all ``U`` characters with ``X``::

    sourceFile.withReader { source ->
        targetFile.withWriter { target ->
            String line
            while( line=source.readLine() ) {
                target << line.replaceAll('U','X')
            }
        }
    }


Here are the most important methods for writing to files:

=================== ==============
Name                Description
=================== ==============
setText             Writes a string value to a file
setBytes            Writes a byte array to a file
write               Writes a string to a file, replacing any existing content
append              Appends a string value to a file without replacing existing content
newWriter           Creates a `Writer <http://docs.oracle.com/javase/7/docs/api/java/io/Writer.html>`_ object that allows you to save text data to a file
newPrintWriter      Creates a `PrintWriter <http://docs.oracle.com/javase/7/docs/api/java/io/PrintWriter.html>`_ object that allows you to write formatted text to a file
newOutputStream     Creates an `OutputStream <http://docs.oracle.com/javase/7/docs/api/java/io/OutputStream.html>`_ object that allows you to write binary data to a file
withWriter          Applies the specified closure to a `Writer <http://docs.oracle.com/javase/7/docs/api/java/io/Writer.html>`_ object, closing it when finished
withPrintWriter     Applies the specified closure to a `PrintWriter <http://docs.oracle.com/javase/7/docs/api/java/io/PrintWriter.html>`_ object, closing it when finished
withOutputStream    Applies the specified closure to an `OutputStream <http://docs.oracle.com/javase/7/docs/api/java/io/OutputStream.html>`_ object, closing it when finished
=================== ==============

Read the Java documentation for the `Writer <http://docs.oracle.com/javase/7/docs/api/java/io/Writer.html>`_,
`PrintWriter <http://docs.oracle.com/javase/7/docs/api/java/io/PrintWriter.html>`_ and
`OutputStream <http://docs.oracle.com/javase/7/docs/api/java/io/OutputStream.html>`_ classes to learn more about
methods available for writing data to files.


List directory content
----------------------

Let's assume that you need to walk through a directory of your choice. You can define the ``myDir`` variable
that points to it::

    myDir = file('any/path')

The simplest way to get a directory list is by using the methods ``list`` or ``listFiles``,
which return a collection of first-level elements (files and directories) of a directory::

    allFiles = myDir.list()
    for( def file : allFiles ) {
        println file
    }

.. note:: The only difference between ``list`` and ``listFiles`` is that the former returns a list of strings, and the latter a
   list of file objects that allow you to access file metadata, e.g. size, last modified time, etc.


The ``eachFile`` method allows you to iterate through the first-level elements only
(just like ``listFiles``). As with other `each-` methods, ``eachFiles`` takes a closure as a parameter::

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
eachFileMatch       Iterates through files and dirs whose names match the given filter. `Read more <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachFileMatch(java.lang.Object,%20groovy.lang.Closure)>`_
eachDirMatch        Iterates through directories whose names match the given filter. `Read more <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachDirMatch(java.lang.Object,%20groovy.lang.Closure)>`_
eachFileRecurse     Iterates through directory elements depth-first. `Read more <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachFileRecurse(groovy.lang.Closure)>`_
eachDirRecurse      Iterates through directories depth-first (regular files are ignored). `Read more <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachDirRecurse(groovy.lang.Closure)>`_
=================== ==================


See also: Channel :ref:`channel-path` method.


Create directories
------------------

Given a file variable representing a nonexistent directory, like the following::

  myDir = file('any/path')

the method ``mkdir`` creates a directory at the given path, returning ``true`` if the directory is created
successfully, and ``false`` otherwise::

   result = myDir.mkdir()
   println result ? "OK" : "Cannot create directory: $myDir"

.. note:: If the parent directories do not exist, the above method will fail and return ``false``.

The method ``mkdirs`` creates the directory named by the file object, including any nonexistent parent directories::

    myDir.mkdirs()


Create links
------------

Given a file, the method ``mklink`` creates a *file system link* for that file using the path specified as a parameter::

  myFile = file('/some/path/file.txt')
  myFile.mklink('/user/name/link-to-file.txt')


Table of optional parameters:

==================  ================
Name                Description
==================  ================
hard                When ``true`` creates a *hard* link, otherwise creates a *soft* (aka *symbolic*) link. (default: ``false``)
overwrite           When ``true`` overwrites any existing file with the same name, otherwise throws a `FileAlreadyExistsException <http://docs.oracle.com/javase/8/docs/api/java/nio/file/FileAlreadyExistsException.html>`_ (default: ``false``)
==================  ================


Copy files
----------

The method ``copyTo`` copies a file into a new file or into a directory, or copie a directory to a new
directory::

  myFile.copyTo('new_name.txt')


.. note:: If the target file already exists, it will be replaced by the new one. Note also that, if the target is
  a directory, the source file will be copied into that directory, maintaining the file's original name.


When the source file is a directory, all its content is copied to the target directory::

  myDir = file('/some/path')
  myDir.copyTo('/some/new/path')


  If the target path does not exist, it will be created automatically.

.. tip:: The ``copyTo`` method mimics the semantics of the Linux command ``cp -r <source> <target>``, with the
         following caveat: While Linux tools often treat paths ending with a slash (e.g. ``/some/path/name/``)
         as directories, and those not (e.g. ``/some/path/name``) as regular files, Nextflow (due to its use of
         the Java files API) views both these paths as the same file system object. If the path exists, it is
         handled according to its actual type (i.e. as a regular file or as a directory). If the path does not
         exist, it is treated as a regular file, with any missing parent directories created automatically.



Move files
----------

You can move a file by using the method ``moveTo``::

  myFile = file('/some/path/file.txt')
  myFile.moveTo('/another/path/new_file.txt')


.. note:: When a file with the same name as the target already exists, it will be replaced by the source. Note
          also that, when the target is a directory, the file will be moved to (or within) that directory,
          maintaining the file's original name.

When the source is a directory, all the directory content is moved to the target directory::

  myDir = file('/any/dir_a')
  myDir.moveTo('/any/dir_b')


Please note that the result of the above example depends on the existence of the target directory. If the target
directory exists, the source is moved into the target directory, resulting in the path::

  /any/dir_b/dir_a

If the target directory does not exist, the source is just renamed to the target name, resulting in the path::

  /any/dir_b


.. tip:: The ``moveTo`` method mimics the semantics of the Linux command ``mv <source> <target>``, with the
         same caveat as that given for ``copyTo``, above.


Rename files
------------

You can rename a file or directory by simply using the ``renameTo`` file method::

  myFile = file('my_file.txt')
  myFile.renameTo('new_file_name.txt')


Delete files
------------

The file method ``delete`` deletes the file or directory at the given path, returning ``true`` if the
operation succeeds, and ``false`` otherwise::

  myFile = file('some/file.txt')
  result = myFile.delete
  println result ? "OK" : "Can delete: $myFile"


.. note:: This method deletes a directory ONLY if it does not contain any files or sub-directories. To
          delete a directory and ALL its content (i.e. removing all the files and sub-directories it may
          contain), use the method ``deleteDir``.


Check file attributes
---------------------

The following methods can be used on a file variable created by using the ``file`` method:

==================  ================
Name                Description
==================  ================
getName             Gets the file name e.g. ``/some/path/file.txt`` -> ``file.txt``
getBaseName         Gets the file name without its extension e.g. ``/some/path/file.txt`` -> ``file``
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

Given a file variable representing a file (or directory), the method ``getPermissions`` returns a
9-character string representing the file's permissions using the
`Linux symbolic notation <http://en.wikipedia.org/wiki/File_system_permissions#Symbolic_notation>`_
e.g. ``rw-rw-r--``::

    permissions = myFile.getPermissions()


Similarly, the method ``setPermissions`` sets the file's permissions using the same notation::

    myFile.setPermissions('rwxr-xr-x')


A second version of the ``setPermissions`` method sets a file's permissions given three digits representing,
respectively, the `owner`, `group` and `other` permissions::

    myFile.setPermissions(7,5,5)


Learn more about `File permissions numeric notation <http://en.wikipedia.org/wiki/File_system_permissions#Numeric_notation>`_.

HTTP/FTP files
--------------

Nextflow provides transparent integration of HTTP/S and FTP protocols for handling remote resources
as local file system objects. Simply specify the resource URL as the argument of the `file` object::

    pdb = file('http://files.rcsb.org/header/5FID.pdb')

Then, you can access it as a local file as described in the previous sections::

    println pdb.text

The above one-liner prints the content of the remote PDB file. Previous sections provide code examples
showing how to stream or copy the content of files.

.. note:: Write and list operations are not supported for HTTP/S and FTP files.

.. _glob: http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
