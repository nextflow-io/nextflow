(your-first-script)=

# Your first script

## Run a pipeline

This script defines two processes. The first splits a string into 6-character chunks, writing each one to a file with the prefix `chunk_`, and the second receives these files and transforms their contents to uppercase letters. The resulting strings are emitted on the `result` channel and the final output is printed by the `view` operator. Copy the following example into your favorite text editor and save it to a file named `tutorial.nf`:

```{literalinclude} snippets/your-first-script.nf
:language: nextflow
```

Execute the script by entering the following command in your terminal:

```console
$ nextflow run tutorial.nf

N E X T F L O W  ~  version 23.10.0
executor >  local (3)
[69/c8ea4a] process > splitLetters   [100%] 1 of 1 ✔
[84/c8b7f1] process > convertToUpper [100%] 2 of 2 ✔
HELLO
WORLD!
```

:::{note}
For versions of Nextflow prior to `22.10.0`, you must explicitly enable DSL2 by adding `nextflow.enable.dsl=2` to the top of the script or by using the `-dsl2` command-line option.
:::

You can see that the first process is executed once, and the second twice. Finally the result string is printed.

It's worth noting that the process `convertToUpper` is executed in parallel, so there's no guarantee that the instance processing the first split (the chunk `Hello`) will be executed before the one processing the second split (the chunk `world!`). Thus, you may very likely see the final result printed in a different order:

```
WORLD!
HELLO
```

:::{tip}
The hexadecimal string, e.g. `22/7548fa`, is the unique hash of a task, and the prefix of the directory where the task is executed. You can inspect a task's files by changing to the directory `$PWD/work` and using this string to find the specific task directory.
:::

(getstarted-resume)=

## Modify and resume

Nextflow keeps track of all the processes executed in your pipeline. If you modify some parts of your script, only the processes that are actually changed will be re-executed. The execution of the processes that are not changed will be skipped and the cached result used instead. This helps a lot when testing or modifying part of your pipeline without having to re-execute it from scratch.

For the sake of this tutorial, modify the `convertToUpper` process in the previous example, replacing the process script with the string `rev $x`, like so:

```nextflow
process convertToUpper {
  input:
    path x
  output:
    stdout

  """
  rev $x
  """
}
```

Then save the file with the same name, and execute it by adding the `-resume` option to the command line:

```bash
nextflow run tutorial.nf -resume
```

It will print output similar to this:

```
N E X T F L O W  ~  version 23.10.0
executor >  local (2)
[69/c8ea4a] process > splitLetters   [100%] 1 of 1, cached: 1 ✔
[d0/e94f07] process > convertToUpper [100%] 2 of 2 ✔
olleH
!dlrow
```

You will see that the execution of the process `splitLetters` is actually skipped (the process ID is the same), and its results are retrieved from the cache. The second process is executed as expected, printing the reversed strings.

:::{tip}
The pipeline results are cached by default in the directory `$PWD/work`. Depending on your script, this folder can take up a lot of disk space. It's a good idea to clean this folder periodically, as long as you know you won't need to resume any pipeline runs.
:::

For more information, see the {ref}`cache-resume-page` page.

(getstarted-params)=

## Pipeline parameters

Pipeline parameters are simply declared by prepending to a variable name the prefix `params`, separated by dot character. Their value can be specified on the command line by prefixing the parameter name with a double dash character, i.e. `--paramName`

For the sake of this tutorial, you can try to execute the previous example specifying a different input string parameter, as shown below:

```bash
nextflow run tutorial.nf --str 'Bonjour le monde'
```

The string specified on the command line will override the default value of the parameter. The output will look like this:

```
N E X T F L O W  ~  version 23.10.0
executor >  local (4)
[8b/16e7d7] process > splitLetters   [100%] 1 of 1 ✔
[eb/729772] process > convertToUpper [100%] 3 of 3 ✔
m el r
edno
uojnoB
```

:::{versionchanged} 20.11.0-edge
Any `.` (dot) character in a parameter name is interpreted as the delimiter of a nested scope. For example, `--foo.bar Hello` will be interpreted as `params.foo.bar`. If you want to have a parameter name that contains a `.` (dot) character, escape it using the back-slash character, e.g. `--foo\.bar Hello`.
:::
