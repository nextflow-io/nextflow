(your-first-script)=

# Your first script

This guide details fundamental skills to run a basic Nextflow pipeline. It includes:

- Running a pipeline
- Modifying and resuming a pipeline
- Configuring a pipeline parameter

<h3>Prerequisites</h3>

You will need the following to get started:

- Nextflow. See {ref}`install-page` for instructions to install or update your version of Nextflow.

## Run a pipeline

You will run a basic Nextflow pipeline that splits a string of text into two files and then converts lowercase letters to uppercase letters. You can see the pipeline here:

```{code-block} groovy
:class: copyable
// Default parameter input
params.str = "Hello world!"

// splitString process
process splitString {
    publishDir "results/lower"
    
    input:
    val x
    
    output:
    path 'chunk_*'

    script:
    """
    printf '${x}' | split -b 6 - chunk_
    """
}

// convertToUpper process
process convertToUpper {
    publishDir "results/upper"
    tag "$y"

    input:
    path y

    output:
    path 'upper_*'

    script:
    """
    cat $y | tr '[a-z]' '[A-Z]' > upper_${y}
    """
}

// Workflow block
workflow {
    ch_str = Channel.of(params.str)     // Create a channel using parameter input
    ch_chunks = splitString(ch_str)     // Split string into chunks and create a named channel
    convertToUpper(ch_chunks.flatten()) // Convert lowercase letters to uppercase letters
}
```

This script defines two processes:

- `splitString`: takes a string input, splits it into 6-character chunks, and writes the chunks to files with the prefix `chunk_`
- `convertToUpper`: takes files as input, transforms their contents to uppercase letters, and writes the uppercase strings to files with the prefix `upper_`

The `splitString` output is emitted as a single element. The `flatten` operator splits this combined element so that each file is treated as a sole element.

The outputs from both processes are published in subdirectories, that is, `lower` and `upper`, in the `results` directory.

To run your pipeline:

1. Create a new file named `main.nf` in your current directory
2. Copy and save the above pipeline to your new file
3. Run your pipeline using the following command:

    ```{code-block}
    :class: copyable
    nextflow run main.nf
    ```

You will see output similar to the following:

```console
 N E X T F L O W   ~  version 24.10.3

Launching `main.nf` [big_wegener] DSL2 - revision: 13a41a8946

executor >  local (3)
[82/457482] splitString (1)           | 1 of 1 ✔
[2f/056a98] convertToUpper (chunk_aa) | 2 of 2 ✔
```

Nextflow creates a `work` directory to store files used during a pipeline run. Each execution of a process is run as a separate task. The `splitString` process is run as one task and the `convertToUpper` process is run as two tasks. The hexadecimal string, for example, `82/457482`, is the beginning of a unique hash. It is a prefix used to identify the task directory where the script was executed.

:::{tip}
Run your pipeline with `-ansi-log false` to see each task printed on a separate line:

```{code-block} bash
:class: copyable
nextflow run main.nf -ansi-log false
```

You will see output similar to the following:

```console
N E X T F L O W  ~  version 24.10.3
Launching `main.nf` [peaceful_watson] DSL2 - revision: 13a41a8946
[43/f1f8b5] Submitted process > splitString (1)
[a2/5aa4b1] Submitted process > convertToUpper (chunk_ab)
[30/ba7de0] Submitted process > convertToUpper (chunk_aa)
```

::: 

(getstarted-resume)=

## Modify and resume

Nextflow tracks task executions in a task cache, a key-value store of previously executed tasks. The task cache is used in conjunction with the work directory to recover cached tasks. If you modify and resume your pipeline, only the processes that are changed will be re-executed. The cached results will be used for tasks that don't change.

You can enable resumability using the `-resume` flag when running a pipeline. To modify and resume your pipeline:

1. Open `main.nf`
2. Replace the `convertToUpper` process with the following:

    ```{code-block} groovy
    :class: copyable
    process convertToUpper {
        publishDir "results/upper"
        tag "$y"

        input:
        path y

        output:
        path 'upper_*'

        script:
        """
        rev $y > upper_${y}
        """
    }
    ```

3. Save your changes
4. Run your updated pipeline using the following command:

    ```{code-block} bash
    :class: copyable
    nextflow run main.nf -resume
    ```

You will see output similar to the following:

```console
 N E X T F L O W   ~  version 24.10.3

Launching `main.nf` [furious_curie] DSL2 - revision: 5490f13c43

executor >  local (2)
[82/457482] splitString (1)           | 1 of 1, cached: 1 ✔
[02/9db40b] convertToUpper (chunk_aa) | 2 of 2 ✔
```

Nextflow skips the execution of the `splitString` process and retrieves the results from the cache. The `convertToUpper` process is executed twice.

See {ref}`cache-resume-page` for more information about Nextflow cache and resume functionality. 

(getstarted-params)=

## Pipeline parameters

Parameters are used to control the inputs to a pipeline. They are declared by prepending a variable name to the prefix `params`, separated by dot character. Parameters can be specified on the command line by prefixing the parameter name with a double dash character, for example, `--paramName`. Parameters specified on the command line override parameters specified in a main script.

You can configure the `str` parameter in your pipeline. To modify your `str` parameter:

1. Run your pipeline using the following command:

    ```{code-block} bash
    :class: copyable
    nextflow run main.nf --str 'Bonjour le monde'
    ```

You will see output similar to the following:

```console
 N E X T F L O W   ~  version 24.10.3

Launching `main.nf` [distracted_kalam] DSL2 - revision: 082867d4d6

executor >  local (4)
[55/a3a700] process > splitString (1)           [100%] 1 of 1 ✔
[f4/af5ddd] process > convertToUpper (chunk_ac) [100%] 3 of 3 ✔
```

The input string is now longer and the `splitString` process splits it into three chunks. The `convertToUpper` process is run three times.

See {ref}`cli-params` for more information about modifying pipeline parameters.

<h2>Next steps</h2>

Your first script is a brief introduction to running pipelines, modifying and resuming pipelines, and pipeline parameters. See [training.nextflow.io](https://training.nextflow.io/) for further Nextflow training modules.
