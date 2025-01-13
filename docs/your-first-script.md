(your-first-script)=

# Your first script

Once you have Nextflow installed, you’re ready to run your first script. This guide details how to run a basic pipeline on the command line. It includes:

- Running a pipeline
- Modifying and and resuming a pipeline
- Configuring a pipeline parameter

**Prerequisites**

You will need the following to get started:

- Nextflow 24.10.0 or later. See {ref}`install-page` for instructions to install or update your version of Nextflow.

:::{note}
This guide uses the second preview of publish targets, introduced in Nextflow 24.10.0. Publish targets are intended to replace the `publishDir` directive.
:::

## Run a pipeline

You will run a basic Nextflow pipeline that splits a string of text into two files, and then converts lowercase letters to uppercase letters. You can see the pipeline here:

```groovy
// Enable output targets
nextflow.preview.output = true

// Default parameter input
params.str = "Hello world!"

// splitString process
process splitString {
    
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
    tag "$y"

    input:
    path y

    output:
    path 'upper_chunk_*'

    script:
    """
    cat $y | tr '[a-z]' '[A-Z]' > upper_${y}
    """
}

// Entry workflow block
workflow {
    
    // Main section
    main:
    ch_str = Channel.of(params.str)     // Create a channel using parameter input
    ch_chunks = splitString(ch_str)     // Split string into chunks and create a named channel
    convertToUpper(ch_chunks.flatten()) // Convert lowercase letters to uppercase letters

    // Publish section
    publish:
    ch_chunks >> 'lower'          // chunk_* files will be published
    convertToUpper.out >> 'upper' // upper_* files will be published

}

// Output block
output {}
```

This script defines two processes:

- `splitString`: takes a string input, splits it into 6-character chunks, and writes the chunks to files with the prefix `chunk_`
- `convertToUpper`: takes files as input, transforms their contents to uppercase letters, and writes the uppercase strings to files with the prefix `upper_`

The `splitString` output is emitted as a single element. The `flatten` operator is used so that this element is split and each file is treated as a sole element.

The outputs from both processes are published in the `results` directory. To run your pipeline:

1. Create a new file named `main.nf` in your current directory
2. Copy and save the above pipeline to your new file
3. Run your pipeline using the following command:

    ```bash
    nextflow run main.nf
    ```

You will see output similar to the following:

```console
 N E X T F L O W   ~  version 24.10.3
Launching `main.nf` [cheeky_mandelbrot] DSL2 - revision: 261a76377c

WARN: WORKFLOW OUTPUT DSL IS A PREVIEW FEATURE - SYNTAX AND FUNCTIONALITY CAN CHANGE IN FUTURE RELEASES
executor >  local (3)
[e6/1f8da1] process > splitString (1)           [100%] 1 of 1 ✔
[fb/d2d170] process > convertToUpper (chunk_ab) [100%] 2 of 2 ✔
```

Nextflow creates a `work` directory to store files used during a pipeline run. Each execution of a processes is run as a separate task. The `splitString` process is run as one task and the `convertToUpper` process is run as two tasks. The hexadecimal string, for example `e6/1f8da1`, is the beginning of a unique hash. It is a prefix used to identify the task directory where the script was executed.

Pipeline outputs are published in the `results` directory.

:::{tip}
Run your pipeline with `-ansilog false` to see each task printed on a separate line:

```bash
nextflow run main.nf -ansilog false
```

You will see output similar to the following:

```console
 N E X T F L O W   ~  version 24.10.3
Launching `main.nf` [trusting_boyd] DSL2 - revision: 082867d4d6
WARN: WORKFLOW OUTPUT DSL IS A PREVIEW FEATURE - SYNTAX AND FUNCTIONALITY CAN CHANGE IN FUTURE RELEASES
[e6/1f8da1] Submitted process > splitString (1)
[77/36dbaf] Submitted process > convertToUpper (chunk_aa)
[fb/d2d170] Submitted process > convertToUpper (chunk_ab)
```

::: 

(getstarted-resume)=

## Modify and resume

Nextflow tracks task executions in a task cache, a key-value store of previously-executed tasks. The task cache is used in conjunction with the work directory to recover cached tasks. If you modify and re-run your pipeline, only the processes that are changed will be re-executed. Cached results will be used for tasks that don't change.

You can enable resumability using the `-resume` flag when running a pipeline. To modify and resume your pipeline:

1. Open `main.nf`
2. Replace the `convertToUpper` process with the following:

    ```groovy
    process convertToUpper {
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

    ```bash
    nextflow run main.nf -resume
    ```

You will see output similar to the following:

```console
 N E X T F L O W   ~  version 24.10.3

Launching `main.nf` [naughty_volta] DSL2 - revision: 082867d4d6

WARN: WORKFLOW OUTPUT DSL IS A PREVIEW FEATURE - SYNTAX AND FUNCTIONALITY CAN CHANGE IN FUTURE RELEASES
executor >  local (2)
[e6/1f8da1] process > splitString (1)           [100%] 1 of 1, cached: 1 ✔
[20/a896f5] process > convertToUpper (chunk_ab) [100%] 2 of 2 ✔
```

Nextflow skips the execution of the process `splitString` process and retrieves the results from the cache. The `convertToUpper` process is executed twice.

See {ref}`cache-resume-page` for more information about Nextflow cache and resume. 

(getstarted-params)=

## Pipeline parameters

Nextflow looks for configuration settings in multiple locations when launched. Configuration settings are merged and, if conflicts exist, override each other in a set order. 

Parameters are a type of configuration. They are declared by prepending a variable name to the prefix `params`, separated by dot character. Parameters can also be specified on the command line by prefixing the parameter name with a double dash character, for example, `--paramName`. Parameters specified on the command line override parameters specified in a pipeline scripts.

You can configure the `str` parameter in your pipeline. To modify your `str` parameter:

1. Run your pipeline using the following command:

    ```bash
    nextflow run main.nf --str 'Bonjour le monde'
    ```

You will see output similar to the following:

```console
 N E X T F L O W   ~  version 24.10.3

Launching `main.nf` [distracted_kalam] DSL2 - revision: 082867d4d6

WARN: WORKFLOW OUTPUT DSL IS A PREVIEW FEATURE - SYNTAX AND FUNCTIONALITY CAN CHANGE IN FUTURE RELEASES
executor >  local (4)
[55/a3a700] process > splitString (1)           [100%] 1 of 1 ✔
[f4/af5ddd] process > convertToUpper (chunk_ac) [100%] 3 of 3 ✔
```

The input string is now longer and the `splitString` process splits it into three chunks. The `convertToUpper` process is run three times.

## Next steps

See [training.nextflow.io](https://training.nextflow.io/) for Nextflow training modules.
