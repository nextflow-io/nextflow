#!/usr/bin/env nextflow

/*
 * Demonstration of the unified `package` directive with one process per
 * provider plugin: conda, pixi, uv, nix, and guix.
 *
 * Run with:
 *   nextflow run tests/nf-test/examples/package-providers-demo.nf -c tests/nf-test/examples/package-providers-demo.config
 *
 * Each provider is supplied by its own `nf-<provider>` plugin and resolves the
 * declared dependency for that single process only, aligning with Nextflow's
 * process-specific dependency model.
 *
 * NOTE: the `package` directive requires the v2 syntax parser:
 *   export NXF_SYNTAX_PARSER=v2
 */

nextflow.enable.dsl=2

// conda provider (nf-conda) -- also serves the conda/mamba/micromamba family
process useConda {
    package "samtools=1.17", provider: "conda"

    output:
    stdout

    script:
    """
    samtools --version | head -1
    """
}

// pixi provider (nf-pixi)
process usePixi {
    package "samtools=1.17", provider: "pixi"

    output:
    stdout

    script:
    """
    samtools --version | head -1
    """
}

// uv provider (nf-uv) -- Python packages via uv virtual environments
process useUv {
    package "numpy", provider: "uv"

    output:
    stdout

    script:
    """
    python -c "import numpy; print('numpy', numpy.__version__)"
    """
}

// nix provider (nf-nix) -- resolves bare names against the `nixpkgs` flake
process useNix {
    package "hello", provider: "nix"

    output:
    stdout

    script:
    """
    hello
    """
}

// guix provider (nf-guix)
process useGuix {
    package "hello", provider: "guix"

    output:
    stdout

    script:
    """
    hello
    """
}

// pak provider (nf-pak) -- R packages via pak
process usePak {
    package "jsonlite", provider: "pak"

    output:
    stdout

    script:
    """
    Rscript -e 'cat("jsonlite", as.character(packageVersion("jsonlite")), "\\n")'
    """
}

// install2.r provider (nf-install2r) -- R packages via littler
process useInstall2r {
    package "jsonlite", provider: "install2r"

    output:
    stdout

    script:
    """
    Rscript -e 'cat("jsonlite", as.character(packageVersion("jsonlite")), "\\n")'
    """
}

workflow {
    useConda()      | view { "conda:     ${it.trim()}" }
    usePixi()       | view { "pixi:      ${it.trim()}" }
    useUv()         | view { "uv:        ${it.trim()}" }
    useNix()        | view { "nix:       ${it.trim()}" }
    useGuix()       | view { "guix:      ${it.trim()}" }
    usePak()        | view { "pak:       ${it.trim()}" }
    useInstall2r()  | view { "install2r: ${it.trim()}" }
}
