#!/usr/bin/env nextflow

/*
 * "Hello <provider>" demonstration for every package-manager provider.
 *
 * Each process prints `Hello <provider>` after resolving a trivial dependency
 * through that provider's `nf-<provider>` plugin. Because the package managers
 * themselves are usually not installed on the host, each process runs inside a
 * base container image that ships the corresponding tool (conda, pixi, uv, nix,
 * guix, R+pak, R+install2.r). Nextflow then drives the provider *inside* that
 * container to build the per-process environment.
 *
 * Run with:
 *   nextflow run tests/nf-test/examples/package-providers-hello.nf -c tests/nf-test/examples/package-providers-hello.config
 *
 * Notes:
 *   - conda / pixi / uv / pak / install2r containerise cleanly.
 *   - nix and guix are best-effort: their store/daemon model is awkward inside a
 *     transient task container, so the cached profile does not survive across
 *     runs (fine for a one-shot hello).
 *
 * NOTE: the `package` directive requires the v2 syntax parser:
 *   export NXF_SYNTAX_PARSER=v2
 */

nextflow.enable.dsl=2

process helloConda {
    container 'continuumio/miniconda3:latest'
    package "cowsay", provider: "conda"

    output:
    stdout

    script:
    """
    cowsay "Hello conda"
    """
}

process helloPixi {
    container 'ghcr.io/prefix-dev/pixi:latest'
    package "cowpy", provider: "pixi"

    output:
    stdout

    script:
    """
    cowpy "Hello pixi"
    """
}

process helloUv {
    container 'ghcr.io/astral-sh/uv:python3.12-bookworm'
    package "cowsay", provider: "uv"

    output:
    stdout

    script:
    """
    python -c "import cowsay; print('Hello uv')"
    """
}

process helloNix {
    container 'nixos/nix:latest'
    package "hello", provider: "nix"

    output:
    stdout

    script:
    """
    hello >/dev/null && echo "Hello nix"
    """
}

process helloGuix {
    container 'metacall/guix:latest'
    package "hello", provider: "guix"

    output:
    stdout

    script:
    """
    hello >/dev/null && echo "Hello guix"
    """
}

process helloPak {
    container 'rocker/r-ver:latest'
    package "jsonlite", provider: "pak"

    output:
    stdout

    script:
    """
    Rscript -e 'library(jsonlite); cat("Hello pak\\n")'
    """
}

process helloInstall2r {
    container 'rocker/r-ver:latest'
    package "jsonlite", provider: "install2r"

    output:
    stdout

    script:
    """
    Rscript -e 'library(jsonlite); cat("Hello install2r\\n")'
    """
}

workflow {
    helloConda()     | view { it.trim() }
    helloPixi()      | view { it.trim() }
    helloUv()        | view { it.trim() }
    helloNix()       | view { it.trim() }
    helloGuix()      | view { it.trim() }
    helloPak()       | view { it.trim() }
    helloInstall2r() | view { it.trim() }
}
