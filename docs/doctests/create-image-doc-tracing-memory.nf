#!/usr/bin/env nextflow

process foo {

    memory '1.5 GB'

    """
    memory_vmem_1GiB_ram_0Gib
    """

}

process bar {

    memory '1.5 GB'

    """
    memory_vmem_1GiB_ram_1Gib
    """

}
