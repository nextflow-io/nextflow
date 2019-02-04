#!/usr/bin/env nextflow
  
  process io_read_write_1G {
    """
    dd if=/dev/zero of=/dev/null bs=1G count=1
    """
  }

  
  process io_read_write_256M {
    """
    dd if=/dev/zero of=/dev/null bs=256M count=1
    """
  }
