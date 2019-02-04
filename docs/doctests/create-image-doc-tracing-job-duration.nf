process stress_1cpu {
  """
  stress -c 1 -t 5
  """
}

process stress_2cpu {
  cpus 2
  """
  stress -c 2 -t 10
  """
}

process stress_100M_RAM {
  memory 110.MB
  """
  stress -m 1 --vm-bytes 100000000 -t 10
  """
}

process stress_1G_RAM {
  memory 1.1.GB
  """
  stress -m 1 --vm-bytes 1000000000 -t 5
  """
}



process io_read_write_1G {
  """
  dd if=/dev/zero of=1giga.img bs=1G count=1
  """
}

process io_read_write_200M {
  cpus 2
  """
  dd if=/dev/zero of=100mega.img bs=100M count=1
  dd if=/dev/zero of=100mega.img bs=100M count=1
  """
}

