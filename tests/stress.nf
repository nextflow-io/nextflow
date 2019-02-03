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

process stress_100mega {
  memory 150.MB
  """
  stress -m 1 --vm-bytes 100000000 -t 10
  """
}

process stress_200mega {
  memory 250.MB
  // note: mem usage should not be aggregated 
  """
  stress -m 1 --vm-bytes 200000000 -t 5
  stress -m 1 --vm-bytes 100000000 -t 5
  """
}

process stress_300mega {
  memory 350.MB
  cpus 2
  // note: two parallel workers of 150MB => 300 MB
  """
  stress -m 2 --vm-bytes 150000000 -t 5
  """
}

process io_write_100mega {
  """
  write.pl file.txt 104857600
  """
}

process io_write_200mega {
  cpus 2
  """
  write.pl file1.txt 104857600 &
  pid=\$!
  write.pl file2.txt 104857600
  read.pl file2.txt
  """
}