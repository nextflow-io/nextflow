process stress_1cpu {
  """
  stress -c 1 -t 5
  """
}

process stress_2cpu {
  """
  stress -c 2 -t 10
  """
}

process stress_100mega {
  """
  stress -m 1 --vm-bytes 100000000 -t 10
  """
}

process stress_200mega {
  // note: mem usage is not aggregated 
  """
  stress -m 1 --vm-bytes 200000000 -t 5
  stress -m 1 --vm-bytes 100000000 -t 5
  """
}

process stress_300mega {
  // note: two parallel workers of 150MB => 300 MB
  """
  stress -m 2 --vm-bytes 150000000 -t 5
  """
}