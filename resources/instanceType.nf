
task {
  echo true
  instanceType 'dx_m2.4xlarge'

  """
  echo Num of procs: \$(nproc)
  cat /proc/meminfo | grep Mem
  """

}

