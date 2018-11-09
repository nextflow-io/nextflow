#!/usr/bin/env nextflow
/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * Verify that variable defined in the configuration file
 * context can be accessed in the process script
 *
 * author Emilio Palumbo
 */

echo true

process foo {
 
  echo true
 
  script:
  t = task.ext.out.join(',')
  """
  echo foo ${t}
  """
 
}

process bar {
 
  echo true
 
  shell:
  t = task.ext.out.join(',')
  
  '''
  echo bar !{t}
  '''
 
}