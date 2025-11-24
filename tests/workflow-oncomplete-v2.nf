#!/usr/bin/env nextflow
/*
 * Copyright 2013-2024, Seqera Labs
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

params.command = 'echo'

process sayHello {
  debug true

  input:
  val x

  script:
  """
  ${params.command} '${x} world!'
  """
}

workflow {
  main:
  channel.of('Bojour', 'Ciao', 'Hello', 'Hola', 'Γεια σου') | sayHello

  onComplete:
  log.info('Workflow completed successfully!')

  onError:
  log.error('Workflow failed.')
}
