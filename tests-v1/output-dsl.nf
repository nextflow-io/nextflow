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
nextflow.preview.output = true

params.save_bam_bai = false

process fastqc {
  input:
  val id

  output:
  tuple val(id), path('*.fastqc.log')

  script:
  """
  echo ${id} > ${id}.fastqc.log
  """
}

process align {
  input:
  val id

  output:
  tuple val(id), path('*.bam')
  tuple val(id), path('*.bai')

  script:
  """
  echo ${id} > ${id}.bam
  echo ${id} | rev > ${id}.bai
  """
}

process quant {
  input:
  val id

  output:
  tuple val(id), path('quant')

  script:
  '''
  mkdir quant
  touch quant/cmd_info.json
  touch quant/lib_format_counts.json
  touch quant/quant.sf
  '''
}

workflow {
  main:
  ids = Channel.of('alpha', 'beta', 'delta')
  ch_fastqc = fastqc(ids)
  (ch_bam, ch_bai) = align(ids)
  ch_quant = quant(ids)

  ch_samples = ch_fastqc
    .join(ch_bam)
    .join(ch_bai)
    .join(ch_quant)
    .map { id, fastqc, bam, bai, quant ->
      [
        id: id,
        fastqc: fastqc,
        bam: params.save_bam_bai ? bam : null,
        bai: params.save_bam_bai ? bai : null,
        quant: quant
      ]
    }

  publish:
  ch_samples >> 'samples'
}

output {
  samples {
    path { sample ->
      def dirs = [
        'bam': 'align',
        'bai': 'align',
        'log': 'fastqc'
      ]
      return { filename ->
        def ext = filename.tokenize('.').last()
        def dir = dirs[ext]
        dir != null
          ? "${dir}/${filename}"
          : "${filename}/${sample.id}"
      }
    }
    index {
      path 'samples.csv'
      header true
      sep ','
    }
  }
}
