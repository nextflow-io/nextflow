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
nextflow.preview.types = true

params {
  save_bam_bai = false
}

process fastqc {
  input:
  id: String

  output:
  record(id: id, fastqc: file('*.fastqc.log'))

  script:
  """
  echo ${id} > ${id}.fastqc.log
  """
}

process align {
  input:
  id: String

  output:
  record(id: id, bam: file('*.bam'), bai: file('*.bai'))

  script:
  """
  echo ${id} > ${id}.bam
  echo ${id} | rev > ${id}.bai
  """
}

process quant {
  input:
  id: String

  output:
  record(id: id, quant: file('quant'))

  script:
  '''
  mkdir quant
  touch quant/cmd_info.json
  touch quant/lib_format_counts.json
  touch quant/quant.sf
  '''
}

process summary {
  input:
  logs: Bag<Path>

  output:
  record(
    report: file('summary_report.html'), 
    data: file('summary_data/data.json'), 
    fastqc: file('summary_data/fastqc.txt')
  )

  script:
  '''
  touch summary_report.html
  mkdir summary_data
  touch summary_data/data.json
  touch summary_data/fastqc.txt
  '''
}

workflow {
  main:
  ids = channel.of('alpha', 'beta', 'delta')
  ch_fastqc = fastqc(ids)
  ch_align = align(ids)
  ch_quant = quant(ids)

  ch_samples = ch_fastqc
    .join(ch_align, by: 'id')
    .join(ch_quant, by: 'id')

  ch_logs = ch_samples
    .map { sample -> sample.fastqc }
    .collect()

  ch_summary = summary(ch_logs)

  publish:
  samples = ch_samples
  summary = ch_summary
}

output {
  samples {
    path { sample ->
      sample.fastqc >> 'log/'
      sample.bam >> (params.save_bam_bai ? 'align/' : null)
      sample.bai >> (params.save_bam_bai ? 'align/' : null)
      sample.quant >> "quant/${sample.id}"
    }
    index {
      path 'samples.csv'
      header true
      sep ','
    }
  }

  summary {
    path { r ->
      r.report >> './'
      r.data >> './'
      r.fastqc >> './'
    }
  }
}
