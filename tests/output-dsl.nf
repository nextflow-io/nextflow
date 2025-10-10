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
nextflow.preview.types = true

params.save_bam_bai = false

process fastqc {
  input:
  id: String

  output:
  tuple(id, file('*.fastqc.log'))

  script:
  """
  echo ${id} > ${id}.fastqc.log
  """
}

process align {
  input:
  id: String

  output:
  bam = tuple(id, file('*.bam'))
  bai = tuple(id, file('*.bai'))

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
  tuple(id, file('quant'))

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
  tuple(file('summary_report.html'), file('summary_data/data.json'), file('summary_data/fastqc.txt'))

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

  ch_logs = ch_samples
    .map { sample -> sample.fastqc }
    .collect()

  summary(ch_logs)

  publish:
  samples = ch_samples
  summary = summary.out
}

output {
  samples {
    path { sample ->
      sample.fastqc >> 'log/'
      sample.bam >> 'align/'
      sample.bai >> 'align/'
      sample.quant >> "quant/${sample.id}"
    }
    index {
      path 'samples.csv'
      header true
      sep ','
    }
  }

  summary {
    path { report, data_json, fastqc_txt ->
      report >> './'
      data_json >> './'
      fastqc_txt >> './'
    }
  }
}
