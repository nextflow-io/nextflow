#!/usr/bin/env nextflow
/*
 * Copyright 2013-2026, Seqera Labs
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

params.data_link_aws = 'seqera://seqeralabs/automated-testing/data-links/aws/1000genomes/README.alignment_data'
params.data_link_google = 'seqera://seqeralabs/automated-testing/data-links/google/nf-core-gcpmegatests/test-data/rnaseq/README'
params.data_link_azure = 'seqera://seqeralabs/automated-testing/data-links/azure/seqeralabs-showcase/results/pipeline_info/nf_core_rnaseq_software_mqc_versions.yml'

process TEST {
    input:
        path(file)
    output:
        stdout
    script:
    """
    ls $file
    """
}

workflow {
    ch = channel.fromList([file(params.data_link_aws), file(params.data_link_google), file(params.data_link_azure)])
    TEST(ch).view()
}
