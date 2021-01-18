/*
 * Copyright 2020, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.cloud.azure.batch

import java.nio.file.Paths

import nextflow.processor.TaskBean
import spock.lang.Specification
import test.TestHelper

// TODO: Adapt the s3 related tests for Azure batch/rclone combination
class AzFileCopyStrategyTest extends Specification {

    def 'should strip out file/folder name from target blob container path'() {
        given:
        def OUTPUTS = ["outputs_*", "final_folder"]
        def TARGET = Paths.get('/data/results')
        def FILE = Paths.get('/some/data/nobel_prize_results.gz')
        def EXIT = Paths.get('/some/path/.exitcode')
        def RUN = Paths.get('/some/data/.command.run')
        def SAS_URL = "sas_url"
        def copy = new AzFileCopyStrategy(Mock(TaskBean), SAS_URL)
        expect:
        copy.touchFile(RUN) == "echo start | aws s3 cp --only-show-errors - s3://some/data/.command.run"
        copy.copyFile("nobel_prize_results.gz", Paths.get("/some/data/nobel_prize_results.gz")) == "nxf_s3_upload nobel_prize_results.gz s3://some/data"
        copy.exitFile(EXIT) == "| aws s3 cp --only-show-errors - s3://some/path/.exitcode || true"
        copy.stageInputFile(FILE, 'foo.txt') == """
                                    downloads+=("nxf_s3_download s3://some/data/nobel_prize_results.gz foo.txt")
                                    """
                .stripIndent().trim()
        copy.getUnstageOutputFilesScript(OUTPUTS, TARGET) == '''
                                        uploads=()
                                        IFS=$'\\n'
                                        for name in $(eval "ls -1d outputs_* final_folder" | sort | uniq); do
                                            uploads+=("nxf_s3_upload '$name' s3://data/results")
                                        done
                                        unset IFS
                                        nxf_parallel "${uploads[@]}"
                                        '''
                .stripIndent().leftTrim()
    }

}
