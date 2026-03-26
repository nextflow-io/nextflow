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

package nextflow.cloud.aws.util

import nextflow.cloud.aws.nio.S3Path
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class S3PathFactoryTest extends Specification {

    def 'should parse s3 paths' () {

        when:
        def path = S3PathFactory.parse(S3_PATH)
        then:
        path instanceof S3Path
        with(path as S3Path) {
            getBucket() == BUCKET
            getKey() == KEY
        }

        when:
        def str = S3PathFactory.getUriString(path)
        then:
        str == S3_PATH


        where:
        S3_PATH                                                 | BUCKET        | KEY
        's3://cbcrg-eu/raw/x_r1.fq'                             | 'cbcrg-eu'    | 'raw/x_r1.fq'
        's3://cbcrg-eu/raw/**_R1*{fastq,fq,fastq.gz,fq.gz}'     | 'cbcrg-eu'    | 'raw/**_R1*{fastq,fq,fastq.gz,fq.gz}'

    }

    def 'should ignore double slashes' () {
        when:
        def path = S3PathFactory.parse('s3://cbcrg-eu/raw//x_r1.fq' )
        then:
        S3PathFactory.getUriString(path) == 's3://cbcrg-eu/raw/x_r1.fq'
    }
}
