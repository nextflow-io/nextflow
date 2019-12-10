/*
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

package nextflow.cloud.google.lifesciences

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GoogleLifeSciencesSubmitRequestTest extends Specification {

    def 'should create wrapper scripts' () {
        when:
        def dir = CloudStorageFileSystem.forBucket("my-bucket").getPath('/work/dir')
        def req = new GoogleLifeSciencesSubmitRequest(workDir: dir)
        then:
        req.stagingScript ==
                    'set -x; { cd /work/dir; gsutil -m -q cp gs://my-bucket/work/dir/.command.run .; bash .command.run nxf_stage; } 2>&1 > /work/dir/.command.log'
        and:
        req.mainScript ==
                    '{ cd /work/dir; bash .command.run; } 2>&1 | tee -a /work/dir/.command.log'
        and:
        req.unstagingScript ==
                    'set -x; { cd /work/dir; bash .command.run nxf_unstage; } 2>&1 | tee -a /work/dir/.command.log; gsutil -m -q cp -R /work/dir/.command.log gs://my-bucket/work/dir/.command.log || true; [[ $GOOGLE_LAST_EXIT_STATUS -gt 0 || $NXF_DEBUG -gt 0 ]] && { gsutil -m -q cp -R /google/ gs://my-bucket/work/dir; } || rm -rf /work/dir'
    }

}
