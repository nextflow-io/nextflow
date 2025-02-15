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

package nextflow.util

import nextflow.cloud.aws.util.S3PathFactory
import spock.lang.Specification
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class S3PathSerializerTest extends Specification {

    def 'should serialise s3 path' () {
        when:
        def path = S3PathFactory.parse('s3://mybucket/file.txt')
        def buffer = KryoHelper.serialize(path)
        then:
        KryoHelper.deserialize(buffer).getClass().getName() == 'nextflow.cloud.aws.nio.S3Path'
        KryoHelper.deserialize(buffer) == S3PathFactory.parse('s3://mybucket/file.txt')
    }

    def 'should serialise s3 path with spaces' () {
        when:
        def path = S3PathFactory.parse('s3://mybucket/file with spaces.txt')
        def buffer = KryoHelper.serialize(path)
        then:
        KryoHelper.deserialize(buffer).getClass().getName() == 'nextflow.cloud.aws.nio.S3Path'
        KryoHelper.deserialize(buffer) == S3PathFactory.parse('s3://mybucket/file with spaces.txt')
    }

}
