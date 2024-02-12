/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.cloud.google.file

import java.nio.file.Path

import nextflow.Global
import nextflow.Session
import nextflow.cloud.google.nio.GsPath
import nextflow.file.FileHelper
import nextflow.util.KryoHelper
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GsPathSerializerTest extends Specification {

    def 'should serialize a google cloud path'() {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [google:[project:'foo', region:'x']]
        }

        when:
        def path = FileHelper.asPath("gs://my-seq/data/ggal/sample.fq")
        def buffer = KryoHelper.serialize(path)
        def copy = (Path)KryoHelper.deserialize(buffer)
        then:
        copy instanceof GsPath
        copy.toUri() == URI.create("gs://my-seq/data/ggal/sample.fq")
        copy.toUriString() == "gs://my-seq/data/ggal/sample.fq"
    }
}
