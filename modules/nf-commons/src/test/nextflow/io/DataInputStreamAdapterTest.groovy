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

package nextflow.io

import java.nio.file.Files

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataInputStreamAdapterTest extends Specification {

    def testDataInput() {

        given:
        def file = Files.createTempFile('test',null).toFile()
        file.text = 'a\nbb\nccc\n'
        def data = new RandomAccessFile(file,'rw')

        when:
        def x = new BufferedReader(new InputStreamReader(new DataInputStreamAdapter(data)))
        then:
        x.readLine() == 'a'
        x.readLine() == 'bb'
        x.readLine() == 'ccc'
        x.readLine() == null

        cleanup:
        file?.delete()
    }

}
