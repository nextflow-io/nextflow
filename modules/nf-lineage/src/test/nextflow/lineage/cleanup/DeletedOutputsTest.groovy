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


package nextflow.lineage.cleanup

import java.nio.file.Files
import java.nio.file.attribute.BasicFileAttributes

import nextflow.file.FileHolder
import nextflow.util.CacheHelper
import nextflow.util.CacheHelper.HashMode
import spock.lang.Specification

class DeletedOutputsTest extends Specification {

    def 'should hash a deleted file the same as the original file' () {
        given:
        def file = Files.createTempFile('test', '.txt')
        file.text = 'hello'
        def attrs = Files.readAttributes(file, BasicFileAttributes)
        def deleted = new DeletedOutputs.DeletedFile(file, attrs.size(), attrs.lastModifiedTime().toMillis(), null)

        expect:
        CacheHelper.hasher([new FileHolder(file)], mode).hash() == CacheHelper.hasher([deleted], mode).hash()

        cleanup:
        Files.deleteIfExists(file)

        where:
        mode << [HashMode.STANDARD, HashMode.LENIENT]
    }

}
