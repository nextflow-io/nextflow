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

import com.google.common.hash.HashCode
import nextflow.exception.ProcessException
import nextflow.file.FileHolder
import nextflow.processor.TaskResolver
import nextflow.processor.TaskRun
import nextflow.util.CacheHelper
import nextflow.util.CacheHelper.HashMode
import spock.lang.Specification
import spock.lang.Timeout

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

    @Timeout(10)
    def 'should fail to restore an input when an upstream restore fails to start' () {
        given:
        def workDir = Files.createTempDirectory('test')
        def a = workDir.resolve('a.txt')
        def b = workDir.resolve('b.txt')
        def resolverA = Stub(TaskResolver) { workDirFor(_) >> { throw new IllegalStateException('boom') } }
        def resolverB = Stub(TaskResolver) { workDirFor(_) >> workDir.resolve('restore-b') }
        def taskA = Stub(TaskRun) { makeCopy() >> Stub(TaskRun) }
        def taskB = Stub(TaskRun) { makeCopy() >> Stub(TaskRun) { getInputFiles() >> [new FileHolder(a)] } }
        def producerA = new DeletedOutputs.Producer(taskA, HashCode.fromInt(1), workDir.resolve('task-a'), resolverA)
        def producerB = new DeletedOutputs.Producer(taskB, HashCode.fromInt(2), workDir.resolve('task-b'), resolverB)
        def deleted = new DeletedOutputs(null, null)
        deleted.files.put(a, new DeletedOutputs.DeletedFile(a, 1, null, producerA))
        deleted.files.put(b, new DeletedOutputs.DeletedFile(b, 1, null, producerB))
        def task = Stub(TaskRun) { getInputFiles() >> [new FileHolder(b)] }

        when:
        deleted.restoreInputs(task)
        then:
        def e = thrown(ProcessException)
        e.cause.message == 'boom'
        deleted.restoring.isEmpty()

        cleanup:
        workDir.deleteDir()
    }

}
