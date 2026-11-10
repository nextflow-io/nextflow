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

package nextflow.processor

import java.nio.file.Files
import java.nio.file.Path

import nextflow.file.FileHolder
import nextflow.util.CacheHelper
import spock.lang.Specification
import spock.lang.TempDir

class TaskStageAsFilesTest extends Specification {

    @TempDir
    Path folder

    private static TaskStageAsFiles staged(String pattern, List<Path> values, String inputName = null) {
        final files = TaskInputResolver.expandWildcards(pattern, values.collect { new FileHolder(it) })
        return new TaskStageAsFiles(pattern, files, inputName)
    }

    private static String hash(List<TaskStageAsFiles> staged, CacheHelper.HashMode mode = CacheHelper.HashMode.STANDARD) {
        return CacheHelper.hasher(staged, mode).hash().toString()
    }

    private Path file(String name, String text) {
        final path = folder.resolve(name)
        Files.createDirectories(path.parent)
        return Files.write(path, text.bytes)
    }

    def 'should hash the pattern and files of each stageAs directive' () {
        given:
        def a = file('a.txt', 'a')
        def b = file('b.txt', 'b')

        expect:
        hash([staged('one.txt', [a]), staged('two.txt', [b])]) == hash([staged('one.txt', [a]), staged('two.txt', [b])])
        and: 'renaming a file changes the hash'
        hash([staged('one.txt', [a]), staged('two.txt', [b])]) != hash([staged('one.txt', [a]), staged('three.txt', [b])])
        and: 'staging under the default name is distinct from not staging'
        hash([staged('*', [a])]) != hash([])
        and: 'swapping the names between files changes the hash'
        hash([staged('one.txt', [a]), staged('two.txt', [b])]) != hash([staged('one.txt', [b]), staged('two.txt', [a])])
        and: 'selecting a different file changes the hash'
        hash([staged('*', [a])]) != hash([staged('*', [b])])
        and: 'the same names in a different iteration order do not change the hash'
        hash([staged('dir/*', [a, b])]) == hash([staged('dir/*', [b, a])])
        and: 'numbered names that follow a different iteration order change the hash'
        hash([staged('dir?/*', [a, b])]) != hash([staged('dir?/*', [b, a])])

        when:
        def before = hash([staged('one.txt', [a])])
        a.text = 'edited with a different size'
        then: 'editing a file changes the hash'
        before != hash([staged('one.txt', [a])])
    }

    def 'should hash staged files according to the cache mode' () {
        given:
        def a1 = file('d1/a.txt', 'same')
        def a2 = file('d2/a.txt', 'same')

        expect: 'the same content at another location has the same deep hash'
        hash([staged('one.txt', [a1])], CacheHelper.HashMode.DEEP) == hash([staged('one.txt', [a2])], CacheHelper.HashMode.DEEP)
        and: 'but not the same standard hash, which includes the path'
        hash([staged('one.txt', [a1])]) != hash([staged('one.txt', [a2])])
    }

    def 'should identify the files of a task input by the input name' () {
        given:
        def a = file('a.txt', 'a')

        when:
        def before = hash([staged('one.txt', [a], 'fasta')])
        a.text = 'edited with a different size'
        then: 'the files are not hashed again, since the input hash covers them'
        before == hash([staged('one.txt', [a], 'fasta')])
        and: 'the pattern and the input are bound together'
        before != hash([staged('two.txt', [a], 'fasta')])
        before != hash([staged('one.txt', [a], 'query')])
        before != hash([staged('one.txt', [a])])
    }

    def 'should hash to a stable value' () {
        given:
        def a = file('a.txt', 'a')

        expect:
        hash([staged('query.fa', [Path.of('/data/a.fa')], 'fasta')]) == '8dc5aeeaeed8d5d609a96c60d8fbec21'
        hash([staged('dir?/*', [a])], CacheHelper.HashMode.DEEP) == '3c1947c9cffaed4f2ac969b5a5c327fd'
    }

}
