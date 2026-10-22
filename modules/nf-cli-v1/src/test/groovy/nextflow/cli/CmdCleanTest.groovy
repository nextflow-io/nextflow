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

package nextflow.cli

import java.nio.file.Files
import java.nio.file.Path

import com.google.common.hash.HashCode
import nextflow.cache.CacheDB
import nextflow.cache.DefaultCacheStore
import nextflow.trace.TraceRecord
import nextflow.util.KryoHelper
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Lorenz Gerber <lorenzottogerber@gmail.com>
 */
class CmdCleanTest extends Specification {

    def 'empty folder' () {

        given:
        def folder = Files.createTempDirectory('test')
        def cleaner = new CmdClean()

        when:
        def result = cleaner.deleteFolder(folder,false)


        then:
        assert(result == true)
        assert(!folder.exists())

        cleanup:
        def exists = folder?.exists()
        if (exists)
            folder?.deleteDir()

    }

    def 'empty folder, keep logs' () {

        given:
        def folder = Files.createTempDirectory('test')
        def cleaner = new CmdClean()

        when:
        def result = cleaner.deleteFolder(folder,true)

        then:
        assert(result == true)
        assert(folder.exists())

        cleanup:
        def exists = folder?.exists()
        if (exists)
            folder?.deleteDir()

    }


    def 'non-empty folder [.command.file name] keep logs' () {

        given:
        def folder = Files.createTempDirectory('test')
        def file = Files.createFile(folder.resolve('.command.test'))
        def cleaner = new CmdClean()

        when:
        def result = cleaner.deleteFolder(file,true)

        then:
        assert(result == true)
        assert(file.getName().equals('.command.test'))

        cleanup:
        def fileExists = file?.exists()
        if (fileExists)
            file.delete()
        def folderExists = folder?.exists()
        if (folderExists)
            folder?.deleteDir()
    }

    def 'should not report a work dir as removed when it cannot be visited' () {
        given:
        def cleaner = new CmdClean(quiet: true)

        when:
        // a missing work dir reaches visitFileFailed, which must record the failure
        def result = cleaner.deleteFolder(Path.of('/missing/work/dir'), false)

        then:
        !result
    }

    private static void putEntry(CacheDB db, HashCode key, Map trace, int refCount=1) {
        db.@store.putEntry(key, KryoHelper.serialize([new TraceRecord(trace).serialize(), null, refCount]))
        db.@store.writeIndex(key, false)
    }

    def 'should remove the unlisted previous attempts of a completed task' () {
        given:
        def folder = Files.createTempDirectory('test')
        def name = '0123456789abcdef0123456789abcd'
        def dirs = ['', '-2', '-3', '-4', '-5'].collect { Files.createDirectories(folder.resolve("work/ab/${name}${it}")) }
        and:
        def key = HashCode.fromInt(100)
        def db = new CacheDB(new DefaultCacheStore(UUID.randomUUID(), 'run_1', folder)).open()
        putEntry(db, key, [status: 'COMPLETED', workdir: dirs[3].toString()], refCount)
        and: 'the failed attempt 2 is listed by a run'
        putEntry(db, CacheDB.attemptKey(key, 2), [status: 'FAILED', workdir: dirs[1].toString()])
        and:
        def cleaner = new CmdClean(quiet: true, force: true)
        cleaner.@currentCacheDb = db

        when:
        cleaner.removeRecord(key, db.getTraceRecord(key), refCount)

        then: 'attempt 1 and 3 are removed along with attempt 4, the others are kept'
        dirs.collect { Files.exists(it) } == EXISTS

        cleanup:
        db?.close()
        folder?.deleteDir()

        where:
        refCount | EXISTS
        1        | [false, true, false, false, true]
        2        | [true, true, true, true, true]
    }

    def 'should not remove the previous attempts of a failed task' () {
        given:
        def folder = Files.createTempDirectory('test')
        def name = '0123456789abcdef0123456789abcd'
        def dirs = ['', '-2'].collect { Files.createDirectories(folder.resolve("work/ab/${name}${it}")) }
        and:
        def key = CacheDB.attemptKey(HashCode.fromInt(100), 2)
        def db = new CacheDB(new DefaultCacheStore(UUID.randomUUID(), 'run_1', folder)).open()
        putEntry(db, key, [status: 'FAILED', workdir: dirs[1].toString()])
        and:
        def cleaner = new CmdClean(quiet: true, force: true)
        cleaner.@currentCacheDb = db

        when:
        cleaner.removeRecord(key, db.getTraceRecord(key), 1)

        then:
        Files.exists(dirs[0])
        !Files.exists(dirs[1])

        cleanup:
        db?.close()
        folder?.deleteDir()
    }

    def 'should not look for previous attempts of a first attempt' () {
        given:
        def cleaner = new CmdClean()

        expect:
        cleaner.previousAttempts(HashCode.fromInt(100), Path.of('/work/ab/0123456789abcdef0123456789abcd')) == []
    }

}
