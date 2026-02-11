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
package nextflow.util

import java.nio.file.Files

import spock.lang.Specification

/**
 * Tests for GlobalCacheLockManager
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class GlobalCacheLockManagerTest extends Specification {

    def 'should acquire lock successfully'() {
        given:
        def tempDir = Files.createTempDirectory('test')
        def workDir = tempDir.resolve('work/')
        Files.createDirectories(workDir)
        def manager = new GlobalCacheLockManager(workDir)
        def hash = 'abcdef1234'

        when:
        def lock = manager.acquire(hash)

        then:
        lock != null
        Files.exists(tempDir.resolve('work/.abcdef1234.lock'))

        cleanup:
        lock?.release()
        tempDir?.deleteDir()
    }

    def 'should fail to acquire lock when already held'() {
        given:
        def tempDir = Files.createTempDirectory('test')
        def workDir = tempDir.resolve('work/')
        Files.createDirectories(workDir)
        def manager = new GlobalCacheLockManager(workDir)
        def hash = 'abcdef1234'

        when:
        def lock1 = manager.acquire(hash)
        def lock2 = manager.acquire(hash)

        then:
        lock1 != null
        lock2 == null

        cleanup:
        lock1?.release()
        tempDir?.deleteDir()
    }

    def 'should release lock successfully'() {
        given:
        def tempDir = Files.createTempDirectory('test')
        def workDir = tempDir.resolve('work/')
        Files.createDirectories(workDir)
        def manager = new GlobalCacheLockManager(workDir)
        def hash = 'abcdef1234'

        when:
        def lock = manager.acquire(hash)
        def lockPath = tempDir.resolve('work/.abcdef1234.lock')

        then:
        lock != null
        Files.exists(lockPath)

        when:
        lock.release()

        then:
        !Files.exists(lockPath)

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should allow acquiring lock after release'() {
        given:
        def tempDir = Files.createTempDirectory('test')
        def workDir = tempDir.resolve('work/')
        Files.createDirectories(workDir)
        def manager = new GlobalCacheLockManager(workDir)
        def hash = 'abcdef1234'

        when:
        def lock1 = manager.acquire(hash)
        lock1.release()
        def lock2 = manager.acquire(hash)

        then:
        lock1 != null
        lock2 != null

        cleanup:
        lock2?.release()
        tempDir?.deleteDir()
    }

    def 'should handle concurrent lock attempts'() {
        given:
        def tempDir = Files.createTempDirectory('test')
        def workDir = tempDir.resolve('work/')
        Files.createDirectories(workDir)
        def manager = new GlobalCacheLockManager(workDir)
        def hash = 'abcdef1234'
        def successCount = 0
        def failureCount = 0

        when:
        def threads = (1..10).collect {
            Thread.start {
                def lock = manager.acquire(hash)
                synchronized(this) {
                    if (lock) {
                        successCount++
                        sleep(10) // Hold lock briefly
                        lock.release()
                    } else {
                        failureCount++
                    }
                }
            }
        }
        threads*.join()

        then:
        successCount == 1
        failureCount == 9

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should be safe to release multiple times'() {
        given:
        def tempDir = Files.createTempDirectory('test')
        def workDir = tempDir.resolve('work/')
        Files.createDirectories(workDir)
        def manager = new GlobalCacheLockManager(workDir)
        def hash = 'abcdef1234'

        when:
        def lock = manager.acquire(hash)
        lock.release()
        lock.release() // Should not throw

        then:
        noExceptionThrown()

        cleanup:
        tempDir?.deleteDir()
    }
}
