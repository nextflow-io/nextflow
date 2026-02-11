/*
 * Copyright 2019, Google Inc
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

package nextflow.cloud.google.util

import com.google.cloud.storage.BucketInfo
import com.google.cloud.storage.StorageOptions
import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Timeout

import java.nio.file.Files
import java.nio.file.Path

import nextflow.util.GlobalCacheLockManager
import nextflow.file.FileHelper
import spock.lang.Specification

/**
 * Tests for GlobalCacheLockManager with Google Cloud Storage file system
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@IgnoreIf({System.getenv('NXF_SMOKE')})
@Requires({System.getenv('GOOGLE_APPLICATION_CREDENTIALS')})
class GsGlobalCacheLockManagerTest extends Specification {

   String createBucket(){
       def bucketName = "nf-gs-test-${UUID.randomUUID()}"
       def bucket =  StorageOptions.newBuilder().build().getService().create(BucketInfo.of(bucketName))
       println("Bucket $bucketName created")
       return bucketName
   }

   String deleteBucket(String bucketName) {
       StorageOptions.newBuilder().build().getService().delete(bucketName)
       println("Bucket $bucketName deleted")
   }



    def 'should acquire lock on Google Cloud Storage'() {
        given:
        def bucket = createBucket()
        def hash = 'abcdef1234'
        def workDir = Path.of("gs://${bucket}/test-locks/")
        def manager = new GlobalCacheLockManager(workDir)


        // Ensure parent directory exists
        Files.createDirectories(workDir)

        // Clean up any existing lock file
        def lockPath = workDir.resolve(".${hash}.lock")

        when:
        def lock = manager.acquire(hash)

        then:
        lock != null
        Files.exists(lockPath)

        cleanup:
        FileHelper.deletePath(workDir.parent)
        if( bucket ) deleteBucket( bucket)

    }

    def 'should fail to acquire lock when already held on GCS'() {
        given:
        def bucket = createBucket()
        def hash = 'abcdef5678'
        def workDir = Path.of("gs://${bucket}/test-locks/")
        def manager = new GlobalCacheLockManager(workDir)


        // Ensure parent directory exists
        Files.createDirectories(workDir)

        // Clean up any existing lock file
        def lockPath = workDir.parent.resolve(".${hash}.lock")

        when:
        def lock1 = manager.acquire(hash)
        def lock2 = manager.acquire(hash)

        then:
        lock1 != null
        lock2 == null

        cleanup:
        FileHelper.deletePath(workDir)
        if( bucket ) deleteBucket( bucket)
    }

    def 'should release lock on GCS'() {
        given:
        def bucket = createBucket()
        def hash = 'abcdef9012'
        def workDir = Path.of("gs://${bucket}/test-locks/")
        def manager = new GlobalCacheLockManager(workDir)

        // Ensure parent directory exists
        Files.createDirectories(workDir)

        // Clean up any existing lock file
        def lockPath = workDir.resolve(".${hash}.lock")

        when:
        def lock = manager.acquire(hash)

        then:
        lock != null
        Files.exists(lockPath)

        when:
        lock.release()

        then:
        !Files.exists(lockPath)

        cleanup:
        FileHelper.deletePath(workDir)
        if( bucket ) deleteBucket( bucket)
    }

    def 'should allow acquiring lock after release on GCS'() {
        given:
        def bucket = createBucket()
        def hash = 'abcdef3456'
        def workDir = Path.of("gs://${bucket}/test-locks/")
        def manager = new GlobalCacheLockManager(workDir)

        // Ensure parent directory exists
        Files.createDirectories(workDir)

        // Clean up any existing lock file
        def lockPath = workDir.resolve(".${hash}.lock")

        when:
        def lock1 = manager.acquire(hash)
        lock1.release()
        def lock2 = manager.acquire(hash)

        then:
        lock1 != null
        lock2 != null

        cleanup:
        FileHelper.deletePath(workDir)
        if( bucket ) deleteBucket( bucket)
    }

    def 'should handle concurrent lock attempts on GCS'() {
        given:
        def bucket = createBucket()
        def hash = 'abcdef7890'
        def workDir = Path.of("gs://${bucket}/test-locks/")
        def manager = new GlobalCacheLockManager(workDir)
        def successCount = 0
        def failureCount = 0

        // Ensure parent directory exists
        Files.createDirectories(workDir)

        // Clean up any existing lock file
        def lockPath = workDir.resolve(".${hash}.lock")

        when:
        def lock = null
        def threads = (1..5).collect { int i ->
            Thread.start {
                def myLock = manager.acquire(hash)
                synchronized(this) {
                    if (myLock) {
                        successCount++
                        lock = myLock
                        sleep(100) // Hold lock briefly
                    } else {
                        failureCount++
                    }
                }
            }
        }
        threads*.join()

        then:
        successCount == 1
        failureCount == 4
        lock != null

        cleanup:
        FileHelper.deletePath(workDir)
        if( bucket ) deleteBucket( bucket)
    }
}
