/*
 * Copyright 2021, Microsoft Corp
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

package nextflow.cloud.azure.nio

import com.azure.storage.blob.BlobServiceClient
import com.azure.storage.blob.BlobServiceClientBuilder
import com.azure.storage.common.StorageSharedKeyCredential
import nextflow.Global
import nextflow.Session
import spock.lang.Shared
import spock.lang.Timeout

import java.nio.file.Files
import java.nio.file.Path

import nextflow.util.GlobalCacheLockManager
import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Specification

/**
 * Tests for GlobalCacheLockManager with Azure Blob Storage file system
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */

@Timeout(30)
@IgnoreIf({System.getenv('NXF_SMOKE')})
@Requires({System.getenv('AZURE_STORAGE_ACCOUNT_NAME') && System.getenv('AZURE_STORAGE_ACCOUNT_KEY')})
class AzGlobalCacheLockManagerTest extends Specification implements AzBaseSpec {

    @Shared
    BlobServiceClient storageClient

    def setupSpec() {
        def bucket = createBucket()
        def accountKey = System.getenv('AZURE_STORAGE_ACCOUNT_KEY')
        def credential = new StorageSharedKeyCredential(accountName, accountKey);
        def  endpoint = String.format(Locale.ROOT, "https://%s.blob.core.windows.net", accountName);
        storageClient = new BlobServiceClientBuilder().endpoint(endpoint).credential(credential).buildClient();
        and:
        Global.session = Mock(Session) { getConfig()>>Map.of() }
    }

    def cleanupSpec() {
        Global.session = null
    }

    def 'should acquire lock on Azure Blob Storage'() {
        given:
        def bucket = createBucket()
        def workDir = Path.of("az://${bucket}/test-locks-work/")
        def manager = new GlobalCacheLockManager(workDir)
        def hash = 'abcdef1234'

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
        if( bucket ) deleteBucket(bucket)
    }

    def 'should fail to acquire lock when already held on Azure'() {
        given:
        def bucket = createBucket()
        def workDir = Path.of("az://${bucket}/test-locks-work/")
        def manager = new GlobalCacheLockManager(workDir)
        def hash = 'abcdef5678'

        // Ensure parent directory exists
        Files.createDirectories(workDir)

        // Clean up any existing lock file
        def lockPath = workDir.resolve(".${hash}.lock")
        Files.deleteIfExists(lockPath)

        when:
        def lock1 = manager.acquire(hash)
        def lock2 = manager.acquire(hash)

        then:
        lock1 != null
        lock2 == null

        cleanup:
        if( bucket ) deleteBucket(bucket)
    }

    def 'should release lock on Azure'() {
        given:
        def bucket = createBucket()
        def workDir = Path.of("az://${bucket}/test-locks-work/")
        def manager = new GlobalCacheLockManager(workDir)
        def hash = 'abcdef9012'

        // Ensure parent directory exists
        Files.createDirectories(workDir)

        // Clean up any existing lock file
        def lockPath = workDir.resolve(".${hash}.lock")
        Files.deleteIfExists(lockPath)

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
        if( bucket ) deleteBucket(bucket)
    }

    def 'should allow acquiring lock after release on Azure'() {
        given:
        def bucket = createBucket()
        def workDir = Path.of("az://${bucket}/test-locks-work/")
        def manager = new GlobalCacheLockManager(workDir)
        def hash = 'abcdef3456'

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
        if( bucket ) deleteBucket(bucket)
    }

    def 'should handle concurrent lock attempts on Azure'() {
        given:
        def bucket = createBucket()
        def workDir = Path.of("az://${bucket}/test-locks-work/")
        def manager = new GlobalCacheLockManager(workDir)
        def hash = 'abcdef7890'
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
        if( bucket ) deleteBucket(bucket)
    }
}
