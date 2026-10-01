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

package test

import java.util.concurrent.ConcurrentLinkedQueue
import java.util.concurrent.CountDownLatch
import java.util.concurrent.TimeUnit
import java.util.concurrent.locks.ReentrantLock

import nextflow.util.LockManager
import spock.lang.Specification
import spock.lang.Timeout
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LockManagerTest extends Specification {

    def 'should reuse the same instance from the pool' () {
        given:
        def manager = new LockManager()

        when:
        def lock1 = manager.acquire(1)
        and:
        lock1.release()

        and:
        def copy = manager.acquire(1)
        then:
        copy.is(lock1)

    }

    def 'should lock on the same key' () {
        given:
        def manager = new LockManager()

        when:
        int counter=0
        def lock1 = manager.acquire(1)
        Thread.start { def l = manager.acquire(1); counter++; l.release() }
        and:
        sleep 100
        then:
        counter==0

        when:
        lock1.release()
        sleep 100
        then:
        counter ==1
    }

    def 'should not lock on different keys' () {
        given:
        def manager = new LockManager()

        when:
        int counter=0
        def lock1 = manager.acquire(1)
        Thread.start { def l = manager.acquire(2); counter++; l.release() }
        and:
        sleep 100
        then:
        counter==1

    }

    def 'should return the handle to the pool once fully released' () {
        given:
        def manager = new LockManager()

        when:
        // release a handle for one key, then acquire a *different* key
        // the released handle can only be handed back for the new key
        // if it was actually returned to the pool on release
        def lock1 = manager.acquire('a')
        lock1.release()
        and:
        def lock2 = manager.acquire('b')
        lock2.release()

        then:
        lock2.is(lock1)
    }

    @Timeout(15)
    def 'should not hand out a retired handle to a thread waiting on the same key' () {
        given:
        def manager = new LockManager()
        def errors = new ConcurrentLinkedQueue<Throwable>()
        def releaseAll = new CountDownLatch(1)
        def bAcquired = new CountDownLatch(1)
        def cAcquired = new CountDownLatch(1)
        def worker = { Object key, CountDownLatch acquired ->
            Thread.startDaemon {
                try {
                    def lock = manager.acquire(key)
                    acquired.countDown()
                    releaseAll.await()
                    lock.release()
                }
                catch (Throwable t) { errors.add(t) }
            }
        }

        when: 'thread B waits for the handle of key `k` held by this thread'
        def h1 = manager.acquire('k')
        def b = worker('k', bAcquired)
        waitUntil { queued(h1, b) }
        and: 'the handle is retired while its lock is still held, keeping B parked on it'
        h1.sync.lock()
        h1.release()
        and: 'thread C recycles the retired handle for another key'
        def c = worker('other', cAcquired)
        waitUntil { queued(h1, c) }
        and: 'this thread acquires key `k` again, getting a new handle'
        def h2 = manager.acquire('k')
        and: 'the waiters are woken up'
        h1.sync.unlock()
        waitUntil { bAcquired.count==0 || (cAcquired.count==0 && queued(h2, b)) }
        then: 'B does not enter `k` while this thread holds it, C holds `other`'
        bAcquired.count == 1
        cAcquired.count == 0
        manager.@entries.get('k').is(h2)
        manager.@entries.get('other').is(h1)

        when:
        h2.release()
        then:
        bAcquired.await(10, TimeUnit.SECONDS)

        when:
        releaseAll.countDown()
        b.join(10_000); c.join(10_000)
        then:
        errors.isEmpty()
        !b.isAlive() && !c.isAlive()
        manager.@entries.isEmpty()
        manager.@pool.every { LockManager.LockHandle it -> !((ReentrantLock)it.sync).isLocked() }
    }

    @Timeout(15)
    def 'should not deadlock releasing a key while creating another key in the same map bin' () {
        given:
        def manager = new LockManager()
        def k1 = new GatedKey()
        def k2 = new GatedKey()
        def yHolds = new CountDownLatch(1)
        def yRelease = new CountDownLatch(1)
        def done = new CountDownLatch(2)

        when: 'thread Y releases key k2 and is paused in the map removal (hashing k2)'
        Thread.startDaemon { def lock = manager.acquire(k2); yHolds.countDown(); yRelease.await(); lock.release(); done.countDown() }
        yHolds.await()
        def gate = k2.arm()
        yRelease.countDown()
        k2.entered.await()
        and: 'thread X creates the entry for k1 which lands in the same map bin as k2'
        def x = Thread.startDaemon { def lock = manager.acquire(k1); lock.release(); done.countDown() }
        waitUntil { x.state == Thread.State.BLOCKED || !x.isAlive() }
        and: 'Y is resumed'
        gate.countDown()

        then:
        done.await(10, TimeUnit.SECONDS)
        manager.@entries.isEmpty()
    }

    /**
     * A key with a constant hash code (so that all instances collide in the same map bin)
     * that can pause the next thread computing its hash code
     */
    static class GatedKey {
        volatile CountDownLatch entered
        volatile CountDownLatch gate

        CountDownLatch arm() {
            entered = new CountDownLatch(1)
            gate = new CountDownLatch(1)
            return gate
        }

        @Override
        int hashCode() {
            final g = gate
            if( g != null ) {
                gate = null
                entered.countDown()
                g.await()
            }
            return 42
        }

        @Override
        boolean equals(Object other) { this.is(other) }
    }

    static private boolean queued(LockManager.LockHandle handle, Thread thread) {
        ((ReentrantLock)handle.sync).hasQueuedThread(thread)
    }

    static private void waitUntil(Closure<Boolean> condition) {
        final deadline = System.currentTimeMillis() + 10_000
        while( !condition.call() ) {
            if( System.currentTimeMillis() > deadline )
                throw new AssertionError("Condition not met within timeout")
            Thread.onSpinWait()
        }
    }
}
