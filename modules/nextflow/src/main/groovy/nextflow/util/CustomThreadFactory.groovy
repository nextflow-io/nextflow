/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import java.util.concurrent.ThreadFactory
import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic

import java.lang.Thread.UncaughtExceptionHandler

import groovy.util.logging.Slf4j

/**
 * A customised thread factory
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CustomThreadFactory implements ThreadFactory {

    private ThreadGroup group

    private AtomicInteger threadNumber = new AtomicInteger(1)

    private UncaughtExceptionHandler exceptionHandler

    private prefix

    CustomThreadFactory(String prefix, UncaughtExceptionHandler exceptionHandler=null) {
        this.prefix = prefix ?: ThrottlingExecutor.simpleName
        this.group = System.getSecurityManager()?.getThreadGroup() ?: Thread.currentThread().getThreadGroup()
        this.exceptionHandler = exceptionHandler
    }


    Thread newThread(Runnable r) {
        final name = "${prefix}-${threadNumber.getAndIncrement()}"
        log.trace "Creating thread: $name"

        def thread = new Thread(group, r, name, 0)
        if (thread.isDaemon())
            thread.setDaemon(false);
        if (thread.getPriority() != Thread.NORM_PRIORITY)
            thread.setPriority(Thread.NORM_PRIORITY)
        if( exceptionHandler )
            thread.setUncaughtExceptionHandler(exceptionHandler)
        return thread
    }
}