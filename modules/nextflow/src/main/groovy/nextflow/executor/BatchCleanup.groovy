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

package nextflow.executor

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.Session

/**
 * Delete cluster jobs in a batch manner
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class BatchCleanup {

    /**
     * Default kill batch size
     */
    int size = 100

    /**
     * The current session object
     */
    Session session

    /*
     * helper class
     */
    private static class ExecutorJobsPair {
        AbstractGridExecutor executor
        List jobIds = []

        ExecutorJobsPair(AbstractGridExecutor executor) {
            this.executor = executor
        }

        String toString() { "executor: ${executor.name}; jobs to kill: ${jobIds.join(',')}" }
    }

    @PackageScope
    Map<String, ExecutorJobsPair> aggregate = new HashMap<>()

    /**
     * Add a job ID to the collection of jobs to delete for the specified executor
     *
     * @param executor
     * @param jobId
     */
    @PackageScope
    void collect( AbstractGridExecutor executor, jobId ) {
        final name = executor.name
        def pair = aggregate.get(name)
        if( !pair ) {
            aggregate.put(name, pair = new ExecutorJobsPair(executor))
            if( !session ) {
                session = executor.session
            }
        }
        pair.jobIds.add(jobId)
    }

    /**
     * Group together jobs for the same executor and kill them in bunch of a given
     * size
     */
    void kill() {
        aggregate.values().each { pair ->

            int batchSize =  session ? session.getExecConfigProp(pair.executor.name, 'killBatchSize', size) as int: size
            kill0(pair.executor, pair.jobIds, batchSize)

        }
    }

    /*
     * Implements the kill operation
     */
    @PackageScope
    void kill0( AbstractGridExecutor executor, List ids, int size ){

        int p = 0
        while( p < ids.size() )  {
            int q = Math.min(p+size, ids.size())
            def chunk = ids[p..q-1]
            executor.killTask(chunk)
            p = q
        }

    }

    /**
     * @return A string listing the IDs of the jobs to be killed
     */
    @Override
    String toString() {
        def result = new StringBuilder()
        result += "${this.class.simpleName}[\n"
        aggregate.values().each { result+="${it}\n" }
        result += "]"
        return result
    }

}
