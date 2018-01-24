/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
