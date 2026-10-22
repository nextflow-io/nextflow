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

import java.nio.file.Path

import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.processor.DefaultTaskCacheStrategy
import nextflow.processor.TaskCacheStrategy
import nextflow.processor.TaskEntry
import nextflow.processor.TaskResolver
import nextflow.processor.TaskRun
/**
 * Resume tasks whose output files were deleted by eager cleanup.
 *
 * A task is resumed optimistically from the lineage of its outputs.
 * If a downstream task needs to be executed, the deleted files are
 * restored by re-executing the tasks that produced them.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class CleanupCacheStrategy implements TaskCacheStrategy {

    private static final TaskCacheStrategy DEFAULT = new DefaultTaskCacheStrategy()

    @Override
    boolean isEnabled(Session session) {
        return DeletedOutputs.get(session) != null
    }

    @Override
    void resolve(TaskRun task, HashCode hash, boolean tryCache, TaskResolver resolver) {
        DEFAULT.resolve(task, hash, tryCache, new Resolver(resolver, DeletedOutputs.current()))
    }

    static class Resolver implements TaskResolver {

        private TaskResolver delegate

        private DeletedOutputs deleted

        Resolver(TaskResolver delegate, DeletedOutputs deleted) {
            this.delegate = delegate
            this.deleted = deleted
        }

        @Override
        TaskEntry entry(HashCode key) {
            return delegate.entry(key)
        }

        /**
         * Resume a task normally, otherwise resume it with the output
         * files that were deleted by eager cleanup.
         */
        @Override
        boolean resume(TaskRun task, HashCode key, Path workDir, TaskEntry entry) {
            if( delegate.resume(task, key, workDir, entry) )
                return true

            // register the outputs before resuming, since downstream tasks
            // are hashed as soon as the outputs are emitted
            final producer = new DeletedOutputs.Producer(task, key, workDir, this)
            final outputs = deleted.register(producer)
            if( !outputs )
                return false
            if( !deleted.isPublishedUpToDate(producer) ) {
                deleted.unregister(outputs)
                return false
            }

            task.deletedOutputs = outputs
            try {
                if( delegate.resume(task, key, workDir, entry) )
                    return true
            }
            finally {
                task.deletedOutputs = null
            }

            deleted.unregister(outputs)
            return false
        }

        /**
         * Restore any deleted input files before launching a task.
         */
        @Override
        void launch(TaskRun task, HashCode key, Path workDir) {
            deleted.restoreInputs(task)
            delegate.launch(task, key, workDir)
        }

        @Override
        Path workDirFor(HashCode key) {
            return delegate.workDirFor(key)
        }
    }

}
