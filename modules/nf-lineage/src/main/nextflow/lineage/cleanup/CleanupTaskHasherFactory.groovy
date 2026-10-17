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

import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import nextflow.processor.TaskHasher
import nextflow.processor.TaskHasherFactory
import nextflow.processor.TaskRun
import nextflow.script.params.InParam
import nextflow.util.CacheHelper
/**
 * Hash tasks whose input files were deleted by eager cleanup,
 * using the lineage metadata of the deleted files, so that the task
 * hash is the same as in the previous run.
 *
 * Tasks that hash input files by content (i.e. `deep` and `sha256`
 * cache modes) are not supported, since the deleted files can't be
 * hashed without restoring them.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class CleanupTaskHasherFactory implements TaskHasherFactory {

    @Override
    TaskHasher create(TaskRun task) {
        final deleted = DeletedOutputs.current()
        if( deleted == null || !deleted.hasDeletedInputs(task) )
            return null
        final mode = task.processor.config.getHashMode()
        if( mode == CacheHelper.HashMode.DEEP || mode == CacheHelper.HashMode.SHA256 )
            return null
        return new Hasher(task, deleted)
    }

    static class Hasher extends TaskHasher {

        private TaskRun task

        private DeletedOutputs deleted

        Hasher(TaskRun task, DeletedOutputs deleted) {
            super(task)
            this.task = task
            this.deleted = deleted
        }

        @Override
        HashCode compute() {
            final inputs = new LinkedHashMap<InParam,Object>(task.inputs)
            task.inputs.replaceAll( (k, v) -> deleted.substitute(v) )
            try {
                return super.compute()
            }
            finally {
                task.inputs.clear()
                task.inputs.putAll(inputs)
            }
        }
    }

}
