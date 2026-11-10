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
package nextflow.processor

import com.google.common.hash.Hasher
import groovy.transform.CompileStatic
import nextflow.file.FileHolder
import nextflow.util.CacheFunnel
import nextflow.util.CacheHelper

/**
 * Models the files staged by a `stageAs` directive declared in the
 * `stage:` block of a typed process.
 */
@CompileStatic
class TaskStageAsFiles implements CacheFunnel {

    /**
     * The resolved file pattern
     */
    final String filePattern

    /**
     * The files staged by the directive
     */
    final Collection<FileHolder> files

    /**
     * The name of the task input whose value is staged by the directive,
     * or {@code null} if the directive stages any other value
     */
    final String inputName

    TaskStageAsFiles(String filePattern, Collection<FileHolder> files, String inputName) {
        this.filePattern = filePattern
        this.files = files
        this.inputName = inputName
    }

    /**
     * The files of a task input are already hashed with the task inputs,
     * so they are identified by the input name to avoid hashing them twice.
     * Any other files are hashed with their stage names, since a wildcard
     * pattern names them in the iteration order of the staged value.
     */
    @Override
    Hasher funnel(Hasher hasher, CacheHelper.HashMode mode) {
        hasher.putUnencodedChars(filePattern)
        hasher.putChar((char)0)
        if( inputName != null )
            return hasher.putUnencodedChars(inputName).putChar((char)0)
        final named = new HashMap<String,FileHolder>(files.size())
        for( final holder : files )
            named.put(holder.stageName, holder)
        return CacheHelper.hasher(hasher, named, mode)
    }

    @Override
    String toString() {
        return "TaskStageAsFiles[filePattern=${filePattern}; inputName=${inputName}; files=${files}]"
    }

}
