/*
 * Copyright 2013-2023, Seqera Labs
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

import groovy.transform.CompileStatic

/**
 * Models a file input directive, which defines the file
 * or set of files to be staged into the task environment.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class TaskFileInput implements PathArityAware {
    private Object value
    private boolean coerceToPath
    private String name
    private Object filePattern

    TaskFileInput(Object value, boolean coerceToPath, String name, Map<String,?> opts) {
        this.value = value
        this.coerceToPath = coerceToPath
        this.name = name
        this.filePattern = opts.stageAs ?: opts.name

        if( opts.arity )
            this.setArity(opts.arity.toString())
    }

    Object getValue(Map ctx) {
        return resolve(ctx, value)
    }

    boolean isPathQualifier() {
        return coerceToPath
    }

    String getName() {
        return name
    }

    String getFilePattern(Map ctx) {
        if( filePattern != null )
            return resolve(ctx, filePattern)

        if( value != null )
            return resolve(ctx, value)

        return filePattern = '*'
    }

    private Object resolve( Map ctx, value ) {
        if( value instanceof GString )
            return value.cloneAsLazy(ctx)

        if( value instanceof Closure )
            return ctx.with(value)

        return value
    }

}
