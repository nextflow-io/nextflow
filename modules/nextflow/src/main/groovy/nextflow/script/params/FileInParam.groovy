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

package nextflow.script.params

import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import nextflow.script.TokenVar

/**
 * Represents a process *file* input parameter
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@InheritConstructors
class FileInParam extends BaseInParam  {

    protected filePattern

    @Override String getTypeName() { 'file' }

    /**
     * Define the file name
     */
    FileInParam name( obj ) {
        if( obj instanceof String ) {
            filePattern = obj
            return this
        }

        if( obj instanceof GString ) {
            filePattern = obj
            return this
        }

        // the ability to pass a closure as file name has been replaced by
        // lazy gstring -- this should be deprecated
        if( obj instanceof Closure ) {
            filePattern = obj
            return this
        }

        throw new IllegalArgumentException()
    }

    String getName() {

        if( bindObject instanceof Map ) {
            def entry = bindObject.entrySet().first()
            return entry?.key
        }

        if( bindObject instanceof GString ) {
            return '__$' + this.toString()
        }

        return super.getName()

    }

    String getFilePattern(Map ctx = null) {

        if( filePattern != null  )
            return resolve(ctx,filePattern)

        if( bindObject instanceof Map ) {
            def entry = bindObject.entrySet().first()
            return resolve(ctx, entry?.value)
        }

        if( bindObject instanceof TokenVar )
            return filePattern = '*'

        if( bindObject != null )
            return resolve(ctx, bindObject)

        return filePattern = '*'
    }

    private resolve( Map ctx, value ) {
        if( value instanceof GString ) {
            value.cloneAsLazy(ctx)
        }

        else if( value instanceof Closure ) {
            return ctx.with(value)
        }

        else
            return value
    }

}
