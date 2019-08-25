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
class FileInParam extends BaseInParam implements PathQualifier {

    protected filePattern

    private boolean pathQualifier

    @Override String getTypeName() { pathQualifier ? 'path' : 'file' }

    @Override String getTypeSimpleName() { getTypeName() + "inparam" }

    /**
     * Define the file name
     */
    FileInParam name( obj ) {
        if(pathQualifier)
            throw new MissingMethodException("name", this.class, [String] as Object[])

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
            assert !pathQualifier
            def entry = bindObject.entrySet().first()
            return entry?.key
        }

        if( bindObject instanceof GString ) {
            return '__$' + this.toString()
        }

        return super.getName()
    }

    @Override
    BaseInParam bind( obj ) {
        if( pathQualifier && obj instanceof Map )
            throw new IllegalArgumentException("Input `path` does not allow such argument: ${obj.entrySet().collect{"${it.key}:${it.value}"}.join(',')}")
        super.bind(obj)
        return this
    }

    String getFilePattern(Map ctx = null) {

        if( filePattern != null  )
            return resolve(ctx,filePattern)

        if( bindObject instanceof Map ) {
            assert !pathQualifier
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

    @Override
    FileInParam setPathQualifier(boolean flag) {
        pathQualifier = flag
        return this
    }

    @Override
    boolean isPathQualifier() { pathQualifier }

    @Override
    FileInParam setOptions(Map<String,?> opts) {
        (FileInParam)super.setOptions(opts)
    }

    /**
     * Defines the `stageAs:` option to define the input file stage name pattern
     *
     * @param value
     *      A string representing the target file name or a file name pattern
     *      ie. containing the star `*` or question mark wildcards
     * @return
     *      The param instance itself
     */
    FileInParam setStageAs(String value) {
        this.filePattern = value
        return this
    }

    FileInParam setName(String value) {
        this.filePattern = value
        return this
    }

}
