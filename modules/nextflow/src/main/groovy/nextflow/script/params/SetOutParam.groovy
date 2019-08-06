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
import nextflow.script.TokenFileCall
import nextflow.script.TokenPathCall
import nextflow.script.TokenStdoutCall
import nextflow.script.TokenValCall
import nextflow.script.TokenVar

/**
 * Model a set of process output parameters
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@InheritConstructors
class SetOutParam extends BaseOutParam implements OptionalParam {

    enum CombineMode implements OutParam.Mode { combine }

    final List<BaseOutParam> inner = []

    String getName() { toString() }

    SetOutParam bind( Object... obj ) {

        for( def item : obj ) {
            if( item instanceof TokenVar )
                create(ValueOutParam).bind(item)

            else if( item instanceof TokenValCall )
                create(ValueOutParam).bind(item.val)

            else if( item instanceof GString )
                create(FileOutParam).bind(item)

            else if( item instanceof TokenStdoutCall || item == '-'  )
                create(StdOutParam).bind('-')

            else if( item instanceof String )
                create(FileOutParam).bind(item)

            else if( item instanceof TokenFileCall )
                // note that 'filePattern' can be a string or a GString
                create(FileOutParam).bind(item.target)

            else if( item instanceof TokenPathCall ) {
                // note that 'filePattern' can be a string or a GString
                create(FileOutParam)
                        .setPathQualifier(true)
                        .setOptions(item.opts)
                        .bind(item.target)
            }


            else
                throw new IllegalArgumentException("Invalid `set` output parameter declaration -- item: ${item}")
        }

        return this
    }

    protected <T extends BaseOutParam> T create(Class<T> type) {
        type.newInstance(binding,inner,index)
    }

    @Override
    void lazyInit() {
        super.lazyInit()
        inner.each { opt ->
            if( opt instanceof FileOutParam ) opt.optional(this.optional)
        }
    }

    SetOutParam mode( def value ) {

        def str = value instanceof String ? value : ( value instanceof TokenVar ? value.name : null )
        if( str ) {
            try {
                this.mode = CombineMode.valueOf(str)
            }
            catch( Exception e ) {
                super.mode(value)
            }
        }

        return this
    }
}
