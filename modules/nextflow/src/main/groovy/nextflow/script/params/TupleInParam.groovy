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
import nextflow.script.TokenEnvCall
import nextflow.script.TokenFileCall
import nextflow.script.TokenPathCall
import nextflow.script.TokenStdinCall
import nextflow.script.TokenValCall
import nextflow.script.TokenVar

/**
 * Models a tuple of input parameters
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@InheritConstructors
class TupleInParam extends BaseInParam {

    protected List<BaseInParam> inner = []

    @Override String getTypeName() { 'set' }

    @Override
    TupleInParam clone() {
        final copy = (TupleInParam)super.clone()
        copy.@inner = new ArrayList<>(inner.size())
        for( BaseInParam p : inner ) {
            copy.@inner.add((BaseInParam)p.clone())
        }
        return copy
    }

    String getName() { '__$'+this.toString() }

    TupleInParam bind(Object... obj ) {

        for( def item : obj ) {

            if( item instanceof TokenVar )
                newItem(ValueInParam).bind(item)

            else if( item instanceof TokenFileCall )
                newItem(FileInParam).bind( item.target )

            else if( item instanceof TokenPathCall ) {
                newItem(FileInParam)
                        .setPathQualifier(true)
                        .setOptions(item.opts)
                        .bind( item.target )
            }

            else if( item instanceof Map )
                newItem(FileInParam).bind(item)

            else if( item instanceof TokenValCall )
                newItem(ValueInParam).bind(item.val)

            else if( item instanceof TokenEnvCall )
                newItem(EnvInParam).bind(item.val)

            else if( item instanceof TokenStdinCall )
                newItem(StdInParam)

            else if( item instanceof GString )
                newItem(FileInParam).bind(item)

            else if( item == '-' )
                newItem(StdInParam)

            else if( item instanceof String )
                newItem(FileInParam).bind(item)

            else
                throw new IllegalArgumentException()
        }

        return this

    }

    private <T extends BaseInParam> T newItem( Class<T> type )  {
        type.newInstance(binding, inner, index)
    }

}
