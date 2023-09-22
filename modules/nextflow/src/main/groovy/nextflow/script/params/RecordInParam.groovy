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

package nextflow.script.params

import groovy.transform.InheritConstructors
import nextflow.NF
import nextflow.script.TokenEnvCall
import nextflow.script.TokenFileCall
import nextflow.script.TokenPathCall
import nextflow.script.TokenStdinCall
import nextflow.script.TokenValCall
import nextflow.script.TokenVar

/**
 * Models a custom record input
 * 
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@InheritConstructors
class RecordInParam extends BaseInParam {

    protected List<String> innerFields = []

    protected List<BaseInParam> inner = []

    @Override
    String getTypeName() { 'record' }

    List<BaseInParam> getInner() { inner }

    @Override
    RecordInParam clone() {
        final copy = (RecordInParam)super.clone()
        copy.@innerFields = new ArrayList<>(innerFields.size())
        copy.@inner = new ArrayList<>(inner.size())
        for( BaseInParam p : inner ) {
            copy.@inner.add((BaseInParam)p.clone())
        }
        return copy
    }

    String getName() { '__$'+this.toString() }

    RecordInParam bind(Map opts, Class clazz) {
        for( def entry : opts ) {
            def item = entry.value
            def nested

            if( item instanceof TokenValCall ) {
                nested = create(ValueInParam).bind(item.val)
            }
            else if( item instanceof TokenEnvCall ) {
                nested = create(EnvInParam).bind(item.val)
            }
            else if( item instanceof TokenStdinCall || item == '-' ) {
                nested = create(StdInParam)
            }
            else if( item instanceof TokenFileCall ) {
                nested = create(FileInParam).bind( item.target )
            }
            else if( item instanceof TokenPathCall ) {
                nested = create(FileInParam)
                        .setPathQualifier(true)
                        .setOptions(item.opts)
                        .bind( item.target )
            }
            else
                throw new IllegalArgumentException("Invalid `record` input parameter declaration -- item: ${item}")

            innerFields.add(entry.key)
            nested.owner = this
        }

        return this

    }

    private <T extends BaseInParam> T create( Class<T> type )  {
        type.newInstance(binding, inner, index)
    }

    List extractValues(Map map) {
        innerFields.collect { map[it] }
    }

}
