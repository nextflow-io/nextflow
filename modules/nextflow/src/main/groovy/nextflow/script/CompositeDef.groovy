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

package nextflow.script

import groovy.transform.CompileStatic

/**
 * Models a composition of process and workflow components
 * in a pipe operation expression
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CompositeDef extends ComponentDef implements ChainableDef {

    private List<ChainableDef> elements = new ArrayList<>(5)

    CompositeDef add(ChainableDef comp) {
        elements.add(comp)
        return this
    }

    String getType() { 'composite' }

    CompositeDef cloneWithName(String name) {
        throw new UnsupportedOperationException()
    }

    @Override
    String getName() {
        return "( ${elements.collect{ it.name }.join(' & ')} )"
    }

    @Override
    Object invoke_a(Object[] args) {
        int i=0
        def result = new ArrayList(elements.size())
        for( def entry : elements )
            result[i++] = entry.invoke_a(args)

        new ChannelOut(ChannelOut.spread(result))
    }

    @Override
    String getSignature() {
        "expression $name"
    }

}
