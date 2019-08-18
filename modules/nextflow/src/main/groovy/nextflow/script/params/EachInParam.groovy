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
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.expression.DataflowExpression
import nextflow.extension.CH
import nextflow.extension.ToListOp
import nextflow.script.TokenFileCall
import nextflow.script.TokenPathCall

/**
 * Represents a process input *iterator* parameter
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@InheritConstructors
@Slf4j
class EachInParam extends BaseInParam {

    @Override String getTypeName() { 'each' }

    private List<BaseInParam> inner = []

    String getName() { '__$'+this.toString() }

    Object clone() {
        final copy = (EachInParam)super.clone()
        copy.@inner = new ArrayList<>(inner.size())
        for( BaseInParam p : inner ) {
            copy.@inner.add((BaseInParam)p.clone())
        }
        return copy
    }

    EachInParam bind( def obj ) {
        final nested = createNestedParam(obj)
        nested.owner = this
        this.bindObject = nested.bindObject
        return this
    }

    protected BaseInParam createNestedParam(obj) {
        if( obj instanceof TokenFileCall ) {
            return new FileInParam(binding, inner, index)
                    .bind(obj.target)
        }
        
        if( obj instanceof TokenPathCall ) {
            return new FileInParam(binding, inner, index)
                    .setPathQualifier(true)
                    .bind(obj.target)
        }

        return new ValueInParam(binding, inner, index)
                .bind(obj)
    }

    InParam getInner() { inner[0] }

    @Override
    protected DataflowReadChannel inputValToChannel( value ) {
        def variable = normalizeToVariable( value )
        super.inputValToChannel(variable)
    }

    @PackageScope
    DataflowReadChannel normalizeToVariable( value ) {
        def result
        if( value instanceof DataflowExpression ) {
            result = value
        }
        else if( CH.isChannel(value) ) {
            def read = CH.getReadChannel(value)
            result = new ToListOp(read).apply()
        }
        else {
            result = new DataflowVariable()
            result.bind(value)
        }

        return result.chainWith { it instanceof Collection || it == null ? it : [it] }
    }

}
