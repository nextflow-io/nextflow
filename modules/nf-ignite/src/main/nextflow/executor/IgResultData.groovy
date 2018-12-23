/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.executor

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import nextflow.util.KryoHelper

/**
 * Models the result of a remote closure task execution
 */
@CompileStatic
@EqualsAndHashCode
class IgResultData implements Serializable {

    private static final long serialVersionUID = - 7200781198107958188L ;

    /**
     * The closure returned value serialized as a byte array
     */
    private byte[] fValueObj

    /**
     * The closure execution context serialized as a byte array
     */
    private byte[] fContextObj

    transient Object value

    transient Map context

    Throwable error

    def getValue() {
        if( !value && fValueObj != null ) {
            value = KryoHelper.deserialize(fValueObj)
        }
        return value
    }

    void setValue( obj ) {
        this.value = obj
        this.fValueObj = KryoHelper.serialize(obj)
    }

    Map getContext() {
        if( context == null && fContextObj != null ) {
            context = (Map)KryoHelper.deserialize(fContextObj)
        }
        return context
    }

    void setContext( Map ctx ) {
        this.context = ctx
        this.fContextObj = KryoHelper.serialize(ctx)
    }

}
