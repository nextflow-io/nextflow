/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
