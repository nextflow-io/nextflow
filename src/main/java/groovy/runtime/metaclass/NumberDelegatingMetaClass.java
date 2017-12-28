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

package groovy.runtime.metaclass;

import groovy.lang.MetaClass;
import groovy.lang.MissingMethodException;
import groovy.lang.MissingPropertyException;
import nextflow.util.Duration;
import nextflow.util.MemoryUnit;

/**
 * Adds support for custom duration and memory unit to number literals e.g.
 * 2.MB, 5.GB or 1.hour, 2.days, etc
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class NumberDelegatingMetaClass extends groovy.lang.DelegatingMetaClass {

    public NumberDelegatingMetaClass(MetaClass delegate) {
        super(delegate);
    }

    public Object getProperty( Object obj, String property ) {

        try {
            return super.getProperty(obj,property);
        }
        catch( MissingPropertyException e ) {
            String unit;
            if( Duration.UNITS.contains((unit=property.toLowerCase())) ) {
                return new Duration(obj.toString() + unit);
            }
            else if(MemoryUnit.UNITS.contains((unit=property.toUpperCase())) ) {
                return new MemoryUnit(obj.toString() + unit);
            }
            // fallback on the original exception
            throw e;
        }
    }

    @Override
    public Object invokeMethod(Object obj, String methodName, Object[] args) {

        try {
            return delegate.invokeMethod(obj, methodName, args);
        }
        catch( MissingMethodException e ) {
            int len = args.length;
            if( len != 1 ) throw e;
            Object operand = args[0];
            Number number = (Number)obj;

            /*
             * extend 'multiply' for custom objects
             */
            if( "multiply".equals(methodName) ) {
                if( operand instanceof Duration ) {
                    return ((Duration)operand).multiply(number);
                }
                else if( operand instanceof MemoryUnit ) {
                    return ((MemoryUnit)operand).multiply(number);
                }
            }
            // default fall back on MissingMethodException
            throw e;
        }

    }

}
