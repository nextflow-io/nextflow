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
