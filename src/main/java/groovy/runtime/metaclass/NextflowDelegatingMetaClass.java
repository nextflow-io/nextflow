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

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import groovy.lang.Closure;
import groovy.lang.MetaClass;
import groovy.lang.MissingMethodException;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowWriteChannel;
import nextflow.extension.DataflowExtensions;
import nextflow.file.FileHelper;
import nextflow.splitter.SplitterFactory;

/**
 * Provides the "dynamic" splitter methods and {@code isEmpty} method for {@link File} and {@link Path} classes.
 *
 * See http://groovy.codehaus.org/Using+the+Delegating+Meta+Class
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class NextflowDelegatingMetaClass extends groovy.lang.DelegatingMetaClass {

    public NextflowDelegatingMetaClass(MetaClass delegate) {
        super(delegate);
    }

    @Override
    public Object invokeMethod(Object obj, String methodName, Object[] args)
    {
        int len = args != null ? args.length : 0;

        /*
         * Implements the 'isEmpty' method for {@link Path} and {@link File} objects.
         * It cannot implemented as a extension method since Path has a private 'isEmpty' method
         *
         * See http://jira.codehaus.org/browse/GROOVY-6363
         */
        if( "isEmpty".equals(methodName) && len==0 ) {
            if( obj instanceof File)
                return FileHelper.empty((File)obj);
            if( obj instanceof Path )
                return FileHelper.empty((Path)obj);
        }
        else if( len>0 && isInto(methodName, obj, args) ) {
            DataflowWriteChannel[] target = new DataflowWriteChannel[len];
            System.arraycopy(args,0,target,0,len);
            DataflowExtensions.into((DataflowReadChannel) obj, target);
            return null;
        }

        /*
         * invoke the method
         */
        try {
            return delegate.invokeMethod(obj, methodName, args);
        }

        /*
         * try to fallback on a splitter method if missing
         */
        catch( MissingMethodException e1 ) {

            // make it possible to invoke some dataflow operator specifying an open array
            // in place of a list object. For example:
            // DataflowReadChannel# separate(final List<DataflowWriteChannel<?>> outputs, final Closure<List<Object>> code)
            // can be invoked as:
            // queue.separate( x, y, z ) { ... }
            if( checkOpenArrayDataflowMethod(NAMES, obj, methodName, args) ) {
                Object[] newArgs = new Object[] { toListOfChannel(args), args[len - 1] };
                return delegate.invokeMethod(obj, methodName, newArgs);
            }

            return SplitterFactory.tryFallbackWithSplitter(obj, methodName, args, e1);
        }
    }

    static private final List<String> NAMES = Arrays.asList("choice","merge","separate");

    /**
     *  Check the following conditions:
     *  <li>method name is the requested one
     *  <li>the target object is a {@link groovyx.gpars.dataflow.DataflowReadChannel}
     *  <li>the all items in the last arguments array are {@link groovyx.gpars.dataflow.DataflowWriteChannel}
     *  <li>the last item in the arguments array is a {@link Closure}
     */
    private static boolean checkOpenArrayDataflowMethod(List<String> validNames, Object obj, String methodName, Object[] args) {
        if( !validNames.contains(methodName) ) return false;
        if( !(obj instanceof DataflowReadChannel)) return false;
        if( args == null || args.length<2 ) return false;
        if( !(args[args.length-1] instanceof Closure) ) return false;
        for( int i=0; i<args.length-1; i++ )
            if( !(args[i] instanceof DataflowWriteChannel )) return false;

        return true;
    }


    /**
     * Given an array of arguments copy the first n-1 to a list of {@link groovyx.gpars.dataflow.DataflowWriteChannel}
     */
    private static List<DataflowWriteChannel> toListOfChannel( Object[] args )  {
        List<DataflowWriteChannel> result = new ArrayList<>(args.length-1);
        for( int i=0; i<args.length-1; i++ ) {
            result.add(i, (DataflowWriteChannel)args[i]);
        }
        return result;
    }

    private static boolean isInto( String name, Object source, Object[] args ) {
        if( !"into".equals(name)) return false;

        if( !(source instanceof DataflowReadChannel) ) return false;

        for( int i=0; i<args.length; i++ )
            if( !(args[i] instanceof DataflowWriteChannel )) return false;

        return true;
    }



}
