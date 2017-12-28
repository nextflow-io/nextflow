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
import java.io.InputStream;
import java.io.Reader;
import java.nio.file.Path;

import groovy.lang.MetaClass;
import groovy.lang.MetaClassRegistry;
import groovyx.gpars.dataflow.DataflowQueue;
import groovyx.gpars.dataflow.DataflowVariable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Register a generic meta-class handler to provide extra "dynamic" extension methods.
 * <p>
 * It add the splitter methods to the following classes:
 * <li>{@link Path}
 * <li>{@link File}
 * <li>{@link String}
 * <li>{@link java.io.InputStream}
 * <li>{@link Reader}
 * <li>{@link groovyx.gpars.dataflow.DataflowVariable}
 * <li>{@link groovyx.gpars.dataflow.DataflowQueue}
 *
 *
 *
 * @link http://www.objectpartners.com/2013/07/30/customizing-mop-in-groovy/
 * @link http://stackoverflow.com/questions/543479/groovy-delegating-metaclass-for-an-interface
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class CustomMetaClassCreationHandle extends MetaClassRegistry.MetaClassCreationHandle {

    static final Logger log = LoggerFactory.getLogger(CustomMetaClassCreationHandle.class);

    protected MetaClass createNormalMetaClass(Class theClass, MetaClassRegistry registry) {
        MetaClass metaClass = super.createNormalMetaClass( theClass, registry );

        if( Number.class.isAssignableFrom(theClass) ) {
            log.trace("Registering number meta-class for: {}", theClass);
             return new NumberDelegatingMetaClass(metaClass);
        }
        else if (isSplitterClass(theClass)) {
            log.trace("Registering custom meta-class for: {}", theClass);
            return new NextflowDelegatingMetaClass(metaClass);
        }

        return metaClass;
    }

    protected boolean isSplitterClass( Class theClass ) {
        return  String.class == theClass ||
                File.class == theClass ||
                Path.class.isAssignableFrom(theClass) ||
                InputStream.class.isAssignableFrom(theClass) ||
                Reader.class.isAssignableFrom(theClass) ||
                DataflowVariable.class.isAssignableFrom(theClass) ||
                DataflowQueue.class.isAssignableFrom(theClass);
    }

}
