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

package groovy.runtime.metaclass;

import java.io.File;
import java.io.Reader;
import java.nio.file.Path;

import groovy.lang.MetaClass;
import groovy.lang.MetaClassRegistry;
import groovyx.gpars.dataflow.DataflowBroadcast;
import groovyx.gpars.dataflow.DataflowReadChannel;
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


    static final private Logger log = LoggerFactory.getLogger(CustomMetaClassCreationHandle.class);

    protected MetaClass createNormalMetaClass(Class theClass, MetaClassRegistry registry) {
        MetaClass metaClass = super.createNormalMetaClass( theClass, registry );

        if( Number.class.isAssignableFrom(theClass) ) {
            log.trace("Registering number meta-class for: {}", theClass);
             return new NumberDelegatingMetaClass(metaClass);
        }
        else if (isExtensionClass(theClass)) {
            log.trace("Registering custom meta-class for: {}", theClass);
            return new NextflowDelegatingMetaClass(metaClass);
        }

        return metaClass;
    }

    protected boolean isExtensionClass(Class theClass) {
        return  File.class == theClass ||
                Path.class.isAssignableFrom(theClass) ||
                DataflowBroadcast.class.isAssignableFrom(theClass) ||
                DataflowReadChannel.class.isAssignableFrom(theClass) ||
                // NOTE: groovy class cannot be referenced explicitly
                // otherwise it creates a circular dependencies causing
                // the compile to crash. Also the class name needs to be
                // resolved statically because this result is cached.
                "nextflow.script.ChannelOut".equals(theClass.getName())
                ;
    }

}
