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

package nextflow.splitter

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
/**
 * Base class for splitter handling binary data
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@CompileStatic
@InheritConstructors
abstract class AbstractBinarySplitter extends AbstractSplitter<InputStream> {

    protected InputStream normalizeSource( obj ) {

        if( obj instanceof InputStream )
            return (InputStream) obj

        if( obj instanceof byte[] )
            return new ByteArrayInputStream((byte[])obj)

        if( obj instanceof CharSequence )
            return new ByteArrayInputStream(obj.toString().bytes)

        if( obj instanceof Path )
            return newInputStream(obj)

        if( obj instanceof File )
            newInputStream(obj.toPath())

        throw new IllegalAccessException("Object of class '${obj.class.name}' does not support 'splitter' methods")

    }



}
