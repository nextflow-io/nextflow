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

package nextflow.file.http

import com.esotericsoftware.kryo.Kryo
import com.esotericsoftware.kryo.Serializer
import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.file.FileHelper

/**
 * Implements Kryo serializer for {@link XPath}
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class XPathSerializer extends Serializer<XPath> {

    @Override
    void write(Kryo kryo, Output output, XPath target) {
        final uri = target.toUri().toString()
        log.trace "XPath serialization > uri=$uri"
        output.writeString(uri)
    }

    @Override
    XPath read(Kryo kryo, Input input, Class<XPath> type) {
        final uri = input.readString()
        log.trace "Path de-serialization > uri=$uri"
        (XPath) FileHelper.asPath(new URI(uri))
    }
}
