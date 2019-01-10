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

package nextflow.util

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 *  A object deserializer that allows you to specify a custom class loader
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
class InputStreamDeserializer extends ObjectInputStream {

    ClassLoader loader

    InputStreamDeserializer( InputStream input, ClassLoader loader ) {
        super(input)
        this.loader = loader
    }

    protected Class<?> resolveClass(ObjectStreamClass desc)
            throws IOException, ClassNotFoundException
    {
        try {
            return super.resolveClass(desc)
        }
        catch( ClassNotFoundException e ) {
            log.debug "Trying to load class: ${desc.getName()}"
            return Class.forName( desc.getName(), false, loader )
        }
    }

    /**
     * Deserialize an object given its bytes representation and a custom class-loader
     *
     * @param bytes
     * @param loader
     * @return
     */
    static <T> T deserialize( byte[] bytes, ClassLoader loader ) {
        assert bytes

        def buffer = new ByteArrayInputStream(bytes)
        (T)new InputStreamDeserializer(buffer, loader).readObject()
    }

}
