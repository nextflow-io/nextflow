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
