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
import java.nio.charset.Charset

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CharsetHelper {

    static Charset getCharset( object ) {

        if( object instanceof Map ) {
            if( object.containsKey('charset') )
                object = (object as Map).charset
            else
                return Charset.defaultCharset()
        }

        if( object instanceof Charset )
            return (Charset)object

        if( object instanceof String && Charset.isSupported(object) )
            return Charset.forName(object)

        if( object != null )
            log.warn "Invalid charset object: $object -- using defualt: ${Charset.defaultCharset()}"

        Charset.defaultCharset()
    }

}
