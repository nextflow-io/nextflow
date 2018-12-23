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
