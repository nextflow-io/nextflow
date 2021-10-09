/*
 * Copyright 2020-2021, Seqera Labs
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
 *
 */

package nextflow.executor

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import nextflow.util.Escape

/**
 * Task output name helper class
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ScriptOutputFiles implements Map<String,Glob> {

    @Delegate
    private LinkedHashMap<String,Glob> target = new LinkedHashMap<>(20)

    /**
     * Model a boolean value representing a glob pattern flag
     *
     * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
     */
    @Canonical
    static class Glob {
        boolean value
        boolean asBoolean() { value }
    }

    static Glob glob(boolean value) { new Glob(value) }

    static ScriptOutputFiles wrap(String... names) {
        final result = new ScriptOutputFiles()
        for( String it : names )
            result.put(it,new Glob(true))
        return result
    }

    /**
     * Escape special characters from file names prepending with a backslash
     *
     * @param name
     *      The file name string to be escaped
     * @param glob
     *      When {@code true} the name is considered a glob pattern, therefore chars like {@code *} and {@core ?}
     *      are not escaped, when {@code false} also these chars are escaped. Note, special chars like blank and round
     *      parenthesis are escaped in any case irrespective of this flag.
     * @return
     *      The name containing escaping special chars to make it shell friendly
     */
    static String shellSpecialChars(String name, Glob glob) {
        return glob ? Escape.path(removeGlobStar(name)) : Escape.wildcards(name)
    }

    /**
     * Trim-trailing-slash
     * @return
     */
    static private String tts(String loc) {
        while( loc && loc.size()>1 && loc.endsWith('/') )
            loc = loc.substring(0,loc.length()-1)
        return loc
    }

    /**
     * Remove a glob star specifier returning the longest parent path
     *
     * <pre>
     * /some/data/&#42;&#42;/file.txt  ->   /some/data
     * /some/data/file.txt     ->   /some/data/file.txt
     * /some&#42;&#42;                 ->   *
     * </pre>
     *
     * @param path
     * @return
     */
    static String removeGlobStar(String path) {

        def p = path?.indexOf('**')
        if( p == -1 )
            return path

        def slash = -1
        for( int i=p-1; i>=0; i-- ) {
            if( path.charAt(i) == '/' as char) {
                slash = i
                break
            }
        }

        if( slash == -1 )
            return '*'

        path.substring(0,slash)
    }

    void putName(String name, boolean glob) {
        target.put(name, new Glob(glob))
    }

    String toShellEscapedNames() {
        final result = new ArrayList(target.size())
        for( String it : target.keySet() ) {
            result.add( tts(shellSpecialChars(it, target[it])) )
        }
        return result.join(' ')
    }
}
