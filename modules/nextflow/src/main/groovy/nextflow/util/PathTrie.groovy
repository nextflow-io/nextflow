/*
 * Copyright 2013-2024, Seqera Labs
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

import java.nio.file.FileSystems
import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.file.FileSystemPathFactory
/**
 * A trie data structure specialised to find the longest common paths
 * in a given list of paths
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class PathTrie {

    private static String PATH_SEP = File.separator

    // and empty string used to mark the end of a path
    private static final String END_PATH = ''

    List<Trie<String>> paths = []

    PathTrie() {}

    /**
     * Get or create a not for the given path item.
     *
     * @param item A path item represented as string
     * @return A {@code Trie<String>} for the given path item
     */
    protected Trie<String> getOrCreate( String item ) {
        def found = paths.find { Trie it -> it.node == item }
        if( !found ) {
            found = new Trie<String>(item)
            paths << found
        }
        return found
    }

    /**
     * Add a path to the trie collection
     *
     * @param path The path to add, it can an absolute or relative path
     */
    void add( String path )  {
        assert path
        if( path=='/' ) return
        
        def splitter = PathSplitter.parse(path)
        def head = splitter.head
        def tail = splitter.tail
        // add an extra entry to mark the end of a path
        if( tail ) tail.add(END_PATH)
        else tail = [END_PATH]
        getOrCreate(head).addPath(tail)
    }

    /**
     * Add a file to the collection of path
     *
     * @param file
     */
    void add( File file ) {
        add(file.toString())
    }

    /**
     * Add a string path to the collection
     *
     * @param path
     */
    void add( Path path ) {
        final str = path.fileSystem == FileSystems.default
                ? path.toString()
                : FileSystemPathFactory.getUriString(path)
        add( str )
    }

    /**
     * Retrieve the list of the longest paths. For example, given:
     *
     * <pre>
     *  /home/data/work
     *  /home/data/work/xx/file_x
     *  /db/data/tutorial
     *  /db/data/xxx
     * </pre>
     *
     * Il returns a list containing two paths:
     *
     * <pre>
     *     /home/data/work
     *     /db/data
     * </pre>
     *
     * @return
     */
    List<String> longest() {

        List<String> result = new LinkedList<String>()
        for( Trie<String> it : paths ) {

            List<String> tokens = it.longest()
            if( !tokens )
                break

            // remove the entry marking the end of the path
            def last = tokens.size()-1
            if( tokens[last] == END_PATH )
                tokens.remove(last)
            result.add(tokens.join(PATH_SEP))
        }

        return result
    }

    List<String> traverse() {
        List<String> result = new LinkedList<String>()
        for( Trie<String> it : paths ) {
            def allTrails = it.traverse(END_PATH)
            for( List<String> t : allTrails ) {
                result.add( t.join('/') )
            }
        }
        return result
    }
}

