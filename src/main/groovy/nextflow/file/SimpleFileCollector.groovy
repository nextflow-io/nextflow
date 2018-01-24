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

package nextflow.file
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ConcurrentMap

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
/**
 *  Helper class used to aggregate values having the same key
 *  to files
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SimpleFileCollector extends FileCollector {

    protected ConcurrentMap<String,Path> cache = new ConcurrentHashMap<>()

    SimpleFileCollector( ) {

    }

    @Override
    SimpleFileCollector add( String name, value ) {

        // -- 1. process the data value, this fetch the header line(s)
        def data = normalizeToStream(value)

        // -- 2. create the file if missing
        Path target = cache.get(name)
        boolean isNew = !target
        if( isNew ) {
            target = Files.createFile(getTempDir().resolve(name))
            cache.put(name, target)
        }

        // -- 3. add the header if required
        def output = Files.newOutputStream(target, APPEND)
        if( isNew )
            appendHeader(data, name, output)

        // -- 4. finally add the data
        appendStream(data, output)
        output.close()
        return this
    }

    /**
     *
     * @return The number of files in the appender accumulator
     */
    int size() {
        cache.size()
    }

    boolean isEmpty() {
        cache.isEmpty()
    }

    boolean containsKey(String key) {
        return cache.containsKey(key)
    }

    Path get(String name) {
        cache.get(name)
    }

    List<Path> getFiles() {
        new ArrayList<Path>(cache.values())
    }

    /**
     * {@inheritDoc}
     */
    @Override
    void saveFile( Closure<Path> closure ) {

        def result = []
        Iterator<Path> itr = cache.values().iterator()
        while( itr.hasNext() ) {
            def item = itr.next()
            def target = closure.call(item.getName())
            result << Files.move(item, target, StandardCopyOption.REPLACE_EXISTING)
            itr.remove()
        }

    }


}