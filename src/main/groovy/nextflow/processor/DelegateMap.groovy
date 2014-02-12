/*
 * Copyright (c) 2012, the authors.
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

package nextflow.processor

import java.nio.file.FileSystemNotFoundException
import java.nio.file.FileSystems
import java.nio.file.Path
import java.nio.file.spi.FileSystemProvider

import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.script.BaseScript
import nextflow.util.BlankSeparatedList
import org.apache.commons.lang.SerializationException
import org.apache.commons.lang.SerializationUtils
/**
 * Map used to delegate variable resolution to script scope
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class DelegateMap implements Map {

    @Delegate
    final private Map<String,Object> holder

    private Script script

    private String name

    private boolean undef


    DelegateMap( TaskProcessor processor, Map holder = null ) {
        this.holder = holder ?: [:]
        this.script = processor.ownerScript
        this.undef = processor.taskConfig.getUndef()
        this.name = processor.name
    }

    @Deprecated
    DelegateMap(BaseScript script, boolean undef = false) {
        this.script = script
        this.holder = [:]
        this.undef = undef
    }


    @Override
    public Object get(Object property) {

        if( holder.containsKey(property) ) {
            return holder.get(property)
        }
        // TODO verify if it could better using "script.getBinding().getVariable()"
        else if ( script && script.hasProperty(property)) {
            return script.getProperty(property?.toString())
        }

        if( undef )
        // so give a chance to the bash interpreted to evaluate it
            return '$' + property
        else
            throw new MissingPropertyException("Unknown variable '$property' -- Make sure you didn't misspell it or define somewhere in the script before use it")

    }

    Object invokeMethod(String name, Object args) {
        script.invokeMethod(name, args)
    }

    public getProperty( String name ) {
        get((String)name)
    }

    public void setProperty( String name, def value ) {
        put(name, value)
    }

    @Override
    public put(String property, Object newValue) {
        holder.put(property, newValue)
    }

    /**
     * The the delegate object to the file specified. It takes care to converts {@code Path} objects
     * (that are not serializable) to objects of type {@code SafePath}
     *
     * @param contextFile The file where store the {@code DelegateMap} instance
     */
    def void save( Path contextFile ) {
        try {
            def copy = safeMap(holder)
            contextFile.bytes = SerializationUtils.serialize(copy)
        }
        catch( SerializationException e ) {
            log.warn "Cannot serialize context map. Cause: ${e.cause} -- Resume will not work on this process"
            log.trace "Serialization failed for map: \n${dumpMap(holder)}\n", e
        }
    }

    /**
     * Replaces the entry which value is a {@code Path} with a {@code SafePath}
     *
     * @param map The map to
     * @return
     */
    @PackageScope
    static Serializable safeMap( Map map ) {
        def result = [:]
        map.each { key, value -> result[key] = safeItem(value) }
        return result
    }

    @PackageScope
    static Serializable safeItem( value ) {

        if( value instanceof Serializable || value == null ) {
            return value
        }

        if( value instanceof Path ) {
            return new SafePath(value)
        }

        if( value instanceof BlankSeparatedList ) {
            return new SafeList(value)
        }

        return value
    }

    @PackageScope
    static String dumpMap( Map map ) {
        def result = []
        result << "[ "
        map.each { key, value -> result << "  '$key':[${value?.class?.name}] = ${value}" }
        result << "]"
        return result.join('\n')
    }

    /**
     * Read the context map from the file specified
     *
     * @param processor The current {@code TaskProcessor}
     * @param contextFile The file used to store the context map
     * @return A new {@code DelegateMap} instance holding the values read from map file
     */
    static DelegateMap read( TaskProcessor processor, Path contextFile ) {
        def map = (Map)SerializationUtils.deserialize( contextFile.bytes )
        map.each {  key, value ->
            if( value instanceof SafePath )
                map[key] = value.toPath()
            else if( value instanceof SafeList )
                map[key] = value.toList()
        }
        new DelegateMap(processor, map)
    }


    /**
     * Wrapper class to serialize a Path storing its URI representation
     *
     * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
     */
    @Slf4j
    @EqualsAndHashCode
    static final class SafePath implements Serializable {

        private static final long serialVersionUID = - 2333294071875114518L ;

        String scheme
        String path

        SafePath(Path target) {
            scheme = target.getFileSystem().provider().getScheme()
            path = target.toString()
        }

        public Path toPath() {

            if( "file".equalsIgnoreCase(scheme) ) {
                return FileSystems.getDefault().getPath(path)
            }

            // try to find provider
            for (FileSystemProvider provider: FileSystemProvider.installedProviders()) {
                if (provider.getScheme().equalsIgnoreCase(scheme)) {
                    return provider.getPath(new URI(path))
                }
            }

            throw new FileSystemNotFoundException("Provider \"$scheme\" not installed");
        }
    }

    /**
     * Holds a list that can be safely serialised
     */
    static final class SafeList implements Serializable {

        private static final long serialVersionUID = 2310939302448910665L ;

        final List target

        SafeList( BlankSeparatedList list ) {
            target = new ArrayList<>(list.size())
            list.eachWithIndex { entry, int i ->
                target[i] = safeItem(entry)
            }
        }

        /**
         * @return back the list as a {@code BlankSeparatedList} instance
         */
        BlankSeparatedList toList() {
            def copy = new ArrayList(target.size())
            target.each {

                if( it instanceof SafePath )
                    copy.add(it.toPath())
                else
                    copy.add(it)
            }

            return new BlankSeparatedList(copy)
        }

    }


}