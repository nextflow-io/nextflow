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

import groovy.util.logging.Slf4j
import nextflow.script.BaseScript
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

    private boolean cacheable


    DelegateMap( TaskProcessor processor, Map holder = null ) {
        this.holder = holder ?: [:]
        this.script = processor.ownerScript
        this.undef = processor.taskConfig.getUndef()
        this.cacheable = processor.isCacheable()
        this.name = processor.name
    }

    @Deprecated
    DelegateMap(BaseScript script, boolean undef = false, cacheable = true) {
        this.script = script
        this.holder = [:]
        this.undef = undef
        this.cacheable = cacheable
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

    public getProperty( String name ) {
        get((String)name)
    }

    public void setProperty( String name, def value ) {
        put(name, value)
    }

    @Override
    public put(String property, Object newValue) {
        if( newValue instanceof Path ) {
            newValue = new SerializablePath(newValue)
            log.debug "Swapping to serializable path in context map: $newValue"
        }
        else if( !(newValue instanceof Serializable) && cacheable ) {
            log.warn "Variable '$property' does not implement the Java Serializable interface -- Resume feature will not work for process '$name'"
        }
        holder.put(property, newValue)
    }


    def void save( Path contextFile ) {
        try {
            contextFile.bytes = SerializationUtils.serialize((Serializable)holder)
        }
        catch( SerializationException e ) {
            log.warn "Cannot serialize context map. Cause: ${e.cause} -- Resume will not work on this process"
        }
    }

    static DelegateMap read( TaskProcessor processor, Path contextFile ) {
        def map = (Map)SerializationUtils.deserialize( contextFile.bytes )
        new DelegateMap(processor, map)
    }



    /**
     * Wrapper class to serialize a Path storing its URI representation
     *
     * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
     */
    @Slf4j
    static class SerializablePath implements Serializable, Path {

        @Delegate
        Path target

        SerializablePath(Path path) {
            this.target = path
        }


        private void writeObject(ObjectOutputStream stream) throws IOException {
            final scheme = target.getFileSystem().provider().getScheme()
            final path = 'file'.equalsIgnoreCase(scheme) ? target.toString() : target.toUri().toString()

            log.trace "Serializing path object -- scheme: $scheme; path: $path"
            stream.writeObject(scheme)
            stream.writeObject(path)
        }

        private void readObject(ObjectInputStream stream) throws IOException, ClassNotFoundException {
            final String scheme = stream.readObject()
            final String path = stream.readObject()
            log.trace "De-serializing path object -- scheme: $scheme; path: $path"

            if( "file".equalsIgnoreCase(scheme) ) {
                target = FileSystems.getDefault().getPath(path)
                return
            }

            // try to find provider
            for (FileSystemProvider provider: FileSystemProvider.installedProviders()) {
                if (provider.getScheme().equalsIgnoreCase(scheme)) {
                    target = provider.getPath(new URI(path))
                    return
                }
            }

            throw new FileSystemNotFoundException("Provider \"" + scheme + "\" not installed");
        }


        boolean equals( Object other ) {
            return target.equals(other)
        }

        int hashCode() {
            return target.hashCode()
        }

        int compareTo(Path other) {
            target.compareTo(other)
        }

        String toString() {
            target.toString()
        }

    }
}