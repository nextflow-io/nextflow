/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.plugin

import java.nio.file.Files
import java.nio.file.Path

import com.google.gson.Gson
import com.google.gson.GsonBuilder
import com.google.gson.JsonSyntaxException
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.util.logging.Slf4j

/**
 * Model the {@code plugins.lock} file. It holds a format version number and a map
 * of plugin fully-qualified ids (ie. {@code id@version}) to the corresponding
 * {@link Entry} carrying the {@code sha512} checksum of the plugin archive.
 *
 * The file is serialised as pretty-printed JSON with a stable (sorted) key order to
 * keep diffs deterministic.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@EqualsAndHashCode(includeFields = true)
@ToString(includeNames = true, includeFields = true)
class PluginLockFile {

    /**
     * The current lock file format version
     */
    static final int CURRENT_VERSION = 1

    /**
     * Represent a single locked plugin entry
     */
    @EqualsAndHashCode
    @ToString(includeNames = true)
    static class Entry {
        String sha512

        Entry() {}

        Entry(String sha512) {
            this.sha512 = sha512
        }
    }

    private int version = CURRENT_VERSION

    private Map<String,Entry> plugins = new TreeMap<String,Entry>()

    int getVersion() { version }

    void setVersion(int value) { this.version = value }

    /**
     * @return An immutable view of the locked plugin entries, keyed by fully-qualified id
     */
    Map<String,Entry> getEntries() {
        return Collections.unmodifiableMap(plugins)
    }

    /**
     * Lookup a locked entry by its fully-qualified id ie. {@code id@version}.
     *
     * @param fqid The plugin fully-qualified id
     * @return The corresponding {@link Entry} or {@code null} if not present
     */
    Entry getEntry(String fqid) {
        return plugins.get(fqid)
    }

    /**
     * Add or update a locked entry.
     *
     * @param fqid The plugin fully-qualified id ie. {@code id@version}
     * @param entry The {@link Entry} to associate with the given id
     * @return The object itself to enable method chaining
     */
    PluginLockFile addEntry(String fqid, Entry entry) {
        if( !fqid )
            throw new IllegalArgumentException("Plugin lock entry id cannot be empty")
        if( entry == null )
            throw new IllegalArgumentException("Plugin lock entry cannot be null")
        plugins.put(fqid, entry)
        return this
    }

    /**
     * @return {@code true} when no plugin entries are held
     */
    boolean isEmpty() {
        return plugins.isEmpty()
    }

    /**
     * Serialise this lock file as pretty-printed JSON with a stable key order.
     *
     * @param path The target file path
     */
    void write(Path path) {
        final json = gson0().toJson(toModel())
        Files.write(path, json.getBytes('UTF-8'))
    }

    private Map<String,Object> toModel() {
        // use a plain map so that only 'version' and 'plugins' are emitted, with
        // plugins held in a TreeMap to guarantee a stable, sorted key order
        final result = new LinkedHashMap<String,Object>()
        result.put('version', version)
        result.put('plugins', new TreeMap<String,Entry>(plugins))
        return result
    }

    private static Gson gson0() {
        return new GsonBuilder().setPrettyPrinting().create()
    }

    /**
     * Read and parse a {@code plugins.lock} file.
     *
     * @param path The lock file path
     * @return A {@link PluginLockFile}; an empty (dormant) instance when the file does not exist
     * @throws IllegalStateException when the file content cannot be parsed
     */
    static PluginLockFile read(Path path) {
        if( path == null || !Files.exists(path) ) {
            log.debug "Plugins lock file does not exist: $path - returning an empty lock"
            return new PluginLockFile()
        }

        final text = new String(Files.readAllBytes(path), 'UTF-8')
        // a blank or freshly `touch`ed file is a valid, empty lock (bootstrap case)
        if( !text.trim() )
            return new PluginLockFile()
        try {
            final model = gson0().fromJson(text, ModelBean)
            if( model == null )
                throw new IllegalStateException("Invalid plugins lock file - empty content: $path")
            final result = new PluginLockFile()
            result.version = model.version
            if( model.plugins )
                result.plugins.putAll(model.plugins)
            return result
        }
        catch( JsonSyntaxException e ) {
            throw new IllegalStateException("Invalid plugins lock file - malformed JSON: $path", e)
        }
    }

    /**
     * Deserialization bean matching the on-disk JSON structure
     */
    static class ModelBean {
        int version
        Map<String,Entry> plugins
    }

}
