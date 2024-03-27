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
 *
 */

package nextflow.script.bundle

import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.file.FileHelper
import nextflow.util.CacheHelper
import nextflow.util.MemoryUnit

/**
 * Model a module bundle
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Canonical
@CompileStatic
class ResourcesBundle {

    public static MemoryUnit MAX_FILE_SIZE = MemoryUnit.of('1MB')
    public static MemoryUnit MAX_BUNDLE_SIZE = MemoryUnit.of('5MB')

    private Path root
    private LinkedHashMap<String,Path> content = new LinkedHashMap<>(100)
    private Path dockerfile
    private Path singularityfile
    private MemoryUnit maxFileSize = MAX_FILE_SIZE
    private MemoryUnit maxBundleSize = MAX_BUNDLE_SIZE
    private String baseDirectory

    ResourcesBundle(Path root) {
        this.root = root
        this.dockerfile = pathIfExists0(root.resolveSibling('Dockerfile'))
        this.singularityfile = pathIfExists0(root.resolveSibling('Singularityfile'))
    }

    ResourcesBundle withMaxFileSize(MemoryUnit mem) {
        this.maxFileSize = mem
        return this
    }

    ResourcesBundle withBundleSize(MemoryUnit mem) {
        this.maxBundleSize = mem
        return this
    }

    Path getRoot() { root }

    Map<String,Path> content() { content }

    static private Path pathIfExists0(Path path) {
        return path?.exists() ? path : null
    }

    ResourcesBundle withPaths(Collection<Path> paths) {
        this.content = new LinkedHashMap<String,Path>(100)
        long totSize = 0
        for( Path it : paths ) {
            final attrs = Files.readAttributes(it, BasicFileAttributes, LinkOption.NOFOLLOW_LINKS)
            if( attrs.size()>maxFileSize.bytes )
                throw new IllegalArgumentException("Module file size cannot be bigger than $maxFileSize - offending file: $it")
            if( !attrs.isDirectory() && !attrs.isRegularFile() )
                throw new IllegalArgumentException("Module bundle does not allow link files - offending file: $it")
            if( attrs.isRegularFile() ) {
                totSize += attrs.size()
                if( totSize>maxBundleSize.bytes )
                throw new IllegalArgumentException("Module total size cannot exceed $maxBundleSize")
            }

            final relPath = root.relativize(it)
            final name = baseDirectory
                    ? Path.of(baseDirectory).resolve(relPath).toString()
                    : relPath.toString()
            content.put(name, it.normalize())
        }
        return this
    }

    Path getDockerfile() {
        return dockerfile
    }

    Path getSingularityfile() {
        return singularityfile
    }

    Set<Path> getPaths() {
        return new HashSet<Path>(content.values())
    }

    @Deprecated
    List<Path> getPathsList() {
        final result = new ArrayList<Path>(content.size())
        for( String name : getEntries() )
            result.add(path(name))
        return result
    }

    Path path(String name) {
        return content.get(name)
    }

    Set<String> getEntries() {
        return new TreeSet<String>(content.keySet())
    }

    boolean hasEntries() {
        return content.size()
    }

    boolean asBoolean() {
        return content.size() || dockerfile || singularityfile
    }

    /**
     * Creates a {@link ResourcesBundle} object populated with the set of files in the root directory
     *
     * @param bundleRoot
     *      The bundle root path
     * @return
     *      An instance of {@link ResourcesBundle} holding the set of files that are container
     *      in the bundle directory
     */
    static ResourcesBundle scan(Path bundleRoot, Map config=[:]) {
        final result = new ResourcesBundle(bundleRoot)
        if( !bundleRoot.exists() )
            return result
        if( !bundleRoot.isDirectory() ) {
            log.warn "Module bundle location is not a directory path: '$bundleRoot'"
            return result
        }
        // setup config
        if( config.maxFileSize )
            result.maxFileSize = config.maxFileSize as MemoryUnit
        if( config.maxBundleSize )
            result.maxBundleSize = config.maxBundleSize as MemoryUnit
        if( config.baseDirectory )
            result.baseDirectory = config.baseDirectory as String
        
        // load bundle files
        final files = new HashSet(10)
        final pattern = config.filePattern as String ?: '**'
        final opts = [type: 'any', hidden: true, relative: false]
        FileHelper.visitFiles(opts, bundleRoot, pattern) { files.add(it) }
        result.withPaths(files)
        return result
    }

    private List fileMeta(String name, Path file) {
        final attrs = Files.readAttributes(file, BasicFileAttributes)
        final regular = attrs.isRegularFile()
        final meta = [
                name,
                regular ? attrs.size() : 0,
                regular ? md5(file, attrs) : 0,
                Integer.toOctalString(file.getPermissionsMode()) ]
        log.trace "Module bundle entry=$meta"
        return meta
    }

    private String md5(Path path, BasicFileAttributes attrs) {
        return attrs.size() <= MAX_FILE_SIZE.bytes
                ? Files.readAllBytes(path).md5()
                : '0'
    }

    String fingerprint() {
        final allMeta = new ArrayList()
        for( String name : getEntries() ) {
            final file = this.path(name)
            allMeta.add(fileMeta(name,file))
        }
        if( dockerfile ) {
            allMeta.add(fileMeta(dockerfile.name, dockerfile))
        }
        if( singularityfile ) {
            allMeta.add(fileMeta(singularityfile.name, singularityfile))
        }

        return CacheHelper.hasher(allMeta).hash().toString()
    }

    final private static List<String> BIN_PATHS = ['bin','usr/bin','usr/local/bin']

    List<Path> getBinDirs() {
        final result = new ArrayList<Path>(10)
        for( Map.Entry<String,Path> it : content ) {
            if( it.key in BIN_PATHS && Files.isDirectory(it.value) && !result.contains(it.value) )
                result.add(it.value)
        }
        // sort to make order predictable
        Collections.sort(result)
        return result
    }
}
