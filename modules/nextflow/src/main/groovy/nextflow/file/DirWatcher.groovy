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

package nextflow.file

import static java.nio.file.StandardWatchEventKinds.ENTRY_CREATE
import static java.nio.file.StandardWatchEventKinds.ENTRY_DELETE
import static java.nio.file.StandardWatchEventKinds.ENTRY_MODIFY
import static java.nio.file.StandardWatchEventKinds.OVERFLOW

import java.nio.file.ClosedWatchServiceException
import java.nio.file.FileSystem
import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.PathMatcher
import java.nio.file.SimpleFileVisitor
import java.nio.file.WatchEvent
import java.nio.file.WatchKey
import java.nio.file.WatchService
import java.nio.file.attribute.BasicFileAttributes

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j

/**
 * Watch the content of a directory for file system events
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class DirWatcher {

    static private Map<String,WatchEvent.Kind<Path>> EVENT_MAP = [
            'create':ENTRY_CREATE,
            'delete':ENTRY_DELETE,
            'modify':ENTRY_MODIFY
    ]

    private Path base

    private FileSystem fs

    private String pattern

    private boolean skipHidden

    private String folder

    private String syntax

    private WatchEvent.Kind<Path>[] watchEvents

    private PathMatcher dirMatcher

    private PathMatcher fileMatcher

    private String fileRule

    private String dirRule

    private Closure onNext

    private Closure onComplete

    private WatchService watcher

    private HashMap<WatchKey,Path> watchedPaths

    private Thread thread

    private volatile boolean terminated

    DirWatcher(String syntax, String folder, String pattern, boolean skipHidden, String events, FileSystem fs) {
        assert syntax in ['regex','glob']
        assert folder.endsWith("/")
        log.debug "Watch service for path=$folder; syntax=$syntax; pattern=$pattern; skipHidden=$skipHidden; events=$events"

        this.syntax = syntax
        this.folder = folder
        this.pattern = pattern
        this.skipHidden = skipHidden
        this.base = fs.getPath(folder)
        this.watchEvents = stringToWatchEvents(events)
        this.fs = fs

        this.fileRule = "$syntax:${folder}${pattern}"
        this.fileMatcher = FileHelper.getPathMatcherFor(fileRule, base.fileSystem)

        def p = pattern.indexOf('/')
        if( p>0 ) {
            dirRule = "$syntax:${pattern.substring(0,p)}"
            dirMatcher = FileHelper.getPathMatcherFor(dirRule, base.fileSystem)
        }
        else if( pattern.contains('**') ) {
            dirRule = "$syntax:**"
            dirMatcher = FileHelper.getPathMatcherFor(dirRule, base.fileSystem)
        }

    }

    DirWatcher setOnComplete(Closure action) {
        this.onComplete = action
        return this
    }

    void terminate() {
        terminated = true
        watcher?.close()
        thread?.join()
    }

    void apply( Closure onNext ) {
        this.onNext = onNext

        if( !base.isDirectory() ) {
            log.warn "Cannot watch a not existing directory: $base -- Make sure that path exists and it is a directory"
            onComplete?.call()
            return
        }

        this.watcher = base.getFileSystem().newWatchService()
        this.watchedPaths = new HashMap<WatchKey,Path>()

        thread = Thread.startDaemon {
            try {
                apply0()
            }
            finally {
                onComplete?.call()
            }
        }

    }

    private void register (Path folder, PathMatcher matcher) {
        assert folder
        assert matcher

        Files.walkFileTree(folder, new SimpleFileVisitor<Path>() {
            @Override
            FileVisitResult preVisitDirectory(Path path, BasicFileAttributes attrs) throws IOException
            {
                def dir = base.relativize(path)
                if( matcher.matches(dir) ) {
                    register(path)
                }
                else {
                    log.trace "Skip watcher for dir=$dir; matcher-rule=$dirRule"
                }
                return FileVisitResult.CONTINUE;
            }
        })

    }

    private void register(Path folder) {
        assert folder
        def key = folder.register(watcher, watchEvents)
        watchedPaths[key] = folder
        log.trace "Register watcher for dir=$folder; events=$watchEvents"
    }

    private void apply0() {

        dirMatcher ? register(base,dirMatcher) : register(base)

        while( !terminated ) {
            // wait for key to be signaled
            try {
                WatchKey key = watcher.take()
                if( !key )
                    continue
                def path = watchedPaths.get(key)
                if( !path ) {
                    log.error "Unknown file watch key: $key"
                    continue
                }

                for (WatchEvent<?> event: key.pollEvents()) {
                    WatchEvent.Kind<?> kind = event.kind();

                    if( kind == OVERFLOW ) {
                        log.debug "Watcher event > OVERFLOW; path=$path"
                        continue
                    }

                    // The filename is the context of the event.
                    Path fileName = (event as WatchEvent<Path>).context();
                    log.trace "Watcher event > $kind; fileName=$fileName; path=$path"
                    Path target = path.resolve(fileName)

                    if (fileMatcher.matches(target) && ( !skipHidden || !Files.isHidden(target) )) {
                        log.trace "Watcher match > $kind; fileName=$fileName; rule=$fileRule"
                        onNext.call(target)
                    }

                    if( kind == ENTRY_CREATE && dirMatcher && Files.isDirectory(target, LinkOption.NOFOLLOW_LINKS) ) {
                        register(target, dirMatcher)
                    }
                }

                // Reset the key -- this step is critical if you want to
                // receive further watch events.  If the key is no longer valid,
                // the directory is inaccessible so exit the loop.
                boolean valid = key.reset();
                if (!valid) {
                    break;
                }
            }
            catch (ClosedWatchServiceException e ) {
                log.debug "Closed watch service for path: $base"
                break
            }
            catch (InterruptedException e) {
                log.debug "Interrupted watch service for path: $base"
                break
            }
            catch (Exception e) {
                log.debug "Exception while watching path: $base", e
                break
            }

        }
    }


    /**
     * Converts a comma separated events string to the corresponding {@code WatchEvent.Kind} instances
     *
     * @param events the list of events to watch
     * @return
     */
    @PackageScope
    static WatchEvent.Kind<Path>[] stringToWatchEvents(String events = null){
        def result = []
        if( !events )
            result << ENTRY_CREATE

        else {
            events.split(',').each { String it ->
                def ev = it.trim().toLowerCase()
                def val = EVENT_MAP[ev]
                if( !val )
                    throw new IllegalArgumentException("Invalid watch event: $it -- Valid values are: ${EVENT_MAP.keySet().join(', ')}")
                result << val
            }
        }

        result as WatchEvent.Kind<Path>[]
    }

}
