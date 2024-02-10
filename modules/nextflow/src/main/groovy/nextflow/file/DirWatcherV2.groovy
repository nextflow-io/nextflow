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

package nextflow.file

import java.nio.file.FileSystem
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.PathMatcher
import java.nio.file.StandardWatchEventKinds
import java.nio.file.WatchEvent

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.util.Duration
import org.apache.commons.io.monitor.FileAlterationListener
import org.apache.commons.io.monitor.FileAlterationMonitor
import org.apache.commons.io.monitor.FileAlterationObserver
/**
 * Implements a directory watcher using polling strategy instead
 * of Java NIO watcher service which is not supported by non-posix
 * native file systems e.g. shared file system, FUSE driver, etc.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class DirWatcherV2 implements DirListener, FileAlterationListener {

    private static final String DEFAULT_INTERVAL = '1s'

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

    private Closure<Path> onNext

    private Closure onComplete

    private FileAlterationObserver observer

    private FileAlterationMonitor monitor

    protected DirWatcherV2() { }

    DirWatcherV2(String syntax, String folder, String pattern, boolean skipHidden, String events, FileSystem fs) {
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

        final interval = System.getenv('NXF_DIRWATCHER_INTERVAL') ?: DEFAULT_INTERVAL
        if( interval!=DEFAULT_INTERVAL )
            log.debug "Dir watcher interval=$interval"
        this.monitor = new FileAlterationMonitor(Duration.of(interval).toMillis())
        this.observer = new FileAlterationObserver(base.toFile());
        observer.addListener(this)
        monitor.addObserver(observer)
    }

    void apply( Closure onNext ) {
        this.onNext = onNext
        if( !base.isDirectory() ) {
            log.warn "Cannot watch a non-existent directory: $base -- Ensure the path is a valid directory"
            onComplete?.call()
            return
        }

        monitor.start()
    }

    @Override
    void onComplete(Closure action) {
        this.onComplete = action
    }

    private notify(WatchEvent.Kind<Path> kind, File file) {
        final target = file.toPath()
        if (fileMatcher.matches(target) && ( !skipHidden || !Files.isHidden(target) )) {
            onNext.call(target)
        }
    }

    void terminate() {
        monitor.stop()
    }

    // -- file listeners

    @Override
    void onFileCreate(File file) {
        log.trace "Dir watcher event onFileCreate=$file"
        if( StandardWatchEventKinds.ENTRY_CREATE in watchEvents ) {
            notify(StandardWatchEventKinds.ENTRY_CREATE, file)
        }
    }

    @Override
    void onFileChange(File file) {
        log.trace "Dir watcher event onFileChange=$file"
        if( StandardWatchEventKinds.ENTRY_MODIFY in watchEvents ) {
            notify(StandardWatchEventKinds.ENTRY_MODIFY, file)
        }
    }

    @Override
    void onFileDelete(File file) {
        log.trace "Dir watcher event onFileDelete=$file"
        if( StandardWatchEventKinds.ENTRY_DELETE in watchEvents ) {
            notify(StandardWatchEventKinds.ENTRY_DELETE, file)
        }
    }

    @Override
    void onDirectoryCreate(File directory) {
        log.trace "Dir watcher event onDirectoryCreate=$directory"
    }

    @Override
    void onDirectoryChange(File directory) {
        log.trace "Dir watcher event onDirectoryChange=$directory"
    }

    @Override
    void onDirectoryDelete(File directory) {
        log.trace "Dir watcher event onDirectoryDelete=$directory"
    }


    @Override
    void onStart(FileAlterationObserver observer) {
        log.trace "Dir watcher event onStart"
    }


    @Override
    void onStop(FileAlterationObserver observer) {
        log.trace "Dir watcher event onStop"
    }
}
