/*
 * Copyright 2020, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
import org.apache.commons.io.monitor.FileAlterationListener
import org.apache.commons.io.monitor.FileAlterationMonitor
import org.apache.commons.io.monitor.FileAlterationObserver

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class DirMonitor implements DirChangeListener<DirMonitor>, FileAlterationListener {

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

    protected DirMonitor() { }

    DirMonitor(String syntax, String folder, String pattern, boolean skipHidden, String events, FileSystem fs) {
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

        this.observer = new FileAlterationObserver(base.toFile());
        this.monitor = new FileAlterationMonitor(1_000)
        observer.addListener(this)
        monitor.addObserver(observer)
    }

    void apply( Closure onNext ) {
        this.onNext = onNext
        if( !base.isDirectory() ) {
            log.warn "Cannot watch a not existing directory: $base -- Make sure that path exists and it is a directory"
            onComplete?.call()
            return
        }

        monitor.start()
    }

    @Override
    DirMonitor setOnComplete(Closure action) {
        this.onComplete = action
        return this
    }

    private notify(WatchEvent.Kind<Path> kind, File file) {
        final target = file.toPath()
        if (fileMatcher.matches(target) && ( !skipHidden || !Files.isHidden(target) )) {
            log.trace "Watcher match > $kind; target=$target; rule=$fileRule"
            onNext.call(target)
        }
    }

    void terminate() {
        monitor.stop()
    }

    // -- file listeners

    @Override
    void onFileCreate(File file) {
        if( StandardWatchEventKinds.ENTRY_CREATE in watchEvents ) {
            notify(StandardWatchEventKinds.ENTRY_CREATE, file)
        }
    }

    @Override
    void onFileChange(File file) {
        if( StandardWatchEventKinds.ENTRY_MODIFY in watchEvents ) {
            notify(StandardWatchEventKinds.ENTRY_MODIFY, file)
        }
    }

    @Override
    void onFileDelete(File file) {
        if( StandardWatchEventKinds.ENTRY_DELETE in watchEvents ) {
            notify(StandardWatchEventKinds.ENTRY_DELETE, file)
        }
    }


    @Override
    void onStart(FileAlterationObserver observer) {

    }

    @Override
    void onDirectoryCreate(File directory) {

    }

    @Override
    void onDirectoryChange(File directory) {

    }

    @Override
    void onDirectoryDelete(File directory) {

    }

    @Override
    void onStop(FileAlterationObserver observer) {

    }
}
