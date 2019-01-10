/*
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
 */

package nextflow.util

import java.nio.channels.FileLock
import java.nio.file.Path
import java.text.DateFormat
import java.text.SimpleDateFormat

import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
/**
 * Manages the history file containing the last 1000 executed commands
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class HistoryFile extends File {

    public static final String FILE_NAME = '.nextflow/history'

    @Lazy
    public static final HistoryFile DEFAULT = { def f=new HistoryFile(); f.parentFile?.mkdirs(); return f } ()

    private static final DateFormat TIMESTAMP_FMT = new SimpleDateFormat('yyyy-MM-dd HH:mm:ss')

    private static final VAL_A = (int)('a' as char)
    private static final VAL_F = (int)('f' as char)
    private static final VAL_0 = (int)('0' as char)
    private static final VAL_9 = (int)('9' as char)

    private HistoryFile() {
        super(FILE_NAME)
    }

    HistoryFile(File file) {
        super(file.toString())
    }

    HistoryFile(Path file) {
        super(file.toString())
    }

    void write( String name, UUID key, String revisionId, args, Date date = null ) {
        assert key
        assert args != null

        withFileLock {
            def timestamp = date ?: new Date()
            def value = args instanceof Collection ? args.join(' ') : args
            this << new Record(timestamp: timestamp, runName: name, revisionId: revisionId, sessionId: key, command: value).toString() << '\n'
        }
    }

    void update( String name, boolean success, Date when = null) {
        assert name

        try {
            withFileLock { update0(name, success ? 'OK' : 'ERR', when) }
        }
        catch( Throwable e ) {
            log.warn "Can't update history file: $this",e
        }
    }

    private void update0( String name, String status, Date when ) {

        long ts = when?.time ?: System.currentTimeMillis()
        def newHistory = new StringBuilder()

        this.readLines().each { line ->
            try {
                def current = line ? Record.parse(line) : null
                if( current?.runName == name ) {
                    current.duration = new Duration( ts - current.timestamp.time )
                    current.status = status
                    newHistory << current.toString() << '\n'
                }
                else {
                    newHistory << line << '\n'
                }
            }
            catch( IllegalArgumentException e ) {
                log.warn("Can't read history file: $this", e)
            }
        }

        // rewrite the history content
        this.setText(newHistory.toString())
    }

    Record getLast() {
        if( !exists() || empty() ) {
            return null
        }

        def line = readLines()[-1]
        try {
            line ? Record.parse(line) : null
        }
        catch( IllegalArgumentException e ) {
            log.warn("Can't read history file: $this", e)
            return null
        }
    }

    void print() {

        if( empty() ) {
            System.err.println '(no history available)'
        }
        else {
            println this.text
        }

    }

    /**
     * Check if a session ID exists in the history file
     *
     * @param uuid A complete UUID string or a prefix of it
     * @return {@code true} if the UUID is found in the history file or {@code false} otherwise
     */
    boolean checkExistsById( String uuid ) {
        findById(uuid)
    }

    boolean checkExistsById( UUID uuid ) {
        findById(uuid.toString())
    }

    boolean checkExistsByName( String name ) {
        try {
            return withFileLock { getByName(name) != null }
        }
        catch( AbortOperationException e ) {
            return false
        }
    }


    /**
     * Lookup a session ID given a `run` name string
     *
     * @param name A name of a pipeline run
     * @return The session ID string associated to the `run` or {@code null} if it's not found
     */
    Record getByName(String name) {
        if( !exists() || empty() ) {
            return null
        }

        def result = readLines().findResult {  String line ->
            try {
                def current = line ? Record.parse(line) : null
                if( current?.runName == name )
                    return current
            }
            catch( IllegalArgumentException e ) {
                log.warn("Can't read history file: $this", e)
                return null
            }
        }

        return result
    }

    List<Record> findById(String id) {
        if( !exists() || empty() ) {
            return null
        }

        def found = (List<Record>)this.readLines().findResults { String line ->
            try {
                def current = line ? Record.parse(line) : null
                if( current && current.sessionId.toString().startsWith(id) ) {
                    return current
                }
            }
            catch( IllegalArgumentException e ) {
                log.warn("Can't read history file: $this", e)
                return null
            }
        }

        def results = found.unique()

        // check the multiple results belong to the same session
        def sessions = results.collect { it.sessionId } .unique()
        if( sessions.size()>1 ) {
            String message = 'Which session ID do you mean?\n'
            sessions.each { message += "    $it\n" }
            throw new AbortOperationException(message)
        }

        return results
    }


    List<Record> findByIdOrName( String str ) {
        if( str == 'last' ) {
            def entry = getLast()
            return entry ? [entry] : Collections.emptyList()
        }

        if( isUuidString(str) )
            return findById(str)

        else {
            def entry = getByName(str)
            return entry ? [entry] : Collections.emptyList()
        }

    }

    @PackageScope
    static boolean isUuidChar(char ch) {
        if( ch == '-' as char )
            return true

        final x = (ch as int)

        if(  x >= VAL_0 && x <= VAL_9 )
            return true

        if( x >= VAL_A && x <= VAL_F )
            return true

        return false
    }

    @PackageScope
    static boolean isUuidString(String str) {
        for( int i=0; i<str.size(); i++ )
            if( !isUuidChar(str.charAt(i)))
                return false
        return true
    }

    List<Record> findAll() {
        if( !exists() || empty() ) {
            return Collections.emptyList()
        }

        def results = this.readLines().findResults {  String line ->
            try {
                line ? Record.parse(line) : null
            }
            catch( IllegalArgumentException e ) {
                log.warn("Can't read history file: $this", e)
                return null
            }
        }

        return results.unique()
    }

    List<Record> findBefore(String idOrName) {
        def matching = findByIdOrName(idOrName)
        if( !matching )
            return Collections.emptyList()

        def firstMatch = false

        return findAll().findResults {
            if( it==matching[0] ) {
                firstMatch = true
                return null
            }

            !firstMatch ? it : null
        }
    }

    List<Record> findAfter(String idOrName) {
        def matching = findByIdOrName(idOrName)
        if( !matching )
            return Collections.emptyList()

        def firstMatch = false
        return findAll().findResults {
            if( it==matching[-1] ) {
                firstMatch = true
                return null
            }

            firstMatch ? it : null
        }
    }

    List<Record> findBut(String idOrName) {
        final matching = findByIdOrName(idOrName)
        final result = findAll()
        result.removeAll(matching)
        return result
    }

    void eachRow(Closure action) {
        this.eachLine { String line ->
            if( line ) {
                action.call( Record.parse(line).toList() )
            }
        }
    }

    void deleteEntry(Record entry) {
        withFileLock {
            deleteEntry0(entry)
        }
    }

    void deleteEntry0(Record entry) {

        def newHistory = new StringBuilder()

        this.readLines().each { line ->
            try {
                def current = line ? Record.parse(line) : null
                if( current != entry ) {
                    newHistory << line << '\n'
                }
            }
            catch( IllegalArgumentException e ) {
                log.warn("Can't read history file: $this", e)
            }
        }

        // rewrite the history content
        this.setText(newHistory.toString())

    }

    @EqualsAndHashCode(includes = 'runName,sessionId')
    static class Record {
        Date timestamp
        Duration duration
        String runName
        String status
        String revisionId
        UUID sessionId
        String command

        Record(String sessionId, String name=null) {
            this.runName = name
            this.sessionId = UUID.fromString(sessionId)
        }

        Record(UUID sessionId, String name=null) {
            this.runName = name
            this.sessionId = sessionId
        }

        protected Record() {}

        List<String> toList() {
            def line = new ArrayList<String>(7)
            line << (timestamp ? TIMESTAMP_FMT.format(timestamp) : '-')
            line << (duration ? duration.toString() : '-')
            line << (runName ?: '-')
            line << (status ?: '-')
            line << (revisionId ?: '-')
            line << (sessionId.toString())
            line << (command ?: '-')
        }

        @Override
        String toString() {
            toList().join('\t')
        }

        static Record parse(String line) {
            def cols = line.tokenize('\t')
            if( cols.size() == 2 )
                return new Record(cols[0])

            if( cols.size()==7 ) {

                return new Record(
                        timestamp: TIMESTAMP_FMT.parse(cols[0]),
                        duration: cols[1] && cols[1] != '-' ? Duration.of(cols[1]) : null,
                        runName: cols[2],
                        status: cols[3] && cols[3] != '-' ? cols[3] : null,
                        revisionId: cols[4],
                        sessionId: UUID.fromString(cols[5]),
                        command: cols[6]
                )
            }

            throw new IllegalArgumentException("Not a valid history entry: `$line`")
        }
    }

    /**
     * Apply the given action by using a file lock
     *
     * @param action The closure implementing the action to be executed with a file lock
     * @return The value returned by the action closure
     */
    private withFileLock(Closure action) {

        def rnd = new Random()
        long ts = System.currentTimeMillis()
        String parent = this.parent ?: new File('.').absolutePath
        def file = new File(parent, "${this.name}.lock".toString())
        def fos = new FileOutputStream(file)
        try {
            Throwable error
            FileLock lock = null

            try {
                while( true ) {
                    lock = fos.getChannel().tryLock()
                    if( lock ) break
                    if( System.currentTimeMillis() - ts < 1_000 )
                        sleep rnd.nextInt(75)
                    else {
                        error = new IllegalStateException("Can't lock file: ${this.absolutePath} -- Nextflow needs to run in a file system that supports file locks")
                        break
                    }
                }
                if( lock ) {
                    return action.call()
                }
            }
            catch( Exception e ) {
                return action.call()
            }
            finally {
                if( lock?.isValid() ) lock.release()
            }

            if( error ) throw error
        }
        finally {
            fos.closeQuietly()
            file.delete()
        }
    }

    Set<String> findAllRunNames() {
        findAll().findResults{ it.runName }
    }

    String generateNextName() {
        return withFileLock { NameGenerator.next(findAllRunNames()) }

    }
}
