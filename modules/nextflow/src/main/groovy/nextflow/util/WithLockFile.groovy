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

import java.nio.channels.FileLock

/**
 * File with a file lock.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class WithLockFile extends File {

    WithLockFile(String filepath){
        super(filepath)
    }

    /**
     * Apply the given action by using a file lock
     *
     * @param action The closure implementing the action to be executed with a file lock
     * @return The value returned by the action closure
     */
    protected withFileLock(Closure action) {

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
}
