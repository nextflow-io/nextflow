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

import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.Paths

import org.codehaus.groovy.runtime.InvokerHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Deprecated
class StagePath {

    final Path target

    StagePath( Path source ) {
        this.target = source
    }

    StagePath( String fileName ) {
        this.target = Paths.get(fileName)
    }

    Object getProperty( String name ) {
        InvokerHelper.getProperty(target,name)
    }

    void setProperty(String property, Object newValue) {
        InvokerHelper.setProperty(target, property, newValue)
    }

    Object invokeMethod(String name, Object args) {
        InvokerHelper.invokeMethod(target, name, args)
    }

    Path getFileName() {
        target.getFileName()
    }

    Path toRealPath(LinkOption... options) {
        return target
    }

    String toString() {
        target.getFileName().toString()
    }
}
