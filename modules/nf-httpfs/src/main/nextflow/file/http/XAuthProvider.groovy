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

package nextflow.file.http

/**
 * Implements a pluggable authentication provider for {@link XFileSystemProvider}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface XAuthProvider {

    /**
     * Implementing class should check whenever accept and authorise the connection
     *
     * @param connection A {@link URLConnection} object instance to be authorised
     * @return {@code true} if the connection object has been authorised or {@code false} otherwise
     */
    boolean authorize( URLConnection connection )

    boolean refreshToken( URLConnection connection )

}
