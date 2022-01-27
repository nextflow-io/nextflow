/*
 * Copyright 2021, Sage-Bionetworks
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

package nextflow.secret

import org.pf4j.ExtensionPoint

/**
 * Model the secret provider interface
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface SecretsProvider extends ExtensionPoint, Closeable {

    boolean activable()

    /**
     * Loads the secretes from the underlying provider store
     *
     * @return the {@link SecretsProvider} instance itself
     */
    SecretsProvider load()

    /**
     * Retrieve a {@code Secret} by the name from the secret provider
     *
     * @param name The secret name
     * @return The {@code Secret} instance with the given name or {@code null} otherwise
     */
    Secret getSecret(String name)

    /**
     * Store a {@code Secret} in the secret provider
     *
     * @param name The {@code Secret} instance to store
     */
    void putSecret(String name, String value)

    /**
     * Remove the secret with the given name
     *
     * @param name The name of the secret to be removed
     *
     */
    void removeSecret(String name)

    /**
     * Retrieve the collection of available secret names
     * @return The of secret names
     */
    Set<String> listSecretsNames()

    /**
     * @return A shell snippet defining the secrets as environment variables
     */
    String getSecretsEnv(List<String> secretNames)

    @Deprecated
    String getSecretsEnv()

}
