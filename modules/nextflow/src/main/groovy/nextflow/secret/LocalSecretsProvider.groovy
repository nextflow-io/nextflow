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

package nextflow.secret

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import com.google.gson.Gson
import com.google.gson.GsonBuilder
import com.google.gson.reflect.TypeToken
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.util.Escape
/**
 * Implements a secrets store that saves secrets into a JSON file save into the
 * nextflow home. The file can be relocated using the env variable {@code NXF_SECRETS_FILE}.
 *
 * The file must only the accessible to the running user and it should accessible via a shared
 * file system when using a batch scheduler based executor
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class LocalSecretsProvider implements SecretsProvider, Closeable {

    final private static String ONLY_OWNER_PERMS = 'rw-------'

    private Map<String,String> env = SysEnv.get()

    private Map<String, Secret> secretsMap

    private Path storeFile

    private boolean modified

    LocalSecretsProvider() {
        storeFile = makeStoreFile()
        log.debug "Secrets store: $storeFile"
        secretsMap = makeSecretsSet()
    }

    protected Path makeStoreFile() {
        final name = env.NXF_SECRETS_FILE ?: 'secrets/store.json'
        final secretFile = name.startsWith('/')
                ? Paths.get(name)
                : Const.APP_HOME_DIR.resolve(name)
        final path = secretFile.parent
        if( path && !path.exists() && !path.mkdirs() )
            throw new IllegalStateException("Cannot create directory '${path}' -- make sure a file with the same name doesn't already exist and you have write permissions")
        return secretFile
    }

    protected Map<String,Secret> makeSecretsSet() {
        new TreeMap<String, Secret>(new Comparator<String>() {
            @Override
            int compare(String o1,String o2) {
                return o1 <=> o2
            }
        })
    }

    @Override
    boolean activable() { true }

    /**
     * Load secrets from the stored json file
     * @return
     */
    @Override
    LocalSecretsProvider load() {
        // load secrets field
        final allSecrets = loadSecrets()
        for( Secret secret : allSecrets )
            secretsMap.put(secret.name, secret)
        return this
    }

    @Override
    Secret getSecret(String name) {
        return secretsMap.get(name)
    }

    @Override
    void putSecret(String name, String value) {
        SecretsHelper.checkName(name)
        putSecret(new SecretImpl(name, value))
    }

    void putSecret(Secret secret) {
        if( !secret.getName() )
            throw new IllegalStateException("Missing secret name")
        if( !secret.getValue() )
            throw new IllegalStateException("Missing secret value")
        secretsMap.put(secret.name, secret)
        modified = true
    }

    @Override
    void removeSecret(String name) {
        if( secretsMap.containsKey(name) ) {
            secretsMap.remove(name)
            modified = true
        }
    }

    Set<String> listSecretsNames() {
        return secretsMap.keySet()
    }

    protected List<Secret> loadSecrets() {
        if( !storeFile.exists() ) {
            return Collections.<Secret>emptyList()
        }
        // make sure the permissions are valid
        if( storeFile.getPermissions() != ONLY_OWNER_PERMS )
            throw new AbortOperationException("Invalid permissions for secret store file: $storeFile - It is required that your secret store file is NOT accessible by others.")
        // read the JSON secrets file
        final type = new TypeToken<ArrayList<SecretImpl>>(){}.getType()
        return new Gson().fromJson(storeFile.getText('utf-8'), type)
    }

    protected void storeSecrets(Collection<Secret> secrets) {
        assert secrets != null
        final parent = storeFile.getParent()
        if( !parent.exists() && !parent.mkdirs() )
            throw new IOException("Unable to create directory: $parent -- Check file system permissions" )
        // save the secrets as JSON file
        final json = new GsonBuilder().setPrettyPrinting().create().toJson(secrets)
        Files.write(storeFile, json.getBytes('utf-8'))
        // allow only the current user to read-write the file
        storeFile.setPermissions(ONLY_OWNER_PERMS)
    }

    @Override
    void close() throws IOException {
        if( modified )
            storeSecrets(secretsMap.values())
    }

    @Override
    String getSecretsEnv(List<String> secretNames) {
        if( !secretNames )
            return null
        // find out if any required secret is missing
        final missing = secretNames - this.listSecretsNames()
        if( missing ) {
            final names = missing.collect(it -> "'$it'").join(', ')
            final msg = missing.size()==1
                    ? "Required secret is missing: $names"
                    : "Required secrets are missing: $names"
            throw new ProcessUnrecoverableException(msg)
        }
        final filter = secretNames.collect(it -> "-e '$it=.*'").join(' ')
        final tmp = makeTempSecretsFile()
        // mac does not allow source an anonymous pipe
        // https://stackoverflow.com/a/32596626/395921
        return tmp ? "source /dev/stdin <<<\"\$(cat <(grep -w $filter $tmp))\"" : null
    }

    @Deprecated
    String getSecretsEnv() {
        return getSecretsEnv(null)
    }

    /**
     * Creates temporary file containing secrets to be included in the task environment
     *
     * Note: the file is expected to be located in the user home and accessible from a shared
     * file system when using a batch scheduler. The file remove once the execution completes.
     *
     * @return The secret shell formatted file of secrets
     */
    @Memoized // cache this method to avoid creating multiple copies of the same file
    protected Path makeTempSecretsFile() {
        if( !secretsMap )
            return null

        final name = ".nf-${UUID.randomUUID().toString()}.secrets"
        final path = storeFile.parent.resolve(name)
        if( path.exists() ) {
            // make sure the file can only be accessed by the owner user
            path.setPermissions(ONLY_OWNER_PERMS)
            // remove it on completion
            path.toFile().deleteOnExit()
            return path
        }

        def result = ''
        for( Secret s : secretsMap.values() ) {
            result += /export ${s.name}="${Escape.variable(s.value)}"/
            result += '\n'
        }
        Files.createFile(path)
        // make sure the file can only be accessed by the owner user
        path.setPermissions(ONLY_OWNER_PERMS)
        path.text = result
        // remove it on completion
        path.toFile().deleteOnExit()
        // return the temp path
        return path
    }
}
