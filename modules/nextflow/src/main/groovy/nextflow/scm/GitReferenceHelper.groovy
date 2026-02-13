/*
 * Copyright 2013-2025, Seqera Labs
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
package nextflow.scm

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.eclipse.jgit.lib.ObjectId
import org.eclipse.jgit.lib.Ref
import org.eclipse.jgit.lib.Repository

/**
 * Class with helper methods for Git references
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class GitReferenceHelper {
    private static final String REMOTE_REFS_ROOT = "refs/remotes/origin/"
    private static final String REMOTE_DEFAULT_HEAD = REMOTE_REFS_ROOT + "HEAD"

    /**
     * Check if a ref has remote changes
     */
    static boolean hasRemoteChange(Ref ref, Map<String, Ref> remote) {
        if( !remote )
            return false

        final remoteRef = remote.get(ref.name)
        if( !remoteRef )
            return false

        return ref.getObjectId() != remoteRef.getObjectId()
    }

    /**
     * Format update information for a remote ref
     */
    static String formatUpdate(Ref remoteRef, int level) {
        final result = new StringBuilder()
        result << 'updates on remote'
        if( level ) {
            result << ' '
            result << formatObjectId(remoteRef.getObjectId(), level == 1)
        }
        return result.toString()
    }

    static boolean isRemoteBranch(Ref ref) {
        return ref.name.startsWith(REMOTE_REFS_ROOT) && ref.name != REMOTE_DEFAULT_HEAD
    }

    static formatObjectId(ObjectId obj, boolean human) {
        return human ? obj.name.substring(0, 10) : obj.name
    }

    static String shortenRefName(String name) {
        if( name.startsWith('refs/remotes/origin/') )
            return name.replace('refs/remotes/origin/', '')

        return Repository.shortenRefName(name)
    }

    static Map refToMap(Ref ref, Map<String, Ref> remote) {
        final entry = new HashMap(2)
        final objId = ref.getPeeledObjectId() ?: ref.getObjectId()
        // the branch or tag name
        entry.name = shortenRefName(ref.name)
        // the local commit it
        entry.commitId = objId.name()
        // the remote commit Id for this branch or tag
        if( remote && hasRemoteChange(ref, remote) ) {
            entry.latestId = remote.get(ref.name).objectId.name()
        }
        return entry
    }

    static boolean isRefInCommits(Ref ref, List<String> commits) {
        if( !commits )
            return false
        String peeledId = ref.getPeeledObjectId()?.name()
        String id = ref.getObjectId()?.name()
        for( String commit : commits ) {
            if( commit.equals(peeledId) || commit.equals(id) )
                return true
        }
        return false
    }


}
