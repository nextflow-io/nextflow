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

import org.eclipse.jgit.lib.ObjectId
import org.eclipse.jgit.lib.Ref
import org.junit.Rule
import spock.lang.Requires
import spock.lang.Specification
import test.TemporaryPath

class GitReferenceHelperTest extends Specification {

    @Rule
    TemporaryPath tempDir = new TemporaryPath()

    def setup() {
        AssetManager.root = tempDir.root.toFile()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'can filter remote branches'() {
        given:
        def folder = tempDir.getRoot()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])
        manager.download()
        def branches = manager.getBranchList()

        when:
        def remote_head = branches.find { it.name == 'refs/remotes/origin/HEAD' }
        then:
        remote_head != null
        !GitReferenceHelper.isRemoteBranch(remote_head)

        when:
        def remote_master = branches.find { it.name == 'refs/remotes/origin/master' }
        then:
        remote_master != null
        GitReferenceHelper.isRemoteBranch(remote_master)

        when:
        def local_master = branches.find { it.name == 'refs/heads/master' }
        then:
        local_master != null
        !GitReferenceHelper.isRemoteBranch(local_master)
    }

    def 'should check hasRemoteChange with no remote'() {
        given:
        def ref = Mock(Ref)

        expect:
        !GitReferenceHelper.hasRemoteChange(ref, null)
        !GitReferenceHelper.hasRemoteChange(ref, [:])
    }

    def 'should check hasRemoteChange when ref not in remote'() {
        given:
        def ref = Mock(Ref) {
            getName() >> 'refs/heads/master'
        }
        def remote = ['refs/heads/other': Mock(Ref)]

        expect:
        !GitReferenceHelper.hasRemoteChange(ref, remote)
    }

    def 'should detect hasRemoteChange when ObjectIds differ'() {
        given:
        def localId = ObjectId.fromString('1234567890123456789012345678901234567890')
        def remoteId = ObjectId.fromString('abcdefabcdefabcdefabcdefabcdefabcdefabcd')

        def localRef = Mock(Ref) {
            getName() >> 'refs/heads/master'
            getObjectId() >> localId
        }
        def remoteRef = Mock(Ref) {
            getObjectId() >> remoteId
        }
        def remote = ['refs/heads/master': remoteRef]

        expect:
        GitReferenceHelper.hasRemoteChange(localRef, remote)
    }

    def 'should not detect hasRemoteChange when ObjectIds are same'() {
        given:
        def sameId = ObjectId.fromString('1234567890123456789012345678901234567890')

        def localRef = Mock(Ref) {
            getName() >> 'refs/heads/master'
            getObjectId() >> sameId
        }
        def remoteRef = Mock(Ref) {
            getObjectId() >> sameId
        }
        def remote = ['refs/heads/master': remoteRef]

        expect:
        !GitReferenceHelper.hasRemoteChange(localRef, remote)
    }

    def 'should format ObjectId in human readable form'() {
        given:
        def objectId = ObjectId.fromString('1234567890123456789012345678901234567890')

        when:
        def result = GitReferenceHelper.formatObjectId(objectId, true)

        then:
        result == '1234567890'
    }

    def 'should format ObjectId in full form'() {
        given:
        def objectId = ObjectId.fromString('1234567890123456789012345678901234567890')

        when:
        def result = GitReferenceHelper.formatObjectId(objectId, false)

        then:
        result == '1234567890123456789012345678901234567890'
    }

    def 'should format update with level 0'() {
        given:
        def objectId = ObjectId.fromString('1234567890123456789012345678901234567890')
        def remoteRef = Mock(Ref) {
            getObjectId() >> objectId
        }

        when:
        def result = GitReferenceHelper.formatUpdate(remoteRef, 0)

        then:
        result == 'updates on remote'
    }

    def 'should format update with level 1'() {
        given:
        def objectId = ObjectId.fromString('1234567890123456789012345678901234567890')
        def remoteRef = Mock(Ref) {
            getObjectId() >> objectId
        }

        when:
        def result = GitReferenceHelper.formatUpdate(remoteRef, 1)

        then:
        result == 'updates on remote 1234567890'
    }

    def 'should format update with level 2'() {
        given:
        def objectId = ObjectId.fromString('1234567890123456789012345678901234567890')
        def remoteRef = Mock(Ref) {
            getObjectId() >> objectId
        }

        when:
        def result = GitReferenceHelper.formatUpdate(remoteRef, 2)

        then:
        result == 'updates on remote 1234567890123456789012345678901234567890'
    }

    def 'should shorten remote ref name'() {
        expect:
        GitReferenceHelper.shortenRefName('refs/remotes/origin/master') == 'master'
        GitReferenceHelper.shortenRefName('refs/remotes/origin/develop') == 'develop'
        GitReferenceHelper.shortenRefName('refs/remotes/origin/feature/test') == 'feature/test'
    }

    def 'should shorten standard ref names'() {
        expect:
        GitReferenceHelper.shortenRefName('refs/heads/master') == 'master'
        GitReferenceHelper.shortenRefName('refs/tags/v1.0.0') == 'v1.0.0'
    }

    def 'should not shorten already short ref names'() {
        expect:
        GitReferenceHelper.shortenRefName('master') == 'master'
        GitReferenceHelper.shortenRefName('develop') == 'develop'
    }

    def 'should convert ref to map without remote'() {
        given:
        def objectId = ObjectId.fromString('1234567890123456789012345678901234567890')
        def ref = Mock(Ref) {
            getName() >> 'refs/heads/master'
            getObjectId() >> objectId
            getPeeledObjectId() >> null
        }

        when:
        def result = GitReferenceHelper.refToMap(ref, null)

        then:
        result.name == 'master'
        result.commitId == '1234567890123456789012345678901234567890'
        !result.containsKey('latestId')
    }

    def 'should convert ref to map with peeled ObjectId'() {
        given:
        def objectId = ObjectId.fromString('1234567890123456789012345678901234567890')
        def peeledId = ObjectId.fromString('abcdefabcdefabcdefabcdefabcdefabcdefabcd')
        def ref = Mock(Ref) {
            getName() >> 'refs/tags/v1.0.0'
            getObjectId() >> objectId
            getPeeledObjectId() >> peeledId
        }

        when:
        def result = GitReferenceHelper.refToMap(ref, null)

        then:
        result.name == 'v1.0.0'
        result.commitId == 'abcdefabcdefabcdefabcdefabcdefabcdefabcd'
        !result.containsKey('latestId')
    }

    def 'should convert ref to map with remote changes'() {
        given:
        def localId = ObjectId.fromString('1234567890123456789012345678901234567890')
        def remoteId = ObjectId.fromString('abcdefabcdefabcdefabcdefabcdefabcdefabcd')

        def ref = Mock(Ref) {
            getName() >> 'refs/heads/master'
            getObjectId() >> localId
            getPeeledObjectId() >> null
        }
        def remoteRef = Mock(Ref) {
            getObjectId() >> remoteId
        }
        def remote = ['refs/heads/master': remoteRef]

        when:
        def result = GitReferenceHelper.refToMap(ref, remote)

        then:
        result.name == 'master'
        result.commitId == '1234567890123456789012345678901234567890'
        result.latestId == 'abcdefabcdefabcdefabcdefabcdefabcdefabcd'
    }

    def 'should convert ref to map without remote changes'() {
        given:
        def sameId = ObjectId.fromString('1234567890123456789012345678901234567890')

        def ref = Mock(Ref) {
            getName() >> 'refs/heads/master'
            getObjectId() >> sameId
            getPeeledObjectId() >> null
        }
        def remoteRef = Mock(Ref) {
            getObjectId() >> sameId
        }
        def remote = ['refs/heads/master': remoteRef]

        when:
        def result = GitReferenceHelper.refToMap(ref, remote)

        then:
        result.name == 'master'
        result.commitId == '1234567890123456789012345678901234567890'
        !result.containsKey('latestId')
    }

    def 'should find ref in commits using ObjectId'() {
        given:
        def objectId = ObjectId.fromString('1234567890123456789012345678901234567890')
        def ref = Mock(Ref) {
            getObjectId() >> objectId
            getPeeledObjectId() >> null
        }
        def commits = ['abcdefabcdefabcdefabcdefabcdefabcdefabcd', '1234567890123456789012345678901234567890']

        expect:
        GitReferenceHelper.isRefInCommits(ref, commits)
    }

    def 'should find ref in commits using peeled ObjectId'() {
        given:
        def objectId = ObjectId.fromString('1234567890123456789012345678901234567890')
        def peeledId = ObjectId.fromString('abcdefabcdefabcdefabcdefabcdefabcdefabcd')
        def ref = Mock(Ref) {
            getObjectId() >> objectId
            getPeeledObjectId() >> peeledId
        }
        def commits = ['abcdefabcdefabcdefabcdefabcdefabcdefabcd', 'fedcbafedcbafedcbafedcbafedcbafedcbafed']

        expect:
        GitReferenceHelper.isRefInCommits(ref, commits)
    }

    def 'should not find ref in commits when not present'() {
        given:
        def objectId = ObjectId.fromString('1234567890123456789012345678901234567890')
        def ref = Mock(Ref) {
            getObjectId() >> objectId
            getPeeledObjectId() >> null
        }
        def commits = ['abcdefabcdefabcdefabcdefabcdefabcdefabcd', 'fedcbafedcbafedcbafedcbafedcbafedcbafed']

        expect:
        !GitReferenceHelper.isRefInCommits(ref, commits)
    }
}
