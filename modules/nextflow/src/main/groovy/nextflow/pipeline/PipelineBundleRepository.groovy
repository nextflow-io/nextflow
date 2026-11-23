/*
 * Copyright 2013-2026, Seqera Labs
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


package nextflow.pipeline

import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.npr.client.RegistryClient
import nextflow.file.FileHelper
import nextflow.file.FileMutex
import nextflow.util.VersionNumber
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.lib.CommitBuilder
import org.eclipse.jgit.lib.Constants
import org.eclipse.jgit.lib.FileMode
import org.eclipse.jgit.lib.ObjectId
import org.eclipse.jgit.lib.ObjectInserter
import org.eclipse.jgit.lib.PersonIdent
import org.eclipse.jgit.lib.Ref
import org.eclipse.jgit.lib.Repository
import org.eclipse.jgit.lib.TreeFormatter

/** Bare Git repository with one tagged commit per bundle release, so the run cache can clone it */
@Slf4j
@CompileStatic
class PipelineBundleRepository {

    static final String BRANCH = 'master'

    // all releases are chained on this branch, so that a fetch of any version gets every tag
    static final String RELEASES = 'releases'

    // fixed so that the same bundle always gives the same commit
    private static final PersonIdent IDENT = new PersonIdent('Nextflow', 'nextflow@localhost', new Date(0), TimeZone.getTimeZone('UTC'))

    /** Add the bundle of a release as a commit tagged with its version, unless that tag exists */
    static void add(Path repo, RegistryClient client, PipelineInfo release) {
        Files.createDirectories(repo.parent)
        final lock = repo.resolveSibling(".${repo.fileName}.lock").toFile()
        final mutex = new FileMutex(target: lock, timeout: '60s', waitMessage: "Another Nextflow instance is adding a release to ${repo} -- please wait", errorMessage: "Unable to acquire exclusive lock after 60s on file: ${lock}")
        try {
            mutex.lock { add0(repo, client, release) }
        }
        finally {
            lock.delete()
        }
    }

    private static void add0(Path repo, RegistryClient client, PipelineInfo release) {
        try( final git = open(repo) ) {
            final db = git.repository
            if( db.findRef(Constants.R_TAGS + release.version) != null )
                return
            final tmp = Files.createTempDirectory('nf-bundle-')
            try {
                final bundle = tmp.resolve('bundle.tar.gz')
                PipelineManager.downloadBundle(client, release, bundle)
                final files = Files.createDirectory(tmp.resolve('files'))
                PipelineBundle.extract(bundle, files)
                final parent = db.resolve(Constants.R_HEADS + RELEASES)
                final commit = db.parseCommit(commit(db, files, parent, "${release.name}@${release.version}\n\nartifact: ${release.artifactChecksum}\n"))
                git.tag().setName(release.version).setAnnotated(false).setObjectId(commit).call()
                git.branchCreate().setName(RELEASES).setForce(true).setStartPoint(commit).call()
                git.branchCreate().setName(BRANCH).setForce(true).setStartPoint(db.parseCommit(latest(db))).call()
            }
            finally {
                FileHelper.deletePath(tmp)
            }
        }
    }

    private static Git open(Path repo) {
        return Files.exists(repo.resolve('HEAD'))
            ? Git.open(repo.toFile())
            : Git.init().setBare(true).setInitialBranch(BRANCH).setDirectory(repo.toFile()).call()
    }

    private static ObjectId commit(Repository db, Path dir, ObjectId parent, String message) {
        try( final inserter = db.newObjectInserter() ) {
            final builder = new CommitBuilder()
            builder.treeId = tree(inserter, dir) ?: inserter.insert(new TreeFormatter())
            if( parent )
                builder.setParentId(parent)
            builder.author = IDENT
            builder.committer = IDENT
            builder.message = message
            final id = inserter.insert(builder)
            inserter.flush()
            return id
        }
    }

    // not `git add`: it would skip files that the bundle's own .gitignore matches
    private static ObjectId tree(ObjectInserter inserter, Path dir) {
        final formatter = new TreeFormatter()
        boolean empty = true
        // Git sorts tree entries by name, comparing a folder as its name followed by '/'
        final entries = Files.list(dir).withCloseable { it.toList() }.sort { Path it -> Files.isDirectory(it, LinkOption.NOFOLLOW_LINKS) ? "${it.fileName}/".toString() : it.fileName.toString() }
        for( Path entry : entries ) {
            final name = entry.fileName.toString()
            if( Files.isDirectory(entry, LinkOption.NOFOLLOW_LINKS) ) {
                final sub = tree(inserter, entry)
                if( sub ) {
                    formatter.append(name, FileMode.TREE, sub)
                    empty = false
                }
            }
            else if( Files.isRegularFile(entry, LinkOption.NOFOLLOW_LINKS) ) {
                final blob = Files.newInputStream(entry).withCloseable { InputStream it -> inserter.insert(Constants.OBJ_BLOB, Files.size(entry), it) }
                formatter.append(name, Files.isExecutable(entry) ? FileMode.EXECUTABLE_FILE : FileMode.REGULAR_FILE, blob)
                empty = false
            }
        }
        return empty ? null : inserter.insert(formatter)
    }

    /** The commit of the highest version, so that the default branch is the latest release */
    private static ObjectId latest(Repository db) {
        final top = db.refDatabase.getRefsByPrefix(Constants.R_TAGS).max { Ref a, Ref b -> compare(Repository.shortenRefName(a.name), Repository.shortenRefName(b.name)) }
        return db.refDatabase.peel(top).peeledObjectId ?: top.objectId
    }

    /** Semver order: a pre-release such as 2.0.0-rc1 comes before 2.0.0 */
    static int compare(String a, String b) {
        final pa = a.split('-', 2)
        final pb = b.split('-', 2)
        final core = new VersionNumber(pa[0]) <=> new VersionNumber(pb[0])
        if( core != 0 )
            return core
        if( pa.length != pb.length )
            return pa.length < pb.length ? 1 : -1
        return pa.length == 1 ? 0 : new VersionNumber(pa[1]) <=> new VersionNumber(pb[1])
    }
}
