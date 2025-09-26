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

package nextflow.cloud.aws.scm.jgit;

import org.eclipse.jgit.errors.TransportException;
import org.eclipse.jgit.lib.*;
import org.eclipse.jgit.revwalk.RevCommit;
import org.eclipse.jgit.revwalk.RevWalk;
import org.eclipse.jgit.transport.BundleWriter;
import org.eclipse.jgit.transport.PushConnection;
import org.eclipse.jgit.transport.RemoteRefUpdate;
import org.eclipse.jgit.util.FileUtils;
import software.amazon.awssdk.core.sync.RequestBody;
import software.amazon.awssdk.services.s3.model.*;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

/**
 * Push connection implementation compatible with git-remote-s3 storage.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class S3PushConnection extends S3BaseConnection implements PushConnection {

    public S3PushConnection(TransportS3 transport) {
        super(transport);
    }

    @Override
    public void push(ProgressMonitor monitor, Map<String, RemoteRefUpdate> refUpdates) throws TransportException {
        Path tmpdir = null;
        try {
            tmpdir = Files.createTempDirectory("s3-remote-git-");
            for (Map.Entry<String, RemoteRefUpdate> entry : refUpdates.entrySet()) {
                pushBranch(entry, tmpdir);
            }
        }catch (IOException e){
            throw new TransportException(transport.getURI(), "Exception fetching branches", e);
        }finally {
            if (tmpdir != null)
                try {
                    FileUtils.delete(tmpdir.toFile(), FileUtils.RECURSIVE);
                }catch (IOException e){
                    throw new TransportException(transport.getURI(), "Exception fetching branches", e);
                }
        }

    }

    private void pushBranch(Map.Entry<String, RemoteRefUpdate> entry, Path tmpdir) throws IOException {
        log.debug("Pushing {} reference", entry.getKey());
        final Ref ref = transport.getLocal().findRef(entry.getKey());
        if( ref == null || ref.getObjectId() == null) {
            throw new IllegalStateException("Branch ${branch} not found");
        }
        S3Object oldObject = checkExistingObjectInBranch(entry.getKey());
        if( oldObject != null && isSameObjectId(oldObject, ref.getObjectId())){
            setUpdateStatus(entry.getValue(), RemoteRefUpdate.Status.UP_TO_DATE);
            return;
        }
        if( oldObject != null && !isCommitInBranch(oldObject, ref)) {
            setUpdateStatus(entry.getValue(), RemoteRefUpdate.Status.REJECTED_REMOTE_CHANGED);
            return;
        }
        log.trace("Generating bundle for branch {} in {}", entry.getKey(), tmpdir);
        Path bundleFile = bundle(ref, tmpdir);
        String objectKey = String.format("%s/%s/%s", key, entry.getKey(), bundleFile.getFileName().toString());

        log.trace("Uploading bundle {} to s3://{}/{}", bundleFile, bucket, objectKey);
        s3.putObject(PutObjectRequest.builder()
                .bucket(bucket)
                .key(objectKey)
                .build(),
            bundleFile);
        if( oldObject != null ){
            log.trace("Deleting old bundle s3://{}/{}",bucket,oldObject.key());
            s3.deleteObject(DeleteObjectRequest.builder()
                .bucket(bucket)
                .key(oldObject.key())
                .build()
            );
        }
        setUpdateStatus(entry.getValue(), RemoteRefUpdate.Status.OK);
        if ( getRef(Constants.HEAD) == null) {
            updateRemoteHead(entry.getKey());
        }
    }

    private void updateRemoteHead(String ref) {
        try {
            s3.headObject(HeadObjectRequest.builder()
                .bucket(bucket)
                .key(key + "/HEAD")
                .build());
        } catch (NoSuchKeyException e) {
            log.debug("No remote default branch. Setting to {}.", ref);
            s3.putObject(PutObjectRequest.builder()
                    .bucket(bucket)
                    .key(key + "/HEAD")
                    .build(),
                RequestBody.fromBytes(ref.getBytes(StandardCharsets.UTF_8)));
        }
    }

    private boolean isSameObjectId(S3Object s3object, ObjectId commitId){
        return BranchData.fromKey(s3object.key()).getObjectId().name().equals(commitId.name());
    }

   private void setUpdateStatus(RemoteRefUpdate update, RemoteRefUpdate.Status status) {
        try {
            Field statusField = RemoteRefUpdate.class.getDeclaredField("status");
            statusField.setAccessible(true);
            statusField.set(update, status);
        } catch (Exception e) {
            throw new RuntimeException("Unable to set status on RemoteRefUpdate", e);
        }
    }

    public boolean isCommitInBranch(S3Object s3Object, Ref branchRef) throws IOException {
        ObjectId commitId = BranchData.fromKey(s3Object.key()).getObjectId();
        try (RevWalk walk = new RevWalk(transport.getLocal())) {
            RevCommit branchTip = walk.parseCommit(branchRef.getObjectId());
            RevCommit targetCommit = walk.parseCommit(commitId);

            // Check if the commit is reachable from the branch tip
            return walk.isMergedInto(targetCommit, branchTip);
        }
    }

    private S3Object checkExistingObjectInBranch(String name) throws TransportException {
        final List<S3Object> list = s3.listObjectsV2(ListObjectsV2Request.builder()
                    .bucket(bucket)
                    .prefix(key + '/' + name + '/')
                    .build()
            ).contents();

        if( list == null || list.isEmpty() ) {
            return null;
        }

        if( list.size() > 1 ){
            throw new TransportException(transport.getURI(), " More than one bundle for " + name);
        }
        return list.get(0);
    }

    @Override
    public void push(ProgressMonitor monitor, Map<String, RemoteRefUpdate> refUpdates, OutputStream out) throws TransportException {
        push(monitor, refUpdates);

    }

    private Path bundle(Ref ref, Path tmpdir) throws IOException {
        final BundleWriter writer = new BundleWriter(transport.getLocal());
        Path bundleFile = tmpdir.resolve(ref.getObjectId().name() +".bundle");
        writer.include(ref);
        try (OutputStream out = new FileOutputStream(bundleFile.toFile())) {
            writer.writeBundle(NullProgressMonitor.INSTANCE, out);
        }
        return bundleFile;
    }

    @Override
    public void close() {

    }
}