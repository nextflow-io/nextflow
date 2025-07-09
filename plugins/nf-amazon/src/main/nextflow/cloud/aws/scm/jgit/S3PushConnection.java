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
import org.eclipse.jgit.lib.NullProgressMonitor;
import org.eclipse.jgit.lib.ProgressMonitor;
import org.eclipse.jgit.lib.Ref;
import org.eclipse.jgit.transport.BundleWriter;
import org.eclipse.jgit.transport.PushConnection;
import org.eclipse.jgit.transport.RemoteRefUpdate;
import org.eclipse.jgit.util.FileUtils;
import software.amazon.awssdk.services.s3.model.PutObjectRequest;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
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
                log.trace("Generating bundle for {} in {} ", entry.getKey(), tmpdir);
                Path bundleFile = bundle(entry.getKey(), tmpdir);
                s3.putObject(PutObjectRequest.builder().bucket(bucket).key(key+'/'+entry.getKey()).build(),bundleFile);
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

    @Override
    public void push(ProgressMonitor monitor, Map<String, RemoteRefUpdate> refUpdates, OutputStream out) throws TransportException {
        push(monitor, refUpdates);

    }

    private Path bundle(String refName, Path tmpdir) throws IOException {
        final Ref ref = transport.getLocal().findRef(refName);
        if( ref == null ) {
            throw new IllegalStateException("Branch ${branch} not found");
        }
        final BundleWriter writer = new BundleWriter(transport.getLocal());
        Path bundleFile = tmpdir.resolve(ref.getObjectId() +".bundle");
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