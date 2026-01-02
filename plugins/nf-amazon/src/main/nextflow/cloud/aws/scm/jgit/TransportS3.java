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
import org.eclipse.jgit.lib.Repository;
import org.eclipse.jgit.transport.*;

import java.util.Collections;
import java.util.Set;

/**
 * JGit transport implementation compatible with git-remote-s3 storage.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class TransportS3 extends Transport {

    public static final TransportProtocol PROTO_S3 = new S3TransportProtocol();

    public TransportS3(Repository local, URIish uri) throws TransportException {
        super(local, uri);
    }

    public TransportS3(URIish uri) throws TransportException {
        super(uri);
    }

    @Override
    public FetchConnection openFetch() throws TransportException {
        return new S3FetchConnection(this);
    }

    @Override
    public PushConnection openPush() throws TransportException {
        return new S3PushConnection(this);
    }

    // Optional: Clean up if needed
    @Override
    public void close() {
        // cleanup resources if needed
    }

    public Repository getLocal(){
        return this.local;
    }

    public static class S3TransportProtocol extends TransportProtocol {
        @Override
        public String getName() {
            return "Amazon S3";
        }

        @Override
        public Set<String> getSchemes() {
            return Collections.singleton("s3");
        }

        @Override
        public boolean canHandle(URIish uri, Repository local, String remoteName) {
            return "s3".equals(uri.getScheme());
        }

        @Override
        public Transport open(URIish uri, Repository local, String remoteName) throws TransportException {
            try {
                return new TransportS3(local, uri);
            } catch (TransportException e) {
                throw e;
            }
        }

        @Override
        public Transport open(URIish uri) throws TransportException {
            try {
                return new TransportS3(uri);
            } catch (TransportException e) {
                throw e;
            }
        }
    }

    public static void register() {
        Transport.register(PROTO_S3);
    }
}

