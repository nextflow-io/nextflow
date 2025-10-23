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

package io.seqera.tower.plugin.scm.jgit;

import org.eclipse.jgit.errors.TransportException;
import org.eclipse.jgit.lib.Repository;
import org.eclipse.jgit.transport.*;

import java.util.Collections;
import java.util.Set;

/**
 * JGit transport implementation for Seqera Platform data-links git-remote storage.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class TransportSeqera extends Transport {

    public static final TransportProtocol PROTO_SEQERA = new SeqeraTransportProtocol();

    public TransportSeqera(Repository local, URIish uri) throws TransportException {
        super(local, uri);
    }

    @Override
    public FetchConnection openFetch() throws TransportException {
        return new SeqeraFetchConnection(this);
    }

    @Override
    public PushConnection openPush() throws TransportException {
        return new SeqeraPushConnection(this);
    }

    @Override
    public void close() {
        // cleanup resources if needed
    }

    public Repository getLocal() {
        return this.local;
    }

    public static class SeqeraTransportProtocol extends TransportProtocol {
        @Override
        public String getName() {
            return "Seqera Platform";
        }

        @Override
        public Set<String> getSchemes() {
            return Collections.singleton("seqera");
        }

        @Override
        public boolean canHandle(URIish uri, Repository local, String remoteName) {
            return "seqera".equals(uri.getScheme());
        }

        @Override
        public Transport open(URIish uri, Repository local, String remoteName) throws TransportException {
            try {
                return new TransportSeqera(local, uri);
            } catch (TransportException e) {
                throw e;
            }
        }
    }

    public static void register() {
        Transport.register(PROTO_SEQERA);
    }
}
