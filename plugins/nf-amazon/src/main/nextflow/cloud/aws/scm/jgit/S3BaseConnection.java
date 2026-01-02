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

import org.eclipse.jgit.lib.*;
import org.eclipse.jgit.transport.Connection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import software.amazon.awssdk.core.ResponseInputStream;
import software.amazon.awssdk.services.s3.S3Client;
import software.amazon.awssdk.services.s3.model.GetObjectRequest;
import software.amazon.awssdk.services.s3.model.ListObjectsV2Request;
import software.amazon.awssdk.services.s3.model.S3Object;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Base class for connections with a git-remote-s3 compatibility
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 *
 */
public class S3BaseConnection implements Connection {

    protected static final Logger log = LoggerFactory.getLogger(S3BaseConnection.class);
    protected TransportS3 transport;
    protected String bucket;
    protected String key;
    protected S3Client s3;
    private final Map<String, Ref> advertisedRefs = new HashMap<String, Ref>();

    public S3BaseConnection(TransportS3 transport) {
        this.transport = transport;
        this.bucket = transport.getURI().getHost();
        this.key = transport.getURI().getPath().substring(1);
        S3GitCredentialsProvider credentials = (S3GitCredentialsProvider) transport.getCredentialsProvider();
        this.s3 = S3Client.builder()
                .region(credentials.getRegion())
                .credentialsProvider(credentials.getAwsCredentialsProvider())
                .build();
        loadRefsMap();
        log.trace("Created S3 Connection for  s3://{}/{}", bucket, key);
    }

    protected String getDefaultBranchRef() throws IOException {

        GetObjectRequest getObjectRequest = GetObjectRequest.builder()
            .bucket(bucket)
            .key(key + "/HEAD")
            .build();

        try (ResponseInputStream<?> inputStream = s3.getObject(getObjectRequest)) {
            return new String(inputStream.readAllBytes(), StandardCharsets.UTF_8);
        }
    }



    private void loadRefsMap() {
        try {
            addRefs("heads");
        } catch (Exception e){
            log.debug("No heads found for s3://{}/{}: {}", bucket, key, e.getMessage());
        }
        try {
            addRefs("tags");
        } catch (Exception e){
            log.debug("No tags found for s3://{}/{}: {}", bucket, key, e.getMessage());
        }
        try {
            final String defaultBranch = getDefaultBranchRef();
            Ref target = advertisedRefs.get(defaultBranch);
            if (target != null)
                advertisedRefs.put(Constants.HEAD, new SymbolicRef(Constants.HEAD, target));
        }   catch (Exception e){
            log.debug("No default refs found for s3://{}/{}: {}", bucket, key, e.getMessage());
        }



    }

    private void addRefs(String refType) {
        final List<S3Object> list = s3.listObjectsV2(ListObjectsV2Request.builder()
            .bucket(bucket)
            .prefix(key + "/refs/" + refType)
            .build()
        ).contents();

        if (list == null || list.isEmpty()) {
            log.debug("No {} refs found for s3://{}/{}", refType, bucket, key);
            return;
        }

        for (S3Object obj : list) {
            addRef(obj);
        }
    }
    private void addRef(S3Object obj) {
        String key = obj.key();
        BranchData branch = BranchData.fromKey(key);
        String type = branch.getType();
        if ("heads".equals(type)) {
            advertisedRefs.put(branch.getRefName(), new ObjectIdRef.PeeledNonTag(Ref.Storage.NETWORK, branch.getRefName(), branch.getObjectId()));
        } else if ("tags".equals(type)) {
            advertisedRefs.put(branch.getRefName(), new ObjectIdRef.Unpeeled(Ref.Storage.NETWORK, branch.getRefName(), branch.getObjectId()));
        }
    }

    @Override
    public Map<String, Ref> getRefsMap() {
        return advertisedRefs;
    }

    @Override
    public Collection<Ref> getRefs() {
        return advertisedRefs.values();
    }

    @Override
    public Ref getRef(String name) {
        return advertisedRefs.get(name);
    }

    @Override
    public void close() { }

    @Override
    public String getMessages() {
        return "";
    }

    @Override
    public String getPeerUserAgent() {
        return "";
    }

    static class BranchData{
        private String type;
        private String simpleName;
        private String refName;
        private ObjectId objectId;

        private BranchData(String type, String simpleName, ObjectId objectId){
            this.type = type;
            this.simpleName = simpleName;
            this.refName = String.format("refs/%s/%s" , type, simpleName);
            this.objectId = objectId;
        }
        public static BranchData fromKey(String key){
            String[] parts = key.split("/");
            if (parts.length < 5) throw new RuntimeException("Incorrect s3 key parts inside the S3-git-remote repository. Key should include the following parts: repo-path/refs/<type>/<branch>/<hash>.bundle");
            final String type = parts[parts.length - 3];
            final String rBranch = parts[parts.length - 2];
            final String sha = parts[parts.length - 1].replace(".bundle", "");
            return new BranchData(type, rBranch, ObjectId.fromString(sha));

        }

        public String getType() {
            return type;
        }

        public String getSimpleName() {
            return simpleName;
        }

        public String getRefName() {
            return refName;
        }

        public ObjectId getObjectId() {
            return objectId;
        }
    }

}