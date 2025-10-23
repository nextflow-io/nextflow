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

import org.eclipse.jgit.errors.UnsupportedCredentialItem;
import org.eclipse.jgit.transport.CredentialItem;
import org.eclipse.jgit.transport.CredentialsProvider;
import org.eclipse.jgit.transport.URIish;

/**
 * JGit credentials provider wrapper for the Seqera Platform API credentials and configuration.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class SeqeraGitCredentialsProvider extends CredentialsProvider {

    private String endpoint;
    private String accessToken;
    private String workspaceId;

    public void setEndpoint(String endpoint) {
        this.endpoint = endpoint;
    }

    public void setAccessToken(String accessToken) {
        this.accessToken = accessToken;
    }

    public void setWorkspaceId(String workspaceId) {
        this.workspaceId = workspaceId;
    }

    public String getEndpoint() {
        return endpoint != null ? endpoint : "https://api.cloud.seqera.io";
    }

    public String getAccessToken() {
        if (accessToken == null) {
            throw new IllegalStateException("Seqera Platform access token not configured");
        }
        return accessToken;
    }

    public String getWorkspaceId() {
        return workspaceId;
    }

    @Override
    public boolean isInteractive() {
        return false;
    }

    @Override
    public boolean supports(CredentialItem... items) {
        return false;
    }

    @Override
    public boolean get(URIish uri, CredentialItem... items) throws UnsupportedCredentialItem {
        return false;
    }
}
