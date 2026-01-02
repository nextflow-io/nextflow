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

import org.eclipse.jgit.errors.UnsupportedCredentialItem;
import org.eclipse.jgit.transport.CredentialItem;
import org.eclipse.jgit.transport.CredentialsProvider;
import org.eclipse.jgit.transport.URIish;
import software.amazon.awssdk.auth.credentials.AwsCredentialsProvider;
import software.amazon.awssdk.auth.credentials.DefaultCredentialsProvider;
import software.amazon.awssdk.regions.Region;
import software.amazon.awssdk.regions.providers.DefaultAwsRegionProviderChain;

/**
 * JGit credentials provider wrapper for the AWS credentialsProvider and other client configuration parameters.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class S3GitCredentialsProvider extends CredentialsProvider {

    private Region region;
    private AwsCredentialsProvider awsCredentialsProvider;

    public void setAwsCredentialsProvider(AwsCredentialsProvider provider){
        this.awsCredentialsProvider = provider;
    }
    public void setRegion(Region region) {
        this.region = region;
    }

    public Region getRegion(){
        if( region != null )
            return region;
        var region = DefaultAwsRegionProviderChain.builder().build().getRegion();
        if( region != null)
            return region;
        return Region.US_EAST_1;
    }

    public AwsCredentialsProvider getAwsCredentialsProvider(){
        return awsCredentialsProvider != null ? awsCredentialsProvider : DefaultCredentialsProvider.builder().build();
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
