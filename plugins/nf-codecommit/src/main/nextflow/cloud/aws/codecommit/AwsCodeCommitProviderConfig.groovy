/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.cloud.aws.codecommit

import groovy.transform.CompileStatic
import nextflow.scm.ProviderConfig
/**
 * A {@link ProviderConfig} specialised for AWS CodeCommit
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AwsCodeCommitProviderConfig extends ProviderConfig {

    AwsCodeCommitProviderConfig(String host) {
        super('codecommit', [platform:'codecommit', server: "https://$host"])
        assert host =~ /git-codecommit\.[a-z0-9-]+\.amazonaws\.com/, "Invalid AWS CodeCommit hostname: '$host'"
    }

    AwsCodeCommitProviderConfig(String name, Map attributes) {
        super(name, attributes)
        assert attributes.platform=='codecommit', "Invalid AWS CodeCommit platform value: '$attributes.platform'"
    }

    String getRegion() {
        final host = getDomain()
        if( !host )
            throw new IllegalStateException("Missing AWS CodeCommit repository name")
        final result = host.tokenize('.')[1]
        if( !result )
            throw new IllegalStateException("Invalid AWS CodeCommit hostname: '${host}'")
        return result
    }

    @Override
    protected String resolveProjectName(String path) {
        assert path
        assert !path.startsWith('/')
        final repoName = path.tokenize('/')[-1]
        if( !repoName )
            throw new IllegalArgumentException("Invalid AWS CodeCommit repository path: $path")
        return "codecommit-$region/$repoName"
    }

    String toString() {
        return "AwsCodeCommitProviderConfig[name=$name; platform=$platform; server=$server]"
    }

}
