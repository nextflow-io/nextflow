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

package nextflow.cloud.google.util

import com.google.auth.oauth2.GoogleCredentials
import nextflow.Session
import nextflow.cloud.google.GoogleOpts
import spock.lang.Specification

class GsStorageOptionsTest extends Specification {

    /**
     * {@link GoogleOpts} with credential resolution disabled. Nothing in this spec is about
     * credentials -- the subject is the timeouts and the retry policy that used to be dropped -- but
     * resolving them means Application Default Credentials, which are AMBIENT: present on a
     * developer machine with `gcloud` configured, absent on a CI runner, where the lookup throws
     * {@code IOException} before any assertion is reached.
     */
    static class NoAdcOpts extends GoogleOpts {
        NoAdcOpts(Map opts) { super(opts) }
        @Override GoogleCredentials getCredentials() { return null }
    }

    def 'applies the configured timeouts and retry policy, not just credentials and project'() {
        given: 'the cache client used to take credentials + project only, dropping these'
        def opts = new NoAdcOpts([
                project: 'my-proj',
                httpConnectTimeout: '11s',
                httpReadTimeout: '22s',
                storage: [retryPolicy: [maxAttempts: 7, multiplier: 3.0d]] ])

        when:
        def built = GsStorageOptions.optionsFor(opts)

        then: 'project id is applied'
        built.getProjectId() == 'my-proj'
        and: 'and so are the transport timeouts the user configured'
        built.getTransportOptions().getConnectTimeout() == 11_000
        built.getTransportOptions().getReadTimeout() == 22_000
        and: 'and the retry policy'
        built.getRetrySettings().getMaxAttempts() == 7
        built.getRetrySettings().getRetryDelayMultiplier() == 3.0d
    }

    private GoogleOpts optsOf(Map google) {
        return GoogleOpts.fromSession(Stub(Session) { getConfig() >> [google: google] })
    }

    def 'userProject is set only for a requester-pays bucket'() {
        expect: 'off by default -- nothing to bill'
        GsStorageOptions.userProject(optsOf([project: 'p'])) == null

        and: 'on when requester-pays is enabled, since the cache must pass it PER REQUEST:'
        and: 'StorageOptions has no userProject setter, so a client-level default is impossible'
        GsStorageOptions.userProject(optsOf([project: 'p', enableRequesterPaysBuckets: true])) == 'p'

        and: 'and null config does not blow up -- the providers resolve it without a session in tests'
        GsStorageOptions.userProject(null) == null
    }
}
