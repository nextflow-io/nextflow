/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
package nextflow.cloud.google.lifesciences

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import nextflow.Session
import nextflow.cloud.google.GoogleSpecification
import nextflow.exception.AbortOperationException
import spock.lang.Shared

class GoogleLifeSciencesExecutorTest extends GoogleSpecification {

    @Shared Path CREDS_FILE
    @Shared Path TEMP_DIR

    def setup() {
        TEMP_DIR = Files.createTempDirectory('test')
        CREDS_FILE = TEMP_DIR.resolve('google.json')
        CREDS_FILE.text = '''
        {
            "project_id": "my-project-123"
        }
        '''
    }

    def cleanup() {
        TEMP_DIR?.deleteDir()
    }

    Map<String,String> ENV() {
        [GOOGLE_APPLICATION_CREDENTIALS: CREDS_FILE.toString()]
    }

    def 'should config and init executor'() {
        given:
        def session = Stub(Session)
        def path = mockGsPath('gs://foo/bar')
        session.bucketDir >> path
        session.binDir >> null
        session.config >> [google: [project: 'my-project-123', region: 'eu-west1', location:'us-central']]
        and:
        def executor = Spy(GoogleLifeSciencesExecutor)
        executor.session = session
        executor.env = ENV()

        when:
        executor.register()
        then:
        1 * executor.initClient() >> Mock(GoogleLifeSciencesHelper)
        and:
        executor.config.location == 'us-central'
        executor.config.regions == ['eu-west1']
        executor.config.project == 'my-project-123'
        and:
        executor.helper != null
    }

    def 'should get project from creds file'() {
        given:
        def session = Stub(Session)
        def path = mockGsPath('gs://foo/bar')
        session.bucketDir >> path
        session.binDir >> null
        session.config >> [google: [region: 'eu-west1', location:'us-central']]
        and:
        def executor = Spy(GoogleLifeSciencesExecutor)
        executor.session = session
        executor.env = ENV()

        when:
        executor.register()
        then:
        1 * executor.initClient() >> Mock(GoogleLifeSciencesHelper)
        and:
        executor.config.location == 'us-central'
        executor.config.regions == ['eu-west1']
        executor.config.project == 'my-project-123'
        and:
        executor.helper != null
    }

    def 'should fail because project does not match'() {
        given:
        def session = Stub(Session)
        def path = mockGsPath('gs://foo/bar')
        session.bucketDir >> path
        session.binDir >> null
        session.config >> [google: [project: 'another-project', region: 'eu-west1', location:'us-central']]
        and:
        def executor = Spy(GoogleLifeSciencesExecutor)
        executor.session = session
        executor.env = ENV()

        when:
        executor.register()
        then:
        def err = thrown(AbortOperationException)
        and:
        err.message.startsWith('Project Id `another-project` declared in the nextflow config file')
    }

    def 'should abort operation when the workdir is not a CloudStoragePath'() {
        given:
        def session = Stub(Session)
        session.workDir = Stub(Path)
        and:
        def executor = new GoogleLifeSciencesExecutor(env: ENV(), session: session)

        when:
        executor.register()

        then:
        def error = thrown(AbortOperationException)
        error.getMessage().startsWith("Executor `google-lifesciences` requires a Google Storage bucket")
    }


    def 'should stop on missing credentials' () {
        given:
        def session = Mock(Session)
        def path = mockGsPath('gs://foo/bar')
        session.bucketDir >> path
        session.binDir >> null
        and:
        def executor = new GoogleLifeSciencesExecutor(env: [:], session: session)

        when:
        executor.register()
        then:
        def err = thrown(AbortOperationException)
        err.message.startsWith('Missing Google credentials')
    }

    def 'should stop on missing bucket' () {
        given:
        def session = Mock(Session)
        def path = Paths.get('/local/dir')
        session.bucketDir >> path
        session.binDir >> null
        and:
        def executor = new GoogleLifeSciencesExecutor(env: ENV(), session: session)

        when:
        executor.register()
        then:
        def err = thrown(AbortOperationException)
        err.message.startsWith('Executor `google-lifesciences` requires a Google Storage bucket to be specified as a working directory')
    }

    def 'should be containerNative'() {
        when:
        def executor = new GoogleLifeSciencesExecutor()
        then:
        executor.isContainerNative()
    }
}
