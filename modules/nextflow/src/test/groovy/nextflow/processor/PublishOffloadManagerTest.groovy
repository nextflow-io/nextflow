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
 */

package nextflow.processor

import nextflow.Session
import nextflow.script.BaseScript
import nextflow.script.ScriptBinding
import nextflow.script.ScriptFile
import spock.lang.Specification

import test.TestHelper


class PublishOffloadManagerTest extends Specification {

    def 'should create task processor'(){
        given:
        def session = new Session();
        def scriptBinding = new ScriptBinding(session: session)
        def script = Stub(BaseScript)
        script.getBinding() >> scriptBinding
        def folder = TestHelper.createInMemTempDir()
        def file = folder.resolve('pipeline.nf'); file.text = 'println "hello"'
        def scriptFile = new ScriptFile(file)
        session.init(scriptFile)
        //session.start()
        session.script = script;
        def poManager = new PublishOffloadManager(session);
        when:
        poManager.init()
        then:
        poManager.publishProcessor != null
        cleanup:
        session.classesDir?.deleteDir()


    }

}
