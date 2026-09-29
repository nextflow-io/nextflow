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
package nextflow.agent.pi

import spock.lang.Specification

/**
 * Pins how the runner declares the image it needs: {@code getDefaultContainer()} is read out of
 * {@code META-INF/nf-agent-pi-image.properties}, which the plugin build generates from
 * {@code build-image.sh} and VERSION.
 *
 * <p>That the generated resource actually reaches the shipped artifacts is
 * {@link PiAgentPackagingTest}'s job. This spec pins the other half: what the runner does when it
 * is there, and what it does when it is NOT - it must degrade to "no declaration", so a jar built
 * outside the Gradle build fails with core's "must declare a container" message rather than
 * with an exception out of a static initialiser during PF4J's extension loading.
 */
class PiAgentImageCoordinateTest extends Specification {

    private static InputStream stream(String content) {
        return new ByteArrayInputStream(content.getBytes('ISO-8859-1'))
    }

    def 'the runner declares the generated coordinate as its default container'() {
        expect: 'the resource the plugin build generates is on the runtime classpath'
        PiAgentRunner.getResourceAsStream(PiAgentRunner.IMAGE_RESOURCE) != null

        and: 'and the runner hands it out, so a pi agent needs no explicit `agent.container`'
        new PiAgentRunner().getDefaultContainer() ==~ $/\S+/nf-agent-pi:\S+/$
    }

    def 'a missing resource declares no container instead of throwing'() {
        given: 'the real lookup for a resource that is genuinely absent, not a stubbed stream'
        final absent = PiAgentRunner.getResourceAsStream('/META-INF/nf-agent-pi-image-absent.properties')

        expect: 'the classpath really has nothing there, so the null below is the case under test'
        absent == null

        and: 'which reads as "this runner declares no image", the shape a stale local jar has'
        PiAgentRunner.parseImageCoordinate(absent) == null
    }

    def 'an unusable resource also reads as no declaration'() {
        expect:
        PiAgentRunner.parseImageCoordinate(stream(content)) == expected

        where:
        content                                      || expected
        'image=reg.example.io/ns/nf-agent-pi:1.2.3'  || 'reg.example.io/ns/nf-agent-pi:1.2.3'
        // the value keeps its colons and slashes: load() splits on the FIRST separator, the `=`
        'image=  reg.example.io/ns/x:1.2.3  '        || 'reg.example.io/ns/x:1.2.3'
        'image='                                     || null
        'image=   '                                  || null
        '# generated header only, no image key'      || null
        ''                                           || null
    }
}
