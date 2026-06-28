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
package nextflow.agent

import spock.lang.Requires
import spock.lang.Specification

/**
 * End-to-end test of the agent SKILLS path through the real OpenAI integration.
 * The distinctive `[SEQ-REPORT v1]` output format is defined ONLY inside the skill
 * content (never in the instruction or prompt), so its presence in the final answer
 * proves the model saw the available-skills catalog, called {@code activate_skill}
 * to read the instructions, and followed them. Skipped when OPENAI_API_KEY is unset.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Requires({ System.getenv('OPENAI_API_KEY') })
class AgentSkillEndToEndTest extends Specification {

    def 'should activate a skill and follow its instructions end-to-end against OpenAI'() {
        given: 'the sequence-report skill (same content as the POC example)'
        def skill = new SkillDescriptor(
            'sequence-report',
            'Format a sequencing or genome-assembly summary as a standardized QC report. Use this skill whenever the user asks for a sequence, assembly, or read-QC summary or verdict.',
            '''When producing a sequencing or assembly summary, format the answer EXACTLY as follows and nothing else:
[SEQ-REPORT v1]
STATUS: <PASS | WARN | FAIL>
METRICS: <comma-separated list of the key metrics you were given>
NOTE: <one short sentence of interpretation>''',
            [])

        and: 'a skills-only request whose instruction does NOT mention the report format'
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini',
            instruction: 'You summarize sequencing and genome-assembly results for bioinformaticians.',
            prompt: 'Summarize this assembly: N50 = 45 kb, total length = 5.1 Mb, GC = 50.8%. Is it good enough to proceed?',
            maxIterations: 5, tools: [], outputSchema: null, inputJson: null,
            toolSpecs: null, dispatch: null, requestTimeoutSeconds: 0, goal: null,
            skills: [skill])

        when:
        def answer = new LangChainAgentRunner().run(req)

        then: 'the answer follows the skill-dictated format — proving activate_skill fired and was followed'
        answer.contains('[SEQ-REPORT v1]')
        answer.contains('STATUS:')
    }
}
