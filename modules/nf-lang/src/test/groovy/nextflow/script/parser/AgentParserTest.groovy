/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
package nextflow.script.parser

import nextflow.script.ast.AgentNode
import spock.lang.Specification
import test.TestUtils

class AgentParserTest extends Specification {

    def setupSpec() {
        TestUtils.beforeSpec()
    }

    def 'should parse a minimal agent definition'() {
        when:
        def script = TestUtils.parse('''\
            nextflow.enable.types = true

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()
                maxIterations 20

                input:
                    val question

                output:
                    val plan

                prompt:
                """
                Question: ${question}
                """
            }
            '''.stripIndent())

        then:
        script.agents.size() == 1
        def node = script.agents[0] as AgentNode
        node.name == 'eval_agent'
        node.inputs.length == 1
        node.inputs[0].name == 'question'
    }
}
