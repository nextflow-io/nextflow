package nextflow.config

import groovy.transform.MapConstructor
import groovy.transform.TupleConstructor
import jdk.nashorn.internal.objects.annotations.Constructor
import nextflow.Global
import nextflow.Session
import nextflow.script.ExecutionStack

import java.nio.file.Path


class ConfigRunPlan {

    @TupleConstructor
    static class RunPlan {
        Boolean required
        Boolean compare_output
        Boolean in_channel_override
    }

    static set

    static getRunPlanForProcess(String processName) {
        def session = Global.session as Session
        def map = session.runPlan
        def res = map.get(processName)
        if (res)
            return res
        def workflow = ExecutionStack.workflowStack().reverse().find {
            map.get(it.name) !=null
        }
        if (workflow)
            return map.get(workflow.name)
        return map.get('all')
    }

    static LinkedHashMap<String, RunPlan> parse(String plan) {
        def res = [:]
        def parse_line = { String line, int num ->
            println(line)
            if (!(line =~ /^\s*(:?#.*)?$/)) {
                def accepted_format = /^\s*(\w+)\s*(T|\*)?\s*(<)?/
                def matches = line =~ accepted_format
                if (!matches)
                    throw new Error("error in worklow plan:\n ${plan} line ${num}:\n ${line}\n"
                            + " line must match following regexp: ${accepted_format}")
                def m = matches[0]
                def process = m[1]
                def _plan = new RunPlan(
                        m[2] != null,
                        m[2] == 'T',
                        m[3] == '<')
                res[process] = _plan
            }
        }
        println("RUN PLAN")
        println(res)
        Path planFile = null
        if (plan !=~ /.*(:?T|\*|<|#|;).*/ && (planFile = nextflow.Nextflow.file(null, plan) as Path).exists())
            planFile.eachLine(parse_line)
        else
            plan.split(/\n|;/).eachWithIndex(parse_line)
        res
    }
}
