def global_var = 'global'

def method() {
    println("Method executed!")
}

workflow {
    params.outdir = 'results'

    def local = 'local'

    workflow.onComplete {
        println("local variable: ${local}")
        println("param output dir: ${params.outdir}")
        method() 
        println("Executed script : ${workflow.scriptFile}")
        // expected failure
        try {
                println("global variable: ${global}")
        }catch (groovy.lang.MissingPropertyException e){
                println("Can't access to a global variable")
        }
	try {
                println("workflow private variable: ${events}")
        }catch (Exception e){
                println("Can't access to a private workflow property")
        }
        if( workflow.success )
            println("Success!")
        else
            println("Failure!")
    }
    
    println("Workflow Executed!")
}
