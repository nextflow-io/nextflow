
workflow {
    FOO()
    BAR()
}

process FOO {
    exec:
    println "FOO: ${task.ext}"
}

process BAR {
    exec:
    println "BAR: ${task.ext}"
}
