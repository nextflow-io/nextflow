nextflow.enable.types = true

// Demonstrates both capability tools together (`tools 'module_run', 'filesystem'`):
// the agent runs the local `word_stats` tool, then writes/reads a report file in its
// sandboxed work dir. Fully offline (exec process, no container). See README.md.

// Counts words/chars; returns a small JSON string the agent reads from the tool result.
process word_stats {
    input:
        text: String
    output:
        stats: String
    exec:
        def words = text.trim().split(/\s+/).length
        def chars  = text.length()
        stats = "{\"words\":${words},\"chars\":${chars}}"
}

agent analyst {
    model 'openai/gpt-5-mini'
    instruction '''\
        You analyse text. Work step by step:
          1. Call the word_stats tool to get the word and character counts.
          2. Write the statistics as a JSON file named "report.json" in the work directory
             using the filesystem write tool.
          3. Read "report.json" back using the filesystem read tool to confirm it exists.
          4. Reply: "Stats: <words> words, <chars> chars. Report written to <absolute path>."
        '''.stripIndent()

    tools 'module_run', 'filesystem'

    input:
        text: String
    output:
        summary: String

    prompt:
    """
    Analyse the following text:
    ${text}
    """
}

workflow {
    analyst(channel.of(
        'The quick brown fox jumps over the lazy dog'
    ))
    .view { s -> "RESULT: ${s}" }
}
