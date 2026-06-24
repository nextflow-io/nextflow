nextflow.enable.types = true

/*
 * ============================================================================
 *  EXAMPLE: module_run + filesystem capability tools
 * ============================================================================
 *
 *  This minimal example demonstrates both agent capability tools together:
 *
 *    tools 'module_run', 'filesystem'
 *
 *  `module_run` — exposes every in-scope process as its OWN tool (here just
 *  the `word_stats` tool). The LLM calls it with `{"text":"..."}` (the tool's
 *  enforced parameters schema) and receives the JSON output back.
 *
 *  `filesystem` — a sandboxed read/write/list/exists tool scoped to the
 *  agent's per-invocation work directory. After calling `module_run` the
 *  agent writes a small JSON report to the work dir via `filesystem`, then
 *  reads it back to verify and includes the file path in its answer.
 *
 *  The process itself is a trivial `exec:` block (no container needed)
 *  so the example runs fully offline without any external registry or
 *  container runtime.
 */

/*
 * Count words and characters in the input text. Returns a small JSON string
 * so the agent can read the stats directly from the word_stats tool result.
 */
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

/*
 * The agent:
 *  1. Calls the `word_stats` tool to compute text statistics.
 *  2. Writes a JSON report file to the work dir via the `filesystem` tool.
 *  3. Reads the file back to confirm it was written.
 *  4. Replies with a one-line summary including the report path.
 */
agent analyst {
    model 'openai/gpt-5-mini'
    instruction '''\
        You analyse text. Work step by step:
          1. Call the word_stats tool with {"text":"<the text>"} to
             get the word and character counts.
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
