// Custom Prism language definition for Nextflow DSL
// Based on Groovy but with Nextflow-specific keywords and patterns

export default function(Prism) {
  Prism.languages.nextflow = {
    'shebang': {
      pattern: /^#!.+/,
      alias: 'comment',
      greedy: true
    },
    'comment': [
      {
        pattern: /(^|[^\\])\/\*[\s\S]*?(?:\*\/|$)/,
        lookbehind: true,
        greedy: true
      },
      {
        pattern: /(^|[^\\:])\/\/.*/,
        lookbehind: true,
        greedy: true
      }
    ],
    'string': [
      // Triple-quoted strings
      {
        pattern: /("""|''')[\s\S]*?\1/,
        greedy: true,
        inside: {
          'interpolation': {
            pattern: /\$\{[^}]+\}|\$[\w.]+/,
            inside: {
              'punctuation': /^\$\{|\}$/,
              rest: null // Will be defined below
            }
          }
        }
      },
      // Dollar slashy strings
      {
        pattern: /\$\/(?:[^\/]|\/[^$])*\/\$/,
        greedy: true,
        inside: {
          'interpolation': {
            pattern: /\$\{[^}]+\}|\$[\w.]+/,
            inside: {
              'punctuation': /^\$\{|\}$/,
              rest: null
            }
          }
        }
      },
      // Slashy strings
      {
        pattern: /\/(?:[^\/\\\r\n]|\\.)*\//,
        greedy: true
      },
      // Regular strings
      {
        pattern: /(["'])(?:\\[\s\S]|(?!\1)[^\\])*\1/,
        greedy: true,
        inside: {
          'interpolation': {
            pattern: /\$\{[^}]+\}|\$[\w.]+/,
            inside: {
              'punctuation': /^\$\{|\}$/,
              rest: null
            }
          }
        }
      }
    ],
    // Nextflow-specific directives and keywords
    'nextflow-directive': {
      pattern: /\b(?:accelerator|afterScript|arch|array|beforeScript|cache|clusterOptions|conda|container|containerOptions|cpus|debug|disk|echo|errorStrategy|executor|ext|fair|label|machineType|maxErrors|maxForks|maxRetries|memory|module|penv|pod|publishDir|queue|resourceLabels|resourceLimits|scratch|secret|spack|stageInMode|stageOutMode|storeDir|tag|time)\b/,
      alias: 'property'
    },
    'nextflow-keyword': {
      pattern: /\b(?:process|workflow|params|input|output|script|shell|exec|when|channel|emit|take|main|tuple|path|val|file|env|stdin|stdout|stderr|include|from|addParams|each|collect|collectFile|combine|concat|count|countBy|cross|distinct|dump|filter|first|flatMap|flatten|groupBy|groupTuple|join|last|map|max|merge|min|mix|multiMap|randomSample|reduce|set|splitCsv|splitFasta|splitFastq|splitJson|splitText|subscribe|take|toInteger|toList|toSortedList|transpose|unique|until|view)\b/,
      alias: 'keyword'
    },
    'keyword': /\b(?:as|assert|break|case|catch|class|const|continue|def|default|do|else|enum|extends|finally|for|goto|if|implements|import|in|instanceof|interface|new|package|return|super|switch|this|throw|throws|trait|try|while|var)\b/,
    'boolean': /\b(?:true|false|null)\b/,
    'number': /\b0x[\da-f]+\b|(?:\b\d+(?:\.\d*)?|\B\.\d+)(?:e[+-]?\d+)?[glidf]?\b/i,
    'operator': {
      pattern: /(^|[^.])(?:~|==?~?|\?[.:]?|\*(?:[.=]|\*=?)?|\.[@&]|\.\.<|\.\.(?!\.)|-[-=>]?|\+[+=]?|!=?|<(?:<=?|=>?)?|>(?:>>?=?|=)?|&[&=]?|\|[|=]?|\/=?|\^=?|%=?)/,
      lookbehind: true
    },
    'punctuation': /\.+|[{}[\];(),:$]/
  };

  // Add interpolation rest
  Prism.languages.nextflow.string.forEach(function(pattern) {
    if (pattern.inside && pattern.inside.interpolation) {
      pattern.inside.interpolation.inside.rest = Prism.languages.nextflow;
    }
  });

  // Alias for common usage
  Prism.languages.nf = Prism.languages.nextflow;
}
