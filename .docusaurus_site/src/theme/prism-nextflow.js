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
      pattern: /\b(?:process|workflow|params|input|output|stage|topic|publish|script|shell|exec|when|channel|emit|take|main|tuple|path|val|file|env|stdin|stdout|stderr|include|from)\b/,
      alias: 'keyword'
    },
    // Standard types used by static typing, including generic types such as
    // `Channel<Path>` and `Map<String, Path>`.
    'class-name': {
      pattern: /\b(?:Bag|Boolean|Channel|Duration|Float|Integer|Iterable|List|Map|MemoryUnit|Path|Record|Set|String|Tuple|Value|VersionNumber)\b(?:<[^<>\n]*(?:<[^<>\n]*>[^<>\n]*)*>)?/,
      inside: {
        'class-name': /\b(?:Bag|Boolean|Channel|Duration|Float|Integer|Iterable|List|Map|MemoryUnit|Path|Record|Set|String|Tuple|Value|VersionNumber)\b/,
        'punctuation': /[<>,]/
      }
    },
    'keyword': /\b(?:as|assert|catch|def|else|enum|if|in|instanceof|new|record|return|throw|try)\b/,
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
