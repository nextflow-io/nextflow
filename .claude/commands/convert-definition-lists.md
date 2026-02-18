---
description: Convert Markdown definition lists to heading format for any documentation page
---

# Convert Definition Lists to Headings

Convert Markdown definition list syntax to heading format across all documentation pages including CLI options, environment variables, configuration settings, and API parameters.

## Task

Search for and convert all definition list entries from the MyST definition list syntax:
```markdown
`term`
: Description text here.
```

To heading format:
```markdown
##### `term`

Description text here.
```

This applies to various documentation contexts:
- CLI options and flags
- Environment variables
- Configuration settings
- API/function parameters
- Type definitions
- Any other documented items using definition lists

## Conversion Rules

### Top-Level Definition Lists (h5 headings)

Use h5 headings (`#####`) for top-level items that are primary documentation entries.

1. **Simple definitions:**
   ```markdown
   `TERM_NAME`
   : Description text here.
   ```
   Converts to:
   ```markdown
   ##### `TERM_NAME`

   Description text here.
   ```

2. **With default values or metadata:**
   ```markdown
   `NXF_OPTS` (`-Xms512m`)
   : JVM options for Nextflow execution.
   ```
   Converts to:
   ```markdown
   ##### `NXF_OPTS` (`-Xms512m`)

   JVM options for Nextflow execution.
   ```

3. **Multi-line descriptions:**
   ```markdown
   `NXF_WORK`
   : Directory where working files are stored.

     Default: `./work`
   ```
   Converts to:
   ```markdown
   ##### `NXF_WORK`

   Directory where working files are stored.

   Default: `./work`
   ```

### Examples by Documentation Type

**Environment Variables:**
```markdown
`NXF_HOME`
: Nextflow home directory
```
→ `##### `NXF_HOME``

**CLI Options:**
```markdown
`-with-docker`
: Enable Docker container execution
```
→ `##### `-with-docker``

**Configuration Settings:**
```markdown
`process.executor`
: Defines the executor to use
```
→ `##### `process.executor``

## Embedded Definition Lists (h6 headings)

Embedded definition lists are indented and represent nested or sub-items within a section. These should be converted to h6 headings to preserve hierarchy.

### Pattern

```markdown
  `nestedItem`
  : Description text here (default: *value*)
```

### Conversion Format

**Primary conversion to h6 heading (default):**
```markdown
###### `nestedItem`

Description text here (default: *value*)
```

### Examples by Context

**Function Parameters:**
```markdown
  `maxDepth: int`
  : Maximum number of directory levels to visit (default: *no limit*)
```
→
```markdown
###### `maxDepth: int`

Maximum number of directory levels to visit (default: *no limit*)
```

**Configuration Sub-options:**
```markdown
  `retry.maxAttempts`
  : Maximum number of retry attempts
```
→
```markdown
###### `retry.maxAttempts`

Maximum number of retry attempts
```

**Type Properties:**
```markdown
  `name: String`
  : The name of the resource
```
→
```markdown
###### `name: String`

The name of the resource
```

### Alternative Formats

For very compact lists or where h6 headings are too prominent:

**Option A: Bold inline format:**
```markdown
- **`item`** - Description text here
```

**Option B: Simple list:**
```markdown
- `item`: Description text here
```

**When to use alternatives:** Only when the documentation context specifically requires a more compact format, or when there are many simple items in quick succession.

## Steps

1. **Search for definition lists** using Grep:
   - Top-level: lines with backtick-wrapped terms followed by `: ` on next line (no leading whitespace)
   - Embedded: indented lines (with leading spaces) with backtick-wrapped terms followed by `: `

2. **Determine heading level** based on indentation and context:
   - **No indentation** → Top-level item → h5 heading (`#####`)
   - **Indented (2-4 spaces)** → Nested item → h6 heading (`######`)

3. **For each match:**
   - Read surrounding context to capture full definition including multi-line descriptions
   - Note any special formatting within descriptions (code blocks, lists, emphasis)

4. **Convert using Edit tool:**
   - **Top-level**: Replace with `##### ` + term, remove `: ` prefix, add blank line
   - **Embedded**: Replace with `###### ` + term, remove indentation and `: ` prefix, add blank line
   - Preserve all formatting within description text

5. **Handle edge cases:**
   - Items with no description (rare but possible)
   - Nested content like code blocks, lists, or tables within descriptions
   - Multiple consecutive definition list items
   - Mixed indentation levels
   - Special characters in terms (backticks, colons, etc.)

## Pattern Matching

### Top-Level Definition Lists (h5)

**Pattern:**
```
`TERM`
: Description
```

**Characteristics:**
- Line starts with backtick (no leading whitespace)
- Term can be: environment variable, CLI option, config setting, type name, etc.
- Next line starts with `: ` (colon + space)
- Description may span multiple lines

**Common forms:**
- `` `ENV_VAR` `` - Environment variables
- `` `-option` `` - CLI flags
- `` `config.setting` `` - Configuration options
- `` `TypeName` `` - Type definitions
- `` `functionName()` `` - Function/method names

### Embedded Definition Lists (h6)

**Pattern:**
```
  `nestedTerm`
  : Description
```

**Characteristics:**
- Line starts with 2-4 spaces of indentation
- Term wrapped in backticks
- Next line indented with `: ` prefix
- Represents nested or sub-items

**Common forms:**
- `` `paramName: type` `` - Function parameters
- `` `property` `` - Object properties
- `` `subOption` `` - Nested configuration options
- `` `field` `` - Data structure fields

## Important Notes

### Formatting Preservation
- Preserve all formatting within descriptions: code blocks, lists, emphasis, links, etc.
- Keep default values in parentheses or emphasis exactly as-is
- Maintain relative indentation within multi-line descriptions

### Heading Levels
- **Top-level (no indent)** → h5 (`#####`) - one blank line between heading and description
- **Embedded (indented)** → h6 (`######`) - remove indentation, add one blank line between heading and description

### Common Uses Across Documentation
- **Environment variables** (env-vars pages)
- **CLI options and flags** (CLI reference pages)
- **Configuration settings** (config reference pages)
- **Function/method parameters** (API documentation)
- **Type properties** (type reference pages)
- **Feature flags** (feature documentation)

### Best Practices
- Work systematically top-to-bottom to avoid missing entries
- Read surrounding context to determine correct heading level
- Verify document structure and hierarchy after conversion
- Test rendering in documentation preview if available
- Use `replace_all: false` in Edit tool to handle each occurrence individually when patterns repeat

### Context Awareness
When uncertain about conversion format, consider:
- Where does this item appear in the document hierarchy?
- Is it a primary entry or a nested sub-item?
- What's the surrounding content structure?
- Would h6 preserve or break the logical flow?
