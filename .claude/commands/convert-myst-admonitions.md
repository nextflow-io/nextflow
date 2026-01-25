---
description: Convert MyST-style version admonitions to JSX tag format
---

# Convert MyST Admonitions to JSX Tags

Convert all version admonitions in MDX files from the old MyST colon-fence format to the new JSX tag-based format.

## Task

Search for and convert all instances of:
- `:::{versionadded} X.Y.Z` → `<AddedInVersion version="X.Y.Z" />`
- `:::{deprecated} X.Y.Z` → `<DeprecatedInVersion version="X.Y.Z" />`
- `:::{versionchanged} X.Y.Z` → `<ChangedInVersion version="X.Y.Z" />`

## Conversion Rules

1. **Simple admonitions (empty body):**
   ```markdown
   :::{versionadded} X.Y.Z
   :::
   ```
   Converts to:
   ```markdown
   <AddedInVersion version="X.Y.Z" />
   ```

2. **Admonitions with description (multi-line body):**
   ```markdown
   :::{deprecated} X.Y.Z
   Some description text here.
   :::
   ```
   Converts to:
   ```markdown
   <DeprecatedInVersion version="X.Y.Z">
   Some description text here.
   </DeprecatedInVersion>
   ```

## Steps

1. Use Grep to find all occurrences of `:::{versionadded}`, `:::{deprecated}`, and `:::{versionchanged}` patterns
2. Read the context around each match to understand the full admonition structure
3. Convert each admonition using Edit tool with appropriate old_string/new_string
4. For duplicate patterns that appear multiple times, include enough surrounding context to make each match unique
5. Verify the file has the correct imports at the top:
   ```javascript
      ```
6. After all conversions, verify no old-format admonitions remain using Grep

## Important Notes

- Preserve all surrounding formatting, spacing, and content exactly as-is
- Handle each occurrence individually with sufficient context for unique matching
- Use self-closing tags (`/>`) for empty admonitions
- Use container tags with content for admonitions with descriptions
- Work systematically from top to bottom of the file to avoid confusion with line number shifts