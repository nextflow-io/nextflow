#!/usr/bin/env python3
"""
EditorConfig enforcement hook for Nextflow development.
This hook applies editorconfig formatting rules after file edits.
"""

import json
import os
import subprocess
import sys
from pathlib import Path


def is_source_file(file_path):
    """Check if the file should be formatted"""
    if not file_path:
        return False

    # Only format source code files
    extensions = {'.groovy', '.java', '.gradle', '.md', '.txt', '.yml', '.yaml', '.json'}
    path = Path(file_path)

    # Skip build directories, .git, etc.
    if any(part.startswith('.') or part == 'build' for part in path.parts):
        return False

    return path.suffix.lower() in extensions


def format_with_editorconfig(file_path):
    """Apply editorconfig formatting to a file"""
    try:
        # Check if eclint is available
        result = subprocess.run(['which', 'eclint'], capture_output=True, text=True)
        if result.returncode != 0:
            print("eclint not found. Installing via npm...", file=sys.stderr)
            install_result = subprocess.run(['npm', 'install', '-g', 'eclint'],
                                          capture_output=True, text=True)
            if install_result.returncode != 0:
                return False, "Failed to install eclint"

        # Apply editorconfig formatting
        format_result = subprocess.run(['eclint', 'fix', file_path],
                                     capture_output=True, text=True)

        if format_result.returncode == 0:
            return True, f"Applied editorconfig formatting to {file_path}"
        else:
            return False, f"eclint failed: {format_result.stderr}"

    except Exception as e:
        return False, f"Error formatting file: {str(e)}"


def main():
    try:
        input_data = json.load(sys.stdin)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON input: {e}", file=sys.stderr)
        sys.exit(1)

    hook_event = input_data.get("hook_event_name", "")
    tool_name = input_data.get("tool_name", "")
    tool_input = input_data.get("tool_input", {})

    # Only process Edit, Write, MultiEdit tools
    if tool_name not in ["Edit", "Write", "MultiEdit"]:
        sys.exit(0)

    file_path = tool_input.get("file_path", "")
    if not file_path or not is_source_file(file_path):
        sys.exit(0)

    # Check if file exists after the edit
    if not os.path.exists(file_path):
        sys.exit(0)

    success, message = format_with_editorconfig(file_path)

    if success:
        # Use JSON output to suppress the normal stdout display
        output = {
            "suppressOutput": True,
            "systemMessage": f"✓ EditorConfig formatting applied to {os.path.basename(file_path)}"
        }
        print(json.dumps(output))
        sys.exit(0)
    else:
        # Non-blocking error - show message but don't fail
        print(f"Warning: {message}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()