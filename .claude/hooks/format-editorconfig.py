#!/usr/bin/env python3
"""
EditorConfig enforcement hook for Nextflow development.
This hook applies editorconfig formatting rules after file edits.
"""

import hashlib
import json
import os
import subprocess
import sys
import time
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


def get_file_hash(file_path):
    """Get SHA256 hash of file contents"""
    try:
        with open(file_path, 'rb') as f:
            return hashlib.sha256(f.read()).hexdigest()
    except Exception:
        return None


def detect_formatting_changes(file_path, before_hash, after_hash):
    """Detect what types of formatting changes were made"""
    if before_hash == after_hash:
        return []

    changes = []

    try:
        # Read the file to analyze changes
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        # Detect common formatting changes
        if content.endswith('\n') and not content.endswith('\n\n'):
            changes.append('final_newline')

        if '\t' in content:
            changes.append('indentation')

        if content != content.rstrip():
            changes.append('trailing_whitespace')

        if '\r\n' in content:
            changes.append('line_endings')

        # If we can't detect specific changes, just note something changed
        if not changes:
            changes.append('formatting')

    except Exception:
        changes = ['formatting']

    return changes


def get_eclint_version():
    """Get eclint version if available"""
    try:
        result = subprocess.run(['eclint', '--version'], capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.strip()
    except Exception:
        pass
    return None


def format_with_editorconfig(file_path):
    """Apply editorconfig formatting to a file"""
    start_time = time.time()
    eclint_installed = False

    try:
        # Get file hash before formatting
        before_hash = get_file_hash(file_path)

        # Check if eclint is available
        result = subprocess.run(['which', 'eclint'], capture_output=True, text=True)
        if result.returncode != 0:
            print("eclint not found. Installing via npm...", file=sys.stderr)
            install_start = time.time()
            install_result = subprocess.run(['npm', 'install', '-g', 'eclint'],
                                          capture_output=True, text=True)
            if install_result.returncode != 0:
                execution_time = time.time() - start_time
                return {
                    'success': False,
                    'error': 'Failed to install eclint',
                    'error_details': install_result.stderr,
                    'execution_time': execution_time,
                    'suggestions': [
                        'Install Node.js and npm if not available',
                        'Check npm global install permissions',
                        'Try: sudo npm install -g eclint'
                    ]
                }
            eclint_installed = True
            print(f"eclint installed successfully in {time.time() - install_start:.1f}s", file=sys.stderr)

        # Apply editorconfig formatting
        format_result = subprocess.run(['eclint', 'fix', file_path],
                                     capture_output=True, text=True)

        # Get file hash after formatting
        after_hash = get_file_hash(file_path)
        execution_time = time.time() - start_time

        if format_result.returncode == 0:
            changes = detect_formatting_changes(file_path, before_hash, after_hash)
            eclint_version = get_eclint_version()

            return {
                'success': True,
                'changes_applied': changes,
                'changes_made': len(changes) > 0,
                'execution_time': execution_time,
                'eclint_installed': eclint_installed,
                'eclint_version': eclint_version,
                'file_changed': before_hash != after_hash
            }
        else:
            return {
                'success': False,
                'error': 'eclint formatting failed',
                'error_details': format_result.stderr,
                'execution_time': execution_time,
                'suggestions': [
                    'Check .editorconfig file syntax',
                    'Verify file permissions',
                    'Review eclint documentation'
                ]
            }

    except Exception as e:
        execution_time = time.time() - start_time
        return {
            'success': False,
            'error': f'Exception during formatting: {str(e)}',
            'execution_time': execution_time,
            'suggestions': [
                'Check file exists and is readable',
                'Verify system permissions',
                'Review hook configuration'
            ]
        }


def main():
    try:
        input_data = json.load(sys.stdin)
    except json.JSONDecodeError as e:
        error_output = {
            "decision": "block",
            "reason": f"EditorConfig hook received invalid JSON input: {e}",
            "suggestions": ["Check Claude Code hook configuration"]
        }
        print(json.dumps(error_output))
        sys.exit(0)

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

    # Apply editorconfig formatting with enhanced reporting
    result = format_with_editorconfig(file_path)
    filename = os.path.basename(file_path)

    if result['success']:
        # Generate rich success message
        if result['changes_made']:
            changes_text = ', '.join(result['changes_applied'])
            time_text = f" ({result['execution_time']:.1f}s)"
            system_message = f"✓ EditorConfig: {changes_text} applied to {filename}{time_text}"
        else:
            time_text = f" ({result['execution_time']:.1f}s)"
            system_message = f"✓ EditorConfig: no changes needed for {filename}{time_text}"

        # Add installation note if eclint was installed
        if result['eclint_installed']:
            system_message += " (eclint auto-installed)"

        output = {
            "suppressOutput": True,
            "systemMessage": system_message,
            "formattingResults": {
                "changesApplied": result['changes_applied'],
                "changesMade": result['changes_made'],
                "executionTime": result['execution_time'],
                "eclintInstalled": result['eclint_installed'],
                "eclintVersion": result.get('eclint_version'),
                "fileChanged": result['file_changed']
            }
        }
        print(json.dumps(output))
        sys.exit(0)
    else:
        # Enhanced error output with structured information
        error_message = f"EditorConfig formatting failed for {filename}: {result['error']}"

        output = {
            "decision": "block",
            "reason": error_message,
            "stopReason": f"Code formatting issues in {filename}",
            "formattingError": {
                "errorType": result['error'],
                "errorDetails": result.get('error_details', ''),
                "executionTime": result['execution_time'],
                "suggestions": result.get('suggestions', []),
                "filename": filename
            }
        }
        print(json.dumps(output))
        sys.exit(0)


if __name__ == "__main__":
    main()