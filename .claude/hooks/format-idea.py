#!/usr/bin/env python3
"""
IntelliJ IDEA formatter hook for Nextflow development.
This hook applies IDEA code formatting to Groovy files after edits.
"""

import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path


def is_groovy_file(file_path):
    """Check if the file is a Groovy source file"""
    if not file_path:
        return False

    path = Path(file_path)

    # Skip build directories, .git, etc.
    if any(part.startswith('.') or part == 'build' for part in path.parts):
        return False

    # Only process Groovy files
    return path.suffix.lower() == '.groovy'


def get_file_hash(file_path):
    """Get SHA256 hash of file contents"""
    try:
        with open(file_path, 'rb') as f:
            return hashlib.sha256(f.read()).hexdigest()
    except Exception:
        return None


def find_idea_sh():
    """Locate idea.sh command"""
    # Check environment variable first
    idea_sh = os.environ.get('IDEA_SH')
    if idea_sh and os.path.exists(idea_sh):
        return idea_sh

    # Common installation paths for macOS
    common_paths = [
        '/Applications/IntelliJ IDEA.app/Contents/MacOS/idea',
        '/Applications/IntelliJ IDEA CE.app/Contents/MacOS/idea',
        '/Applications/IntelliJ IDEA Ultimate.app/Contents/MacOS/idea',
        # Common paths for Linux
        '/usr/local/bin/idea.sh',
        '/opt/idea/bin/idea.sh',
        '/snap/intellij-idea-community/current/bin/idea.sh',
        '/snap/intellij-idea-ultimate/current/bin/idea.sh',
        # Check user's home for JetBrains Toolbox installations
        os.path.expanduser('~/.local/share/JetBrains/Toolbox/apps/IDEA-U/bin/idea.sh'),
        os.path.expanduser('~/.local/share/JetBrains/Toolbox/apps/IDEA-C/bin/idea.sh'),
    ]

    for path in common_paths:
        if os.path.exists(path):
            return path

    # Try to find in PATH
    try:
        result = subprocess.run(['which', 'idea'], capture_output=True, text=True)
        if result.returncode == 0:
            idea_path = result.stdout.strip()
            if os.path.exists(idea_path):
                return idea_path
    except Exception:
        pass

    return None


def format_with_idea(file_path, idea_sh_path):
    """Apply IDEA formatting to a Groovy file"""
    start_time = time.time()

    try:
        # Get file hash before formatting
        before_hash = get_file_hash(file_path)

        # Run IDEA formatter
        # Using format.sh script which is typically alongside idea.sh
        idea_dir = os.path.dirname(idea_sh_path)
        format_script = None

        # Try to find format.sh in the same directory
        possible_format_scripts = [
            os.path.join(idea_dir, 'format.sh'),
            os.path.join(os.path.dirname(idea_dir), 'bin', 'format.sh'),
        ]

        for script in possible_format_scripts:
            if os.path.exists(script):
                format_script = script
                break

        # If no format.sh found, use idea.sh with format command
        if format_script:
            command = [format_script, '-allowDefaults', file_path]
        else:
            command = [idea_sh_path, 'format', '-allowDefaults', file_path]

        format_result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=30
        )

        # Get file hash after formatting
        after_hash = get_file_hash(file_path)
        execution_time = time.time() - start_time

        # Check if formatting succeeded
        # IDEA formatter typically returns 0 on success
        if format_result.returncode == 0:
            return {
                'success': True,
                'changes_made': before_hash != after_hash,
                'execution_time': execution_time,
                'file_changed': before_hash != after_hash,
                'idea_path': idea_sh_path,
                'command': ' '.join(command)
            }
        else:
            return {
                'success': False,
                'error': 'IDEA formatting failed',
                'error_details': format_result.stderr or format_result.stdout,
                'execution_time': execution_time,
                'command': ' '.join(command),
                'suggestions': [
                    'Verify IDEA installation and configuration',
                    'Check file permissions',
                    'Ensure IDEA is not already running the formatter'
                ]
            }

    except subprocess.TimeoutExpired:
        execution_time = time.time() - start_time
        return {
            'success': False,
            'error': 'IDEA formatting timed out after 30 seconds',
            'execution_time': execution_time,
            'suggestions': [
                'Check if IDEA is responsive',
                'Try closing and reopening IDEA',
                'Check system resources'
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
            "reason": f"IDEA formatter hook received invalid JSON input: {e}",
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
    if not file_path or not is_groovy_file(file_path):
        sys.exit(0)

    # Check if file exists after the edit
    if not os.path.exists(file_path):
        sys.exit(0)

    # Find IDEA installation
    idea_sh_path = find_idea_sh()
    if not idea_sh_path:
        # If IDEA is not found, silently skip (don't block the workflow)
        # This allows the hook to be optional
        output = {
            "suppressOutput": True,
            "systemMessage": "⚠ IDEA formatter: idea.sh not found (set IDEA_SH environment variable)",
            "formattingResults": {
                "skipped": True,
                "reason": "IDEA not found"
            }
        }
        print(json.dumps(output))
        sys.exit(0)

    # Apply IDEA formatting
    result = format_with_idea(file_path, idea_sh_path)
    filename = os.path.basename(file_path)

    if result['success']:
        # Generate success message
        if result['changes_made']:
            time_text = f" ({result['execution_time']:.1f}s)"
            system_message = f"✓ IDEA formatter: applied to {filename}{time_text}"
        else:
            time_text = f" ({result['execution_time']:.1f}s)"
            system_message = f"✓ IDEA formatter: no changes needed for {filename}{time_text}"

        output = {
            "suppressOutput": True,
            "systemMessage": system_message,
            "formattingResults": {
                "changesMade": result['changes_made'],
                "executionTime": result['execution_time'],
                "ideaPath": result['idea_path'],
                "fileChanged": result['file_changed']
            }
        }
        print(json.dumps(output))
        sys.exit(0)
    else:
        # Enhanced error output
        error_message = f"IDEA formatting failed for {filename}: {result['error']}"

        output = {
            "decision": "block",
            "reason": error_message,
            "stopReason": f"Code formatting issues in {filename}",
            "formattingError": {
                "errorType": result['error'],
                "errorDetails": result.get('error_details', ''),
                "executionTime": result['execution_time'],
                "suggestions": result.get('suggestions', []),
                "filename": filename,
                "command": result.get('command', '')
            }
        }
        print(json.dumps(output))
        sys.exit(0)


if __name__ == "__main__":
    main()
