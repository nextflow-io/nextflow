#!/usr/bin/env python3
"""
Test runner hook for Nextflow development.
This hook runs appropriate tests when source files or test files are edited.
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path


def extract_test_info(file_path):
    """Extract module and test class info from file path"""
    path = Path(file_path)

    # Check if it's in a module directory
    module_match = None
    for part in path.parts:
        if part in ['nextflow', 'nf-commons', 'nf-httpfs', 'nf-lang', 'nf-lineage']:
            module_match = part
            break
        # Handle plugin modules
        if part.startswith('nf-') and 'plugins' in path.parts:
            module_match = f"plugins:{part}"
            break

    if not module_match:
        return None, None, None

    # Extract the class/package info
    src_parts = list(path.parts)

    # Find where the package structure starts
    package_start_idx = -1
    for i, part in enumerate(src_parts):
        if part in ['groovy', 'java'] and i > 0 and src_parts[i-1] in ['main', 'test']:
            package_start_idx = i + 1
            break

    if package_start_idx == -1:
        return module_match, None, None

    # Get package and class name
    package_parts = src_parts[package_start_idx:-1]
    package = '.'.join(package_parts) if package_parts else None

    class_name = path.stem

    return module_match, package, class_name


def determine_test_command(file_path):
    """Determine the appropriate test command based on the file being edited"""
    path = Path(file_path)

    # Only process Groovy and Java files in modules
    if path.suffix not in ['.groovy', '.java']:
        return None

    # Must be in modules or plugins directory
    if 'modules' not in path.parts and 'plugins' not in path.parts:
        return None

    module, package, class_name = extract_test_info(file_path)
    if not module or not class_name:
        return None

    # If it's already a test file, run it directly
    if class_name.endswith('Test'):
        test_pattern = f"*{class_name}"
        return f"./gradlew :{module}:test --tests \"{test_pattern}\""

    # If it's a source file, look for corresponding test
    test_class = f"{class_name}Test"
    test_pattern = f"*{test_class}"

    return f"./gradlew :{module}:test --tests \"{test_pattern}\""


def run_test_command(command):
    """Execute the test command"""
    try:
        print(f"Running: {command}")

        result = subprocess.run(command, shell=True, capture_output=True,
                              text=True, timeout=180)  # 3 minute timeout

        return result.returncode, result.stdout, result.stderr

    except subprocess.TimeoutExpired:
        return 1, "", "Test execution timed out after 3 minutes"
    except Exception as e:
        return 1, "", f"Error running tests: {str(e)}"


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
    if not file_path:
        sys.exit(0)

    # Determine test command
    test_command = determine_test_command(file_path)
    if not test_command:
        # Not a file we can test
        sys.exit(0)

    # Run the tests
    returncode, stdout, stderr = run_test_command(test_command)

    if returncode == 0:
        # Tests passed
        lines = stdout.split('\n')
        test_results = [line for line in lines if 'test' in line.lower() and ('passed' in line.lower() or 'success' in line.lower())]

        message = f"✓ Tests passed for {os.path.basename(file_path)}"
        if test_results:
            # Show a summary of the last few relevant lines
            summary = '\n'.join(test_results[-3:]) if len(test_results) > 3 else '\n'.join(test_results)
            message += f"\n{summary}"

        output = {
            "suppressOutput": True,
            "systemMessage": message
        }
        print(json.dumps(output))
        sys.exit(0)
    else:
        # Tests failed - show error to Claude for potential fixing
        error_msg = f"Tests failed for {os.path.basename(file_path)}:\n"

        # Extract useful error information
        if stderr:
            error_msg += f"Error output:\n{stderr[:500]}\n"

        if stdout:
            # Look for test failure information in stdout
            lines = stdout.split('\n')
            failure_lines = [line for line in lines
                           if any(keyword in line.lower()
                                for keyword in ['failed', 'error', 'exception', 'assertion'])]

            if failure_lines:
                error_msg += f"Test failures:\n" + '\n'.join(failure_lines[-5:])

        # Use JSON output to provide feedback to Claude
        output = {
            "decision": "block",
            "reason": error_msg[:1000] + ("..." if len(error_msg) > 1000 else "")
        }
        print(json.dumps(output))
        sys.exit(0)


if __name__ == "__main__":
    main()