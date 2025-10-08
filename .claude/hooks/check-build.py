#!/usr/bin/env python3
"""
Build check hook for Nextflow development.
This hook runs a quick compilation check to catch syntax errors.
"""

import json
import os
import subprocess
import sys
import time
from pathlib import Path


def run_build_check():
    """Run a quick build check"""
    start_time = time.time()
    
    try:
        # Run make compile for quick syntax checking
        result = subprocess.run(
            ['make', 'compile'],
            capture_output=True,
            text=True,
            timeout=120  # 2 minute timeout
        )
        
        execution_time = time.time() - start_time
        
        if result.returncode == 0:
            return {
                'success': True,
                'execution_time': execution_time,
                'message': f"✓ Build check passed ({execution_time:.1f}s)"
            }
        else:
            # Extract error information
            error_lines = []
            for line in result.stdout.split('\n') + result.stderr.split('\n'):
                if any(keyword in line.lower() for keyword in ['error', 'failed', 'exception']):
                    error_lines.append(line)
            
            error_summary = '\n'.join(error_lines[-10:]) if error_lines else result.stderr[-500:]
            
            return {
                'success': False,
                'execution_time': execution_time,
                'error': f"Build check failed:\n{error_summary}",
                'suggestions': [
                    'Check for syntax errors in modified files',
                    'Review compilation errors above',
                    'Run `make compile` manually for full output'
                ]
            }
            
    except subprocess.TimeoutExpired:
        execution_time = time.time() - start_time
        return {
            'success': False,
            'execution_time': execution_time,
            'error': 'Build check timed out after 2 minutes',
            'suggestions': [
                'Check if build is hung',
                'Try running `make compile` manually'
            ]
        }
    except Exception as e:
        execution_time = time.time() - start_time
        return {
            'success': False,
            'execution_time': execution_time,
            'error': f'Exception during build check: {str(e)}',
            'suggestions': [
                'Verify make is installed',
                'Check build system configuration'
            ]
        }


def main():
    try:
        input_data = json.load(sys.stdin)
    except json.JSONDecodeError as e:
        error_output = {
            "decision": "block",
            "reason": f"Build check hook received invalid JSON input: {e}",
            "suggestions": ["Check Claude Code hook configuration"]
        }
        print(json.dumps(error_output))
        sys.exit(0)

    hook_event = input_data.get("hook_event_name", "")
    
    # Only run on Stop and SubagentStop events
    if hook_event not in ["stop-hook", "subagent-stop-hook"]:
        sys.exit(0)

    # Run build check
    result = run_build_check()
    
    if result['success']:
        output = {
            "suppressOutput": True,
            "systemMessage": result['message']
        }
        print(json.dumps(output))
        sys.exit(0)
    else:
        # Block with error information
        output = {
            "decision": "block",
            "reason": result['error'],
            "stopReason": "Build compilation failed",
            "buildError": {
                "executionTime": result['execution_time'],
                "suggestions": result.get('suggestions', [])
            }
        }
        print(json.dumps(output))
        sys.exit(0)


if __name__ == "__main__":
    main()
