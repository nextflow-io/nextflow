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
import time
import xml.etree.ElementTree as ET
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


def parse_gradle_test_output(stdout, stderr):
    """Parse Gradle test output to extract structured information"""
    results = {
        'build_successful': False,
        'build_time': None,
        'tests_executed': False,
        'test_summary': {},
        'failure_summary': []
    }

    # Check if build was successful
    if 'BUILD SUCCESSFUL' in stdout:
        results['build_successful'] = True
        # Extract build time
        build_time_match = re.search(r'BUILD SUCCESSFUL in (.+)', stdout)
        if build_time_match:
            results['build_time'] = build_time_match.group(1)

    # Look for test execution indicators
    if 'Test worker' in stdout or 'test results' in stdout.lower():
        results['tests_executed'] = True

    # Parse test failure information from stderr and stdout
    combined_output = stdout + '\n' + stderr

    # Look for common test failure patterns
    failure_patterns = [
        r'(\w+Test) > (.+) FAILED',
        r'(\w+Test) > (.+) STANDARD_ERROR',
        r'FAILURE: Build failed with an exception',
        r'(.+): (.+Exception): (.+)',
        r'AssertionError: (.+)',
        r'NullPointerException',
        r'ComparisonFailure: (.+)'
    ]

    for line in combined_output.split('\n'):
        for pattern in failure_patterns:
            match = re.search(pattern, line)
            if match:
                results['failure_summary'].append({
                    'line': line.strip(),
                    'type': 'test_failure' if 'Test' in line else 'build_error',
                    'match_groups': match.groups()
                })

    return results


def parse_test_xml_results(module, test_class):
    """Parse XML test results for detailed statistics"""
    # Common path patterns for test results
    xml_paths = [
        f"modules/{module}/build/test-results/test/TEST-*.xml",
        f"build/test-results/test/TEST-*.xml",
        f"plugins/{module.replace('plugins:', '')}/build/test-results/test/TEST-*.xml"
    ]

    results = {
        'total_tests': 0,
        'passed': 0,
        'failed': 0,
        'skipped': 0,
        'execution_time': 0.0,
        'test_methods': [],
        'failed_tests': []
    }

    for xml_pattern in xml_paths:
        try:
            import glob
            xml_files = glob.glob(xml_pattern)
            for xml_file in xml_files:
                if test_class and test_class not in xml_file:
                    continue

                try:
                    tree = ET.parse(xml_file)
                    root = tree.getroot()

                    if root.tag == 'testsuite':
                        results['total_tests'] += int(root.get('tests', 0))
                        results['failed'] += int(root.get('failures', 0))
                        results['skipped'] += int(root.get('skipped', 0))
                        results['execution_time'] += float(root.get('time', 0))

                        for testcase in root.findall('testcase'):
                            test_name = testcase.get('name')
                            test_time = float(testcase.get('time', 0))

                            if testcase.find('failure') is not None:
                                failure = testcase.find('failure')
                                results['failed_tests'].append({
                                    'name': test_name,
                                    'time': test_time,
                                    'message': failure.get('message', ''),
                                    'type': failure.get('type', '')
                                })
                            else:
                                results['test_methods'].append({
                                    'name': test_name,
                                    'time': test_time,
                                    'status': 'passed'
                                })

                except ET.ParseError:
                    continue
                except Exception:
                    continue

        except Exception:
            continue

    results['passed'] = results['total_tests'] - results['failed'] - results['skipped']
    return results


def generate_failure_suggestions(gradle_output, xml_results, file_path):
    """Generate intelligent suggestions based on test failures"""
    suggestions = []

    # Analyze failure patterns
    if xml_results['failed_tests']:
        for failed_test in xml_results['failed_tests']:
            failure_type = failed_test.get('type', '').lower()
            failure_message = failed_test.get('message', '').lower()

            if 'nullpointerexception' in failure_type:
                suggestions.append(f"Check null handling in {failed_test['name']}")
            elif 'assertionerror' in failure_type or 'assertion' in failure_message:
                suggestions.append(f"Review test assertions in {failed_test['name']}")
            elif 'comparionfailure' in failure_type:
                suggestions.append(f"Check expected vs actual values in {failed_test['name']}")
            elif 'timeout' in failure_message:
                suggestions.append(f"Investigate performance issues in {failed_test['name']}")

    # Analyze Gradle output for build issues
    if gradle_output['failure_summary']:
        for failure in gradle_output['failure_summary']:
            if 'compilation' in failure['line'].lower():
                suggestions.append("Fix compilation errors before running tests")
            elif 'dependency' in failure['line'].lower():
                suggestions.append("Check project dependencies and classpath")

    # File-specific suggestions based on the file being edited
    filename = os.path.basename(file_path)
    if 'Cache' in filename:
        suggestions.append("Consider cache initialization and cleanup")
    elif 'Config' in filename:
        suggestions.append("Verify configuration parameter handling")
    elif 'Plugin' in filename:
        suggestions.append("Check plugin lifecycle and extension points")

    # Generic suggestions if no specific ones found
    if not suggestions:
        suggestions = [
            "Review recent changes for potential issues",
            "Check test data setup and teardown",
            "Verify mock objects and test dependencies"
        ]

    return suggestions[:3]  # Limit to 3 most relevant suggestions


def run_test_command(command, file_path, module, test_class):
    """Execute the test command with enhanced result parsing"""
    start_time = time.time()

    try:
        print(f"Running: {command}")

        result = subprocess.run(command, shell=True, capture_output=True,
                              text=True, timeout=180)  # 3 minute timeout

        execution_time = time.time() - start_time

        # Parse Gradle output for structured information
        gradle_output = parse_gradle_test_output(result.stdout, result.stderr)

        # Parse XML test results for detailed test information
        xml_results = parse_test_xml_results(module, test_class)

        # Generate intelligent suggestions for failures
        suggestions = []
        if result.returncode != 0:
            suggestions = generate_failure_suggestions(gradle_output, xml_results, file_path)

        return {
            'returncode': result.returncode,
            'stdout': result.stdout,
            'stderr': result.stderr,
            'execution_time': execution_time,
            'gradle_output': gradle_output,
            'xml_results': xml_results,
            'suggestions': suggestions,
            'command': command,
            'module': module,
            'test_class': test_class
        }

    except subprocess.TimeoutExpired:
        execution_time = time.time() - start_time
        return {
            'returncode': 1,
            'stdout': "",
            'stderr': "Test execution timed out after 3 minutes",
            'execution_time': execution_time,
            'gradle_output': {'build_successful': False, 'tests_executed': False},
            'xml_results': {},
            'suggestions': [
                "Check for infinite loops or blocking operations",
                "Review test timeout configurations",
                "Consider breaking large tests into smaller units"
            ],
            'command': command,
            'module': module,
            'test_class': test_class
        }
    except Exception as e:
        execution_time = time.time() - start_time
        return {
            'returncode': 1,
            'stdout': "",
            'stderr': f"Error running tests: {str(e)}",
            'execution_time': execution_time,
            'gradle_output': {'build_successful': False, 'tests_executed': False},
            'xml_results': {},
            'suggestions': [
                "Check Gradle installation and configuration",
                "Verify file permissions and paths",
                "Review system resources and availability"
            ],
            'command': command,
            'module': module,
            'test_class': test_class
        }


def main():
    try:
        input_data = json.load(sys.stdin)
    except json.JSONDecodeError as e:
        error_output = {
            "decision": "block",
            "reason": f"Test runner hook received invalid JSON input: {e}",
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
    if not file_path:
        sys.exit(0)

    # Determine test command
    test_command = determine_test_command(file_path)
    if not test_command:
        # Not a file we can test
        sys.exit(0)

    # Extract test information for enhanced reporting
    module, package, class_name = extract_test_info(file_path)
    test_class = f"{class_name}Test" if class_name and not class_name.endswith('Test') else class_name
    filename = os.path.basename(file_path)

    # Run the tests with enhanced parsing
    result = run_test_command(test_command, file_path, module, test_class)

    if result['returncode'] == 0:
        # Tests passed - generate rich success message
        xml_results = result['xml_results']
        gradle_output = result['gradle_output']

        # Create detailed success message
        if xml_results.get('total_tests', 0) > 0:
            test_count = xml_results['total_tests']
            test_time = xml_results.get('execution_time', result['execution_time'])
            build_time = gradle_output.get('build_time', f"{result['execution_time']:.1f}s")

            system_message = f"✓ Tests passed: {test_count} tests in {test_class or 'module'} ({build_time})"

            # Add method details if available
            if xml_results.get('test_methods'):
                method_names = [method['name'] for method in xml_results['test_methods'][:3]]
                if len(xml_results['test_methods']) > 3:
                    method_names.append(f"+ {len(xml_results['test_methods']) - 3} more")
                system_message += f"\n  Methods: {', '.join(method_names)}"

        else:
            # Fallback for when XML parsing doesn't work
            build_time = gradle_output.get('build_time', f"{result['execution_time']:.1f}s")
            system_message = f"✓ Tests passed for {filename} ({build_time})"

        output = {
            "suppressOutput": True,
            "systemMessage": system_message,
            "testResults": {
                "success": True,
                "totalTests": xml_results.get('total_tests', 0),
                "passed": xml_results.get('passed', 0),
                "failed": xml_results.get('failed', 0),
                "skipped": xml_results.get('skipped', 0),
                "executionTime": result['execution_time'],
                "buildTime": gradle_output.get('build_time'),
                "testClass": test_class,
                "module": module,
                "gradleTask": result['command'],
                "testMethods": xml_results.get('test_methods', []),
                "buildSuccessful": gradle_output.get('build_successful', True)
            }
        }
        print(json.dumps(output))
        sys.exit(0)
    else:
        # Tests failed - provide detailed error analysis
        xml_results = result['xml_results']
        gradle_output = result['gradle_output']
        suggestions = result['suggestions']

        # Create comprehensive error message
        if xml_results.get('failed_tests'):
            failed_count = len(xml_results['failed_tests'])
            total_count = xml_results.get('total_tests', failed_count)
            error_message = f"Tests failed: {failed_count}/{total_count} tests failed in {test_class or filename}"

            # Add specific failure details
            failure_details = []
            for failed_test in xml_results['failed_tests'][:2]:  # Show first 2 failures
                failure_type = failed_test.get('type', 'Unknown')
                failure_message = failed_test.get('message', '').split('\n')[0][:100]
                failure_details.append(f"  • {failed_test['name']}: {failure_type}")
                if failure_message:
                    failure_details.append(f"    {failure_message}")

            if len(xml_results['failed_tests']) > 2:
                failure_details.append(f"    + {len(xml_results['failed_tests']) - 2} more failures")

        elif not gradle_output.get('build_successful', False):
            error_message = f"Build failed for {filename}"
            failure_details = []

            # Extract build failure information
            if gradle_output.get('failure_summary'):
                for failure in gradle_output['failure_summary'][:2]:
                    failure_details.append(f"  • {failure['line'][:100]}")

            if result['stderr']:
                error_lines = [line.strip() for line in result['stderr'].split('\n')
                             if line.strip() and any(keyword in line.lower()
                                                   for keyword in ['error', 'failed', 'exception'])]
                failure_details.extend([f"  • {line[:100]}" for line in error_lines[:2]])

        else:
            error_message = f"Test execution failed for {filename}"
            failure_details = [f"  • {result['stderr'][:200]}" if result['stderr'] else "  • Unknown test failure"]

        # Combine error message with details
        full_error_message = error_message
        if failure_details:
            full_error_message += "\n" + "\n".join(failure_details)

        # Add execution context
        execution_context = f"\nCommand: {result['command']}\nDuration: {result['execution_time']:.1f}s"

        output = {
            "decision": "block",
            "reason": (full_error_message + execution_context)[:1500] + ("..." if len(full_error_message + execution_context) > 1500 else ""),
            "stopReason": f"Test failures in {filename}",
            "testFailure": {
                "success": False,
                "filename": filename,
                "testClass": test_class,
                "module": module,
                "totalTests": xml_results.get('total_tests', 0),
                "failed": xml_results.get('failed', 0),
                "passed": xml_results.get('passed', 0),
                "executionTime": result['execution_time'],
                "buildSuccessful": gradle_output.get('build_successful', False),
                "failedTests": xml_results.get('failed_tests', [])[:3],  # Limit to 3
                "suggestions": suggestions,
                "gradleTask": result['command'],
                "errorType": "test_failure" if xml_results.get('failed_tests') else "build_failure"
            }
        }
        print(json.dumps(output))
        sys.exit(0)


if __name__ == "__main__":
    main()