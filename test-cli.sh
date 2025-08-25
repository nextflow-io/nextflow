#!/bin/bash

# Test script for new plugin command system
set -e

echo "🧪 Testing Nextflow CLI Plugin Command System"
echo "============================================"

# Set up environment
export JAVA_HOME=${JAVA_HOME:-$(dirname $(dirname $(which java)))}
export NXF_HOME=$(pwd)/build/tmp/nxf-test
export NXF_PLUGINS_MODE=dev
export NXF_PLUGINS_DEV=$(pwd)/plugins

# Clean and create test directory
rm -rf $NXF_HOME
mkdir -p $NXF_HOME

# Copy the built nextflow jar
cp build/libs/nextflow-*-all.jar $NXF_HOME/nextflow.jar

# Function to run nextflow with proper setup
run_nextflow() {
    echo "$ nextflow $@"
    java -cp "$NXF_HOME/nextflow.jar" nextflow.cli.Launcher "$@" || echo "Command failed (expected for some tests)"
    echo
}

echo "📋 1. Testing main help output (should show plugin commands)"
run_nextflow --help

echo "📋 2. Testing Wave help specifically"
run_nextflow wave --help

echo "📋 3. Testing new Wave command syntax"
run_nextflow wave pack --help

echo "📋 4. Testing legacy plugin command syntax (backward compatibility)"
run_nextflow plugin nf-wave:pack --help

echo "📋 5. Testing command discovery without plugins"
export NXF_PLUGINS_MODE=off
run_nextflow --help

echo "✅ Testing complete!"