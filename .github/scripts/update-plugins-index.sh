#!/usr/bin/env bash
set -e

# change to the project root
cd "$(dirname "$0")/../.."

echo "Updating plugins index"

./gradlew plugins:publishIndex

echo "Done"
