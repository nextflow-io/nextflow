#!/usr/bin/env bash
# Fetch the pinned tree-sitter Groovy grammar and compile it into groovy.dylib.
#
# The grammar is NOT vendored: it is fetched at a pinned, immutable commit
# (git's SHA is a content hash, so the build is deterministic) into a gitignored
# build dir, then compiled with the system C compiler. Run once; the resulting
# groovy.dylib is reused by ast-grep via sgconfig.yml. Portable: compiles native
# on each host (Mach-O on macOS, ELF on Linux) — no committed per-platform binary.
#
#   ./build.sh
#   ast-grep scan -c .config/astgrep/sgconfig.yml
set -euo pipefail
cd "$(dirname "$0")"

# dekobon/tree-sitter-groovy, release v0.2.2, pinned by immutable commit SHA.
GRAMMAR_REPO="https://github.com/dekobon/tree-sitter-groovy.git"
GRAMMAR_COMMIT="8c70dc650be128673efd6f53bdc0ffba3ba4ce7d"
BUILD_DIR=".grammar-build"
OUT="groovy.dylib"

if [ ! -d "$BUILD_DIR/.git" ]; then
  rm -rf "$BUILD_DIR"
  git clone --quiet "$GRAMMAR_REPO" "$BUILD_DIR"
fi
git -C "$BUILD_DIR" fetch --quiet --tags origin
git -C "$BUILD_DIR" checkout --quiet "$GRAMMAR_COMMIT"

# Assert the checkout resolved to the exact pinned SHA (immutable content hash),
# not a moved tag/branch. If this ever fails, determinism is compromised — stop.
RESOLVED="$(git -C "$BUILD_DIR" rev-parse HEAD)"
if [ "$RESOLVED" != "$GRAMMAR_COMMIT" ]; then
  echo "error: expected $GRAMMAR_COMMIT but got $RESOLVED" >&2
  exit 1
fi

SRC="$BUILD_DIR/src"
[ -f "$SRC/parser.c" ] || { echo "error: $SRC/parser.c missing" >&2; exit 1; }

# scanner.c is plain C (not C++); compile parser + external scanner together.
# Drop -std=c11 if the compiler rejects it (clang defaults suffice).
gcc -std=c11 -shared -fPIC -O2 -I "$SRC" \
  "$SRC/parser.c" "$SRC/scanner.c" \
  -o "$OUT"

echo "built $OUT from ${GRAMMAR_COMMIT:0:7}"
