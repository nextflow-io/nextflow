#!/usr/bin/env python3
"""Validate every .data.json file under a directory against a JSON Schema.

Usage:
    validate.py --schema PATH --root DIR [--errors-per-file N] [--summary-only]
"""
import argparse
import json
import sys
from pathlib import Path

from jsonschema import Draft202012Validator


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--schema', required=True, type=Path, help='Path to JSON Schema file')
    ap.add_argument('--root', required=True, type=Path, help='Directory to scan for .data.json')
    ap.add_argument('--errors-per-file', type=int, default=3, help='Max errors to print per failing file')
    ap.add_argument('--summary-only', action='store_true', help='Skip per-file failure detail')
    args = ap.parse_args()

    schema = json.load(open(args.schema))
    validator = Draft202012Validator(schema)

    files = sorted(args.root.rglob('.data.json'))
    print(f'Schema: {args.schema}')
    print(f'Root:   {args.root}')
    print(f'Found {len(files)} .data.json file(s)\n')

    total_fail = 0
    total_parse = 0
    for f in files:
        try:
            doc = json.load(open(f))
        except json.JSONDecodeError as e:
            total_parse += 1
            if not args.summary_only:
                print(f'PARSE  {f.relative_to(args.root)}: {e}')
            continue
        errors = list(validator.iter_errors(doc))
        if errors:
            total_fail += 1
            if not args.summary_only:
                print(f'FAIL   {f.relative_to(args.root)}  ({len(errors)} errors)')
                for e in errors[:args.errors_per_file]:
                    path = '/'.join(str(p) for p in e.absolute_path) or '(root)'
                    print(f'         at {path}: {e.message[:160]}')

    print(f'\n--- summary ---')
    print(f'  files:        {len(files)}')
    print(f'  passed:       {len(files) - total_fail - total_parse}')
    print(f'  failed:       {total_fail}')
    print(f'  parse errors: {total_parse}')
    return 1 if (total_fail or total_parse) else 0


if __name__ == '__main__':
    sys.exit(main())
