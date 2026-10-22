#!/usr/bin/env python3
#
# Copyright 2013-2026, Seqera Labs
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
"""
Clean up the Google Cloud Storage buckets used by the Nextflow tests.

For each --path gs://BUCKET/PREFIX, delete the objects under it older than --max-age-days.
Objects are listed and deleted page by page.

Dry run unless --delete is given.
Credentials: GOOGLE_APPLICATION_CREDENTIALS (service account key file).
Requires the `google-cloud-storage` package.
"""
import argparse
import sys
from datetime import datetime, timedelta, timezone

from google.cloud import storage

PAGE_SIZE = 1000   # max results per list call
BATCH_SIZE = 100   # max requests per batch call


def parse_path(value):
    if not value.startswith('gs://'):
        raise argparse.ArgumentTypeError(f'not a gs:// path: {value}')
    bucket, _, prefix = value[len('gs://'):].partition('/')
    if not prefix:
        raise argparse.ArgumentTypeError(f'refusing to clean a whole bucket: {value}')
    return bucket, prefix


def clean_path(client, bucket, prefix, cutoff, delete):
    name = f'gs://{bucket}/{prefix}'
    total = size = failed = 0
    blobs = client.list_blobs(bucket, prefix=prefix, page_size=PAGE_SIZE)
    for page_num, page in enumerate(blobs.pages, start=1):
        stale = [b for b in page if b.updated < cutoff]
        if delete:
            for i in range(0, len(stale), BATCH_SIZE):
                chunk = stale[i:i + BATCH_SIZE]
                try:
                    with client.batch(raise_exception=False):
                        for blob in chunk:
                            blob.delete()
                except Exception as e:
                    failed += len(chunk)
                    print(f'  failed batch in {name}: {e}', file=sys.stderr)
        total += len(stale)
        size += sum(b.size or 0 for b in stale)
        print(f'{name} page {page_num}: {len(stale)} stale objects (total {total}, {size / 1e9:.2f} GB)', flush=True)
    print(f'{name}: {total - failed} objects, {size / 1e9:.2f} GB; {failed} failed', flush=True)
    return failed


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--path', action='append', required=True, type=parse_path, metavar='gs://BUCKET/PREFIX',
                        help='path whose old objects are deleted (repeatable), e.g. gs://rnaseq-nf/scratch/')
    parser.add_argument('--max-age-days', type=int, default=7, help='only delete objects older than this (default: 7)')
    parser.add_argument('--delete', action='store_true', help='delete for real (default is a dry run)')
    args = parser.parse_args()

    client = storage.Client()
    cutoff = datetime.now(timezone.utc) - timedelta(days=args.max_age_days)
    print(f'[{"DELETING" if args.delete else "DRY RUN"}] objects older than {args.max_age_days} days', flush=True)

    failed = sum(clean_path(client, bucket, prefix, cutoff, args.delete) for bucket, prefix in args.path)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
