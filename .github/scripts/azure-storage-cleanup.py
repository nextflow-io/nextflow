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
Clean up the Azure storage account used by the Nextflow tests:

  - delete the blobs older than --max-age-days in the given containers, page by page
    (the containers are kept)
  - delete the `nf-az-blob-*` containers older than --max-age-days, leaked by tests
    that failed or timed out before their cleanup

Dry run unless --delete is given.
Credentials: AZURE_STORAGE_ACCOUNT_NAME and AZURE_STORAGE_ACCOUNT_KEY.
Requires the `azure-storage-blob` package.
"""
import argparse
import os
import sys
from datetime import datetime, timedelta, timezone

from azure.storage.blob import BlobServiceClient

PAGE_SIZE = 5000          # max results per list call
BATCH_SIZE = 256          # max sub-requests per blob batch call
TEST_CONTAINER_PREFIX = 'nf-az-blob-'


def clean_container(service, name, cutoff, delete):
    client = service.get_container_client(name)
    total = failed = size = 0
    for page_num, page in enumerate(client.list_blobs(results_per_page=PAGE_SIZE).by_page(), start=1):
        blobs = [b for b in page if b.last_modified < cutoff]
        names = [b.name for b in blobs]
        size += sum(b.size or 0 for b in blobs)
        if delete:
            for i in range(0, len(names), BATCH_SIZE):
                chunk = names[i:i + BATCH_SIZE]
                responses = client.delete_blobs(*chunk, delete_snapshots='include', raise_on_any_failure=False)
                for blob_name, resp in zip(chunk, responses):
                    # 404 means already gone, which is fine
                    if resp.status_code not in (202, 404):
                        failed += 1
                        print(f'  failed {resp.status_code}: {name}/{blob_name}', file=sys.stderr)
        total += len(names)
        print(f'{name} page {page_num}: {len(names)} stale blobs (total {total}, {size / 1e9:.2f} GB)', flush=True)
    print(f'{name}: {total - failed} blobs older than {cutoff:%Y-%m-%d}, {size / 1e9:.2f} GB; {failed} failed', flush=True)
    return failed


def delete_test_containers(service, cutoff, delete):
    found = failed = 0
    for container in service.list_containers(name_starts_with=TEST_CONTAINER_PREFIX):
        if container.last_modified >= cutoff:
            continue
        found += 1
        print(f'test container {container.name} (last modified {container.last_modified:%Y-%m-%d})', flush=True)
        if delete:
            try:
                service.delete_container(container.name)
            except Exception as e:
                failed += 1
                print(f'  failed: {container.name}: {e}', file=sys.stderr)
    print(f'{TEST_CONTAINER_PREFIX}*: {found - failed} containers older than {cutoff:%Y-%m-%d}; {failed} failed', flush=True)
    return failed


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--container', action='append', default=[], metavar='CONTAINER',
                        help='container whose old blobs are deleted (repeatable)')
    parser.add_argument('--max-age-days', type=int, default=7,
                        help=f'only delete blobs and {TEST_CONTAINER_PREFIX}* containers older than this (default: 7)')
    parser.add_argument('--delete', action='store_true', help='delete for real (default is a dry run)')
    args = parser.parse_args()

    account = os.environ['AZURE_STORAGE_ACCOUNT_NAME']
    service = BlobServiceClient(
        account_url=f'https://{account}.blob.core.windows.net',
        credential=os.environ['AZURE_STORAGE_ACCOUNT_KEY'])
    print(f'[{"DELETING" if args.delete else "DRY RUN"}] storage account {account}', flush=True)

    cutoff = datetime.now(timezone.utc) - timedelta(days=args.max_age_days)
    failed = sum(clean_container(service, name, cutoff, args.delete) for name in args.container)
    failed += delete_test_containers(service, cutoff, args.delete)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
