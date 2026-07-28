# Use GitLab-provided links to paginate repository trees

- Authors: Gavin Elder
- Status: proposed
- Date: 2026-07-28
- Tags: scm gitlab pagination repository

## Summary

Nextflow should paginate GitLab repository-tree responses by following the
API-provided `Link` header. Offset pagination provides the broadest
compatibility, while keyset pagination is the preferred longer-term strategy
for large repositories.

This ADR proposes a direction only. It does not include or depend on any code
changes.

## Context

The repository directory traversal API was introduced in
[PR #6430](https://github.com/nextflow-io/nextflow/pull/6430). The GitLab
implementation uses:

```text
GET /api/v4/projects/:id/repository/tree
```

GitLab returns 20 entries by default and allows at most 100 entries per page.
The current `GitlabRepositoryProvider.listDirectory()` implementation makes
one request and parses one response, so it can return only the first 20 files
and directories.

Branches and tags already use Nextflow's shared
`invokeAndResponseWithPaging()` helper. That helper generates `page=N`
parameters and stops after 100 pages. At GitLab's default page size, branch and
tag results are therefore limited to 2,000 entries.

No open upstream Nextflow issue currently tracks GitLab repository-tree
pagination. Existing open GitLab issues concern authentication, subgroup
handling, and script selection rather than directory traversal.

## Problem statement

`listDirectory()` is expected to return the repository entries that match the
requested path and depth. For GitLab repositories containing more than 20
matching entries, the result is incomplete without warning.

Reusing the shared page-number helper would improve the immediate limit, but
would introduce a different truncation boundary:

- `per_page=100` combined with the helper's 100-page guard returns at most
  10,000 entries;
- the helper constructs subsequent URLs instead of following those returned by
  GitLab;
- it cannot use GitLab's keyset `page_token`; and
- it discovers completion by requesting an additional empty page rather than
  using response metadata.

Large monorepos can exceed 10,000 recursive tree entries, so silently moving
the limit from 20 to 10,000 does not provide complete traversal semantics.

## Goals

- Return complete GitLab repository-tree listings.
- Use pagination URLs supplied by GitLab.
- Avoid silently returning partial results.
- Continue to support self-managed GitLab installations where practical.
- Preserve existing authentication, retry, response-size, and error handling.
- Keep pagination state local to a single invocation.
- Make the behavior testable without requiring a live GitLab instance.

## Non-goals

- Implementing the proposal as part of this ADR.
- Changing `RepositoryProvider.listDirectory()` depth semantics.
- Streaming entries to callers; the current API returns a complete list.
- Changing pagination for other Git providers.
- Changing branch and tag pagination in the initial implementation.
- Guaranteeing snapshot consistency when a branch moves during traversal. A
  commit SHA should be used when a stable snapshot is required.

## Decision drivers

- **Correctness:** callers should not receive an apparently successful but
  incomplete tree.
- **Compatibility:** Nextflow is used with GitLab.com and self-managed GitLab
  versions that may not all support tree keyset pagination.
- **API conformance:** GitLab recommends following its pagination links instead
  of constructing subsequent URLs.
- **Maintainability:** GitLab-specific cursor behavior should not be forced
  into the existing page-number helper.
- **Safety:** an optional request or entry limit should fail explicitly rather
  than truncate silently.

## GitLab pagination options

### Offset pagination

GitLab's broadly supported pagination form is:

```text
?per_page=100&page=1
```

Responses can include a `Link` header with a URL whose relation is `next`.
Following that URL avoids assumptions about parameter ordering and termination.

Offset pagination remains subject to the GitLab instance's configured maximum
offset, so it may not be sufficient for very large trees.

### Keyset pagination

The repository-tree endpoint also supports:

```text
?pagination=keyset&per_page=100
```

The next URL in the `Link` header contains a `page_token`. Keyset pagination
does not depend on increasingly large offsets and is preferred for large
collections.

Repository-tree keyset pagination was introduced in GitLab 17.1. It should not
be required unconditionally until Nextflow's minimum supported self-managed
GitLab version is established.

## Proposed decision

Nextflow should use GitLab's `Link` header as the authoritative source of the
next repository-tree page.

The implementation should be delivered in two stages:

1. **Compatibility-first:** request `per_page=100`, use offset pagination, and
   follow each `rel="next"` link until none is returned.
2. **Keyset pagination:** use `pagination=keyset` where supported, with a
   deliberate fallback to offset pagination for older GitLab installations.

Fallback must occur before results from the first mode are accepted. Nextflow
must not switch pagination strategies midway through a listing.

## Potential implementation

### Expose response headers to repository providers

The current `RepositoryProvider.invoke()` API returns only a response body.
`GitlabRepositoryProvider` therefore cannot inspect pagination headers.

A future implementation could add a protected, header-aware method such as:

```groovy
protected HttpResponse<byte[]> invokeResponse(String api)
```

It would use the existing HTTP client and apply the same:

- authentication;
- retry handling;
- status validation; and
- `NXF_GIT_RESPONSE_MAX_LENGTH` validation.

`invokeBytes()` could delegate to this method and return its body, preserving
the existing interface. Alternatively, Nextflow could define a small response
value object containing the body, status, URI, and headers instead of exposing
`HttpResponse` to provider implementations.

The response should be returned directly. Storing a continuation URL in a
mutable provider field would introduce unnecessary shared state and possible
concurrency problems.

### Add linked pagination for GitLab trees

A GitLab-specific helper could:

1. send the initial repository-tree request;
2. parse and accumulate the JSON array;
3. extract the URL whose relation is `next` from the `Link` header;
4. request that URL verbatim; and
5. stop when no next relation is present.

Conceptually:

```groovy
List<Map> result = []
String nextUrl = initialUrl

while (nextUrl) {
    final response = invokeResponse(nextUrl)
    result.addAll(parseTreeEntries(response.body()))
    nextUrl = findNextLink(response.headers())
}
```

The link parser should:

- handle multiple comma-separated relations;
- match the relation value instead of assuming link order;
- treat an absent header as the end of pagination;
- preserve encoded query parameters such as `page_token`; and
- reject or deliberately resolve unexpected non-absolute URLs.

After pagination, the existing depth filtering, `RepositoryEntry` conversion,
and name sorting can remain unchanged.

### Limits and failure behavior

The implementation should not inherit the shared helper's fixed 100-page
limit. If Nextflow needs a configurable request, entry, memory, or response-size
limit, exceeding that limit should raise a clear exception. Returning a partial
tree as if it were complete is not acceptable.

## Alternatives considered

### Use the shared page-number helper

Adding `per_page=100` and calling `invokeAndResponseWithPaging()` is a small
change and raises the limit from 20 to 10,000 entries. It is not recommended
because it still truncates, does not use response metadata, and cannot support
keyset pagination.

### Remove the shared helper's page limit

This would allow more offset pages but would continue to generate pagination
URLs and remain subject to GitLab's maximum offset.

### Increase the shared helper's page limit

This only moves the truncation boundary and preserves silent partial results.

### Derive `page_token` from the final entry

This would couple Nextflow to GitLab cursor internals. The API-provided next
link is authoritative and should be followed instead.

### Clone or download the repository

A local tree walk avoids API pagination but changes the storage and performance
model of `listDirectory()` and defeats provider-native directory traversal.

## Consequences

### Positive

- GitLab directory traversal can return complete results.
- Pagination follows GitLab's documented protocol.
- Keyset pagination can support large trees efficiently.
- Explicit failure avoids hidden data loss when a configured limit is reached.
- Response-header access may support other provider APIs in the future.

### Negative

- The base HTTP abstraction needs a small extension.
- `Link` parsing adds provider-specific behavior and tests.
- Supporting old and new GitLab versions requires a compatibility policy.
- The API still accumulates the complete tree in memory.
- Offset mode can still encounter a server-configured maximum offset.

## Testing considerations

A future implementation should include unit tests for:

- a single page without a `Link` header;
- multiple pages using `rel="next"`;
- `prev`, `next`, and `last` links in varying orders;
- keyset links containing an encoded `page_token`;
- offset links containing `page=N`;
- malformed or non-HTTP next links;
- more than 100 simulated pages;
- an error on an intermediate page without returning partial results;
- fallback from unsupported keyset pagination before accepting results; and
- preservation of `ref`, `path`, `recursive`, authentication, and
  `per_page=100` through URLs returned by GitLab.

Conditional integration coverage can continue to use
`NXF_GITLAB_ACCESS_TOKEN`.

## Open questions

- What is the oldest self-managed GitLab version Nextflow intends to support?
- Should the first implementation use offset pagination only?
- How should unsupported keyset pagination be detected reliably?
- Should linked pagination become a reusable base abstraction or remain
  GitLab-specific?
- Should branch and tag pagination subsequently migrate to the same mechanism?
- Should a dedicated upstream issue be opened before implementation?
- Should entry or request limits be configurable, and which exception should
  report that a complete listing cannot be returned?

## Implementation status

No implementation is included with this ADR. If accepted, code and test changes
should be submitted separately and linked from this document.

## References

- [GitLab repository tree API](https://docs.gitlab.com/api/repositories/#list-repository-tree)
- [GitLab REST API pagination](https://docs.gitlab.com/api/rest/#pagination)
- [Nextflow directory traversal PR #6430](https://github.com/nextflow-io/nextflow/pull/6430)
- [Repository directory traversal ADR](20250929-repository-directory-traversal.md)
