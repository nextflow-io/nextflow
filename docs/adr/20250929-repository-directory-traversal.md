# ADR: Repository Directory Traversal API

**Date**: 2025-09-29  
**Status**: Accepted  
**Context**: Need for standardized directory listing across Git hosting providers

## Decision

Add a `listDirectory(String path, int depth)` method to the `RepositoryProvider` abstraction so that directory traversal works the same way on every Git hosting platform.

## Context

Nextflow needs to explore repository directory structures on several Git hosting providers (GitHub, GitLab, Bitbucket, Azure DevOps, Gitea) without cloning the whole repository. Each provider offers different API capabilities and constraints for directory listing.

## Technical Implementation

### Core Algorithm

All providers follow a consistent pattern:
1. **Path Resolution**: Normalize path to provider API format
2. **Strategy Selection**: Choose a recursive or iterative approach based on API capabilities
3. **HTTP Request**: Execute provider-specific API calls
4. **Response Processing**: Parse to standardized `RepositoryEntry` objects
5. **Depth Filtering**: Apply client-side limits when APIs lack precise depth control

### API Strategy Classification

**Strategy A: Native Recursive (GitHub, GitLab, Azure)**
- Single HTTP request with recursive parameters
- Server-side tree traversal
- Performance: O(1) API calls

**Strategy B: Iterative Traversal (Bitbucket Server, Gitea)**
- Multiple HTTP requests per directory level  
- Client-side recursion management
- Performance: O(n) API calls where n = number of directories

**Strategy C: Limited Support (Bitbucket Cloud)**
- Single-level listing only
- Throws exceptions for depth > 1

### Provider Implementation Details

| Provider | Endpoint | Recursive Support | Performance |
|----------|----------|-------------------|-------------|
| GitHub | `/git/trees/{sha}?recursive=1` | Native | Single call |
| GitLab | `/repository/tree?recursive=true` | Native | Single call |
| Azure | `/items?recursionLevel=Full` | Native | Single call |
| Bitbucket Server | `/browse/{path}` | Manual iteration | Multiple calls |
| Gitea | `/contents/{path}` | Manual iteration | Multiple calls |
| Bitbucket Cloud | `/src/{commit}/{path}` | None | Unsupported |

### HTTP API Constraints

- **Rate Limiting**: 60-5000 requests/hour depending on provider and authentication
- **Response Size**: Controlled by `NXF_GIT_RESPONSE_MAX_LENGTH` environment variable
- **Timeouts**: 60-second connect timeout across all providers
- **Authentication**: Required for private repositories and higher rate limits

## Consequences

### Positive
- **Unified Interface**: Consistent API across all Git hosting providers
- **Performance Optimization**: Uses native recursive APIs where available
- **Graceful Degradation**: Falls back to iterative traversal when needed
- **Error Resilience**: Handles partial failures and API limitations

### Negative
- **Provider Inconsistency**: Performance varies widely between providers
- **API Rate Limits**: Providers that need several calls hit their rate limit sooner
- **Memory Usage**: Large directory structures loaded entirely into memory

### Neutral
- **Complexity**: Abstraction layer adds code complexity but improves maintainability
- **Testing**: Each provider implementation needs its own tests

## Implementation Notes

- Local Git repositories use JGit TreeWalk instead of an HTTP API
- Client-side depth filtering keeps behavior the same across providers
- Error handling varies by provider: some return empty lists, others throw exceptions
- Later work could add caching keyed on the commit SHA, and pagination

With this method in place, Nextflow can list repository contents on any supported Git host, using the fastest API each one offers.