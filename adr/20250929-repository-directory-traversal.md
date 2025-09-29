 # ADR: Repository Directory Traversal API

**Date**: 2025-09-29  
**Status**: Accepted  
**Context**: Need for standardized directory listing across Git hosting providers

## Decision

Introduce a `listDirectory(String path, int depth)` method to the `RepositoryProvider` abstraction to enable unified directory traversal across different Git hosting platforms.

## Context

Nextflow requires the ability to explore repository directory structures across multiple Git hosting providers (GitHub, GitLab, Bitbucket, Azure DevOps, Gitea) without full repository clones. Each provider has different API capabilities and constraints for directory listing operations.

## Technical Implementation

### Core Algorithm

All providers follow a consistent pattern:
1. **Path Resolution**: Normalize path to provider API format
2. **Strategy Selection**: Choose recursive vs iterative approach based on API capabilities
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
| GitHub | `/git/trees/{sha}?recursive=1` | Native | Optimal |
| GitLab | `/repository/tree?recursive=true` | Native | Optimal |
| Azure | `/items?recursionLevel=Full` | Native | Optimal |
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
- **Provider Inconsistency**: Performance varies significantly between providers
- **API Rate Limits**: Multiple calls required for some providers may hit limits faster
- **Memory Usage**: Large directory structures loaded entirely into memory

### Neutral
- **Complexity**: Abstraction layer adds code complexity but improves maintainability
- **Testing**: Comprehensive test coverage required for each provider implementation

## Implementation Notes

- Local Git repositories use JGit TreeWalk for optimal performance
- Client-side depth filtering ensures consistent behavior across providers
- Error handling varies by provider: some return empty lists, others throw exceptions
- Future enhancements could include caching based on commit SHA and pagination support

This decision enables Nextflow to efficiently explore repository structures regardless of the underlying Git hosting platform, with automatic optimization based on each provider's API capabilities.