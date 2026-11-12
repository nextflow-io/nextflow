# Runner image for the `pi` agent runner. This image IS the distribution unit of the
# nf-agent-pi runtime: the plugin jar carries no proxy binary and no node_modules, so an
# agent selecting the `pi` runner must run in this image (`agent.container`).
#
# Build context is this directory. The Nextflow release publishes this image: release.sh step 1
# runs `make release-agent-image`, i.e. the script beside this file, which sets up a
# multi-platform builder and verifies the pushed manifest. That script is also the escape hatch
# for publishing a tag by hand:
#
#   plugins/nf-agent-pi/build-image.sh push
#
# The tag is the nf-agent-pi plugin version (see the VERSION file in this directory) so an
# agent's runtime is pinned as reproducibly as its model, and it is immutable: `push` is a
# no-op when the tag already exists, so bump VERSION to publish new image content.

# The proxy is cross-compiled ON the build host rather than emulated: with CGO_ENABLED=0 Go
# targets any GOOS/GOARCH natively, so pinning this stage to $BUILDPLATFORM keeps the amd64
# leg of a multi-arch build off QEMU. TARGETOS/TARGETARCH are supplied by BuildKit.
#
# Digest-pinned like the runtime base below, so the proxy the release publishes is built by a
# toolchain this repo names rather than by whatever `1.25-bookworm` points at that day. This is
# the OCI *index* digest (`docker buildx imagetools inspect golang:1.25-bookworm`), not a
# per-architecture manifest digest: a manifest digest would pin the stage to one build host.
FROM --platform=$BUILDPLATFORM golang:1.25-bookworm@sha256:3b4a11519ad929d1e1d261a12cff056f0c85b735253d7d861346b9c6f8b36437 AS agent-rpc
ARG TARGETOS
ARG TARGETARCH

WORKDIR /src
COPY agent-rpc/go.mod agent-rpc/go.sum ./
RUN go mod download
COPY agent-rpc/main.go ./
RUN CGO_ENABLED=0 GOOS=${TARGETOS} GOARCH=${TARGETARCH} \
    go build -trimpath -ldflags='-s -w' -o /agent-rpc .

# The harness is compiled into a single executable that embeds the Bun runtime, so the runtime
# image needs neither Node nor node_modules. Bun cross-compiles for any target, so like the
# proxy this stage runs on $BUILDPLATFORM. `bun install` reads package-lock.json directly.
FROM --platform=$BUILDPLATFORM oven/bun:1.4.4-debian@sha256:93d5032fe06d163427822dc1598f420363759bb18155c9cc90eb1be05a64a19e AS harness
ARG TARGETARCH

WORKDIR /src
COPY package.json package-lock.json ./
RUN bun install --frozen-lockfile --production --ignore-scripts
COPY harness/runner.mjs harness/entry.mjs ./harness/
RUN case "${TARGETARCH}" in \
        amd64) target=bun-linux-x64 ;; \
        arm64) target=bun-linux-arm64 ;; \
        *) echo "unsupported TARGETARCH: ${TARGETARCH}" >&2; exit 1 ;; \
    esac \
    && bun build harness/entry.mjs --compile --minify --target=${target} --outfile=/harness/runner
# Pi passes every image a `read` returns through photon, and a compiled binary cannot embed its
# wasm, so ship it next to the binary, where pi looks for it (`getFallbackWasmPaths`). Without it,
# every image is dropped, not only ones that need resizing.
RUN wasm=$(find node_modules -name photon_rs_bg.wasm) \
    && if [ "$(printf '%s\n' "$wasm" | grep -c .)" -ne 1 ]; then \
        echo "expected exactly one photon_rs_bg.wasm, found: ${wasm:-none}" >&2; exit 1; \
    fi \
    && cp "$wasm" /harness/photon_rs_bg.wasm

FROM debian:bookworm-slim@sha256:99c0b25c46d0214d5a2091092e7908907dcaecd6a962cde0ac700acee222b6f3

RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
COPY --from=harness /harness/runner /opt/nf-agent-pi/runner
COPY --from=harness /harness/photon_rs_bg.wasm /opt/nf-agent-pi/photon_rs_bg.wasm
COPY --from=agent-rpc /agent-rpc /usr/local/bin/agent-rpc

# Canonical Nextflow executors own the container command/entrypoint so they can
# launch the task wrapper (and Fusion when enabled). The generated agent
# script explicitly execs /usr/local/bin/agent-rpc as its first process.
ENTRYPOINT []
