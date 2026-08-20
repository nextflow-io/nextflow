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

FROM node:24.4.1-bookworm-slim@sha256:36ae19f59c91f3303c7a648f07493fe14c4bd91320ac8d898416327bacf1bbfa

WORKDIR /opt/nf-agent-pi
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
COPY package.json package-lock.json ./
RUN npm ci --ignore-scripts --omit=dev && npm cache clean --force
COPY harness/runner.mjs ./runner.mjs
COPY --from=agent-rpc /agent-rpc /usr/local/bin/agent-rpc

# Canonical Nextflow executors own the container command/entrypoint so they can
# launch the task wrapper (and Fusion when enabled). The generated agent
# script explicitly execs /usr/local/bin/agent-rpc as its first process.
ENTRYPOINT []
