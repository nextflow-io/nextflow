#!/bin/bash
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
# Build and publish the `pi` agent runner image for linux/amd64 and linux/arm64.
#
# nf-agent-pi ships no runtime of its own, so this image IS the distribution unit: an
# agent selecting the `pi` runner runs in it. The Nextflow release publishes it, through this
# script: release.sh step 1 runs `make release-agent-image`, which runs the Gradle task
# :plugins:nf-agent-pi:releaseImageIfNotExists, which runs `push` below. That is a no-op when
# the tag is already published, so the release is safe to re-run. An ordinary build never
# touches docker. Running this script by hand remains the escape hatch, and is how a tag gets
# published outside a release.
#
#   plugins/nf-agent-pi/build-image.sh build              # both arches, publishes nothing
#   plugins/nf-agent-pi/build-image.sh build -l           # host arch only, into local docker
#   plugins/nf-agent-pi/build-image.sh push -r <registry> # build, push, verify the manifest
#   plugins/nf-agent-pi/build-image.sh ref                # print the reference; no docker needed
#
# A single-arch image is the failure this script exists to prevent: pulling one on the other
# architecture fails with `no matching manifest for linux/amd64 in the manifest list entries`.
set -u

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PLATFORMS=linux/amd64,linux/arm64
BUILDER=nf-agent-pi
OUT_DIR=$SCRIPT_DIR/build/image

# THE image coordinate. Everything else derives from it: the reference this script builds and
# pushes, the coordinate the examples pin - see examples/agents/pi-runner.config - and the
# "<DEFAULT_REGISTRY>/<IMAGE>:<VERSION>" string the plugin jar embeds and asks for at run time
# (build.gradle reads these two lines by anchored match; keep them plain assignments). This is
# the same registry and namespace the Nextflow release already publishes `nextflow/nextflow` to
# (see docker/Makefile), so the release job's existing login covers it and no new credential
# exists for the runner image. Overridden for a private build by NF_AGENT_PI_REGISTRY or -r.
DEFAULT_REGISTRY=public.cr.seqera.io/nextflow
REGISTRY=${NF_AGENT_PI_REGISTRY:-$DEFAULT_REGISTRY}
IMAGE=nf-agent-pi
TAG=""
LOAD=0              # build -l: single-arch into the local image store
KEEP=0              # build -k: do not delete the verification archive
FORCE=0             # push -f: overwrite a tag that is already published

usage() {
  cat <<'TXT'
Build and publish the multi-arch nf-agent-pi runner image.

  usage: plugins/nf-agent-pi/build-image.sh <command> [options]

commands:
  build     build every platform locally; publishes nothing
  push      build every platform, push it, then verify the pushed manifest
  ref       print the fully resolved image reference and exit; needs no docker

options (every command):
  -r  registry/namespace; default public.cr.seqera.io/nextflow, the coordinate the
      release publishes. Overridable here or with NF_AGENT_PI_REGISTRY; pass -r '' to
      build under a bare local name (which cannot be pushed)
  -t  image tag; default the plugin VERSION file
  -i  image name; default nf-agent-pi
  -P  platforms; default linux/amd64,linux/arm64
  -h  this help

build only:
  -l  build just this host's architecture and load it into the local docker image
      store, for a local test run; a multi-platform result cannot be loaded
  -k  keep the OCI archive the build writes; it is deleted by default

push only:
  -f  build and push even when the tag is already published, overwriting it

`build` writes an OCI archive under plugins/nf-agent-pi/build/image/, which proves the
build works and is then deleted again - it weighs ~200 MB and nothing consumes it. Pass
-k to inspect it. Layer cache stays in the dedicated buildx builder either way, so
re-runs are fast.

Tag releases with the plugin VERSION, which is the default, and keep the tag immutable,
so an agent's runtime is pinned as reproducibly as its model. `push` therefore publishes
only if the tag does not exist yet: an already published tag is a no-op that exits 0, so
a failed release can be re-run, and overwriting one takes an explicit -f. "Already
published" means its manifest carries every requested platform - a tag that exists but is
single-arch is a failure, not a skip.
TXT
  exit "${1:-0}"
}

die() { echo "ERROR: $*" >&2; exit 1; }

# ---------------------------------------------------------------- command line

CMD=${1:-}
case $CMD in
  build|push|ref) shift ;;
  help|-h|--help) usage 0 ;;
  '')           echo "Missing command" >&2; usage 2 ;;
  *)            echo "Unknown command: $CMD" >&2; usage 2 ;;
esac

while getopts ':r:t:i:P:lkfh' opt; do
  case $opt in
    r) REGISTRY=$OPTARG ;;
    t) TAG=$OPTARG ;;
    i) IMAGE=$OPTARG ;;
    P) PLATFORMS=$OPTARG ;;
    l) LOAD=1 ;;
    k) KEEP=1 ;;
    f) FORCE=1 ;;
    h) usage 0 ;;
    :) echo "Missing argument for -$OPTARG" >&2; usage 2 ;;
    ?) echo "Unknown option -$OPTARG" >&2; usage 2 ;;
  esac
done
shift $((OPTIND - 1))
(( $# )) && { echo "Unexpected argument: $1" >&2; usage 2; }

if [[ $CMD == push ]]; then
  (( LOAD )) && die "-l builds a single architecture and cannot be pushed; use build -l"
  (( KEEP )) && die "-k applies to build, which writes an archive; push writes none"
elif (( FORCE )); then
  die "-f overwrites an already published tag and applies to push, which publishes; $CMD does not"
fi

# A stale DOCKER_DEFAULT_PLATFORM silently narrows what buildx produces, which is exactly
# the single-arch image this script exists to avoid. --platform below is authoritative.
# The note goes to stderr, not stdout: `ref` writes ONE machine-readable line that callers
# capture with $(...), and a diagnostic mixed into it would be read as part of the reference.
if [[ -n ${DOCKER_DEFAULT_PLATFORM:-} ]]; then
  echo "note: ignoring DOCKER_DEFAULT_PLATFORM=$DOCKER_DEFAULT_PLATFORM for this build" >&2
  unset DOCKER_DEFAULT_PLATFORM
fi

[[ -f $SCRIPT_DIR/VERSION ]] || die "no VERSION file beside $0"
[[ -n $TAG ]] || TAG=$(tr -d '[:space:]' < "$SCRIPT_DIR/VERSION")
[[ -n $TAG ]] || die "VERSION is empty and no -t tag was given"

if [[ -n $REGISTRY ]]; then
  REF="${REGISTRY%/}/$IMAGE:$TAG"
else
  # only reachable via an explicit `-r ''`, which asks for a bare, unpushable name
  [[ $CMD == push ]] && die "push needs a registry, but -r was given an empty value"
  REF="$IMAGE:$TAG"
fi

# The one place the coordinate is composed. `ref` is the read side of it: the plugin build embeds
# this string in the jar (see build.gradle), so the image the jar asks for is the image this
# script pushes. It exits BEFORE the docker checks below on purpose - resolving a name needs no
# docker, and the build that cross-checks the embedded coordinate must not require one.
if [[ $CMD == ref ]]; then
  echo "$REF"
  exit 0
fi

command -v docker >/dev/null || die "docker is not on PATH"
docker buildx version >/dev/null 2>&1 || die "docker buildx is not available"

# ------------------------------------------------------------------- commands

ensure_builder() {
  # The default `docker` driver cannot emit a multi-platform manifest - it would silently
  # build one architecture. A docker-container builder can, so use a dedicated one.
  docker buildx inspect "$BUILDER" >/dev/null 2>&1 && return 0
  echo "==> creating buildx builder '$BUILDER' (docker-container driver)"
  docker buildx create --name "$BUILDER" --driver docker-container --bootstrap >/dev/null \
    || die "could not create the '$BUILDER' builder"
}

cmd_build() {
  local platforms=$PLATFORMS archive="" out=()

  if (( LOAD )); then
    # --load takes a single image; ask for the host's platform explicitly so the result is
    # runnable here whatever DOCKER_DEFAULT_PLATFORM said.
    local host_arch
    host_arch=$(docker version --format '{{.Server.Arch}}' 2>/dev/null) \
      || die "docker daemon is not reachable"
    platforms=linux/$host_arch
    out=( --load )
    echo "==> building $platforms and loading $REF into the local image store"
  else
    # A previous run's archive is never reused - drop it before writing another one, so
    # repeated verification runs do not accumulate.
    rm -rf "$OUT_DIR"
    mkdir -p "$OUT_DIR" || die "cannot create $OUT_DIR"
    archive=$OUT_DIR/$IMAGE-$TAG.oci.tar
    out=( --output "type=oci,dest=$archive" )
    echo "==> building $platforms to $archive (not published)"
  fi

  docker buildx build --builder "$BUILDER" --pull -t "$REF" \
    --platform "$platforms" "${out[@]}" "$SCRIPT_DIR" || die "build failed"

  if (( LOAD )); then
    echo "==> OK - $REF is in the local image store"
    echo "    point the examples at it:  agent.container = '$REF'"
    return 0
  fi

  echo "==> OK - $REF built for $platforms"
  # The archive is a by-product of proving the build works, not a deliverable: the exit
  # code above already carries that answer. Keep it only when the point was to inspect it.
  if (( KEEP )); then
    echo "    archive kept: $archive ($(du -h "$archive" | cut -f1))"
  else
    rm -rf "$OUT_DIR"
    echo "    archive removed; -k keeps it, or use the push command to publish"
  fi
}

# Assert that the manifest published under $REF carries every requested platform, printing it.
# Returns non-zero, quietly, when the manifest cannot be read at all - "absent", "unreachable" and
# "unauthenticated" are indistinguishable here, and the caller decides what that means. A manifest
# that IS readable but incomplete is fatal wherever it is found: publishing one under the release
# coordinate is the failure this script exists to prevent, so it must not be reachable by pushing
# one, nor by SKIPPING one that a half-finished push - or a hand-run `push -P linux/amd64` - left
# behind. Hence the same check on both paths.
verify_manifest() {
  local manifest missing=() p wanted
  manifest=$(docker buildx imagetools inspect "$REF" 2>/dev/null) || return 1
  IFS=',' read -ra wanted <<< "$PLATFORMS"
  for p in "${wanted[@]}"; do
    # arm64 is reported as linux/arm64 or linux/arm64/v8 depending on the base image
    grep -Eq "Platform:[[:space:]]+${p}(/v[0-9]+)?$" <<< "$manifest" || missing+=( "$p" )
  done
  if (( ${#missing[@]} )); then
    echo "$manifest" >&2
    die "pushed manifest is missing: ${missing[*]} - $REF is published but incomplete; republish it with 'push -f', or bump VERSION, before releasing against it"
  fi
  echo "$manifest" | grep -E '^(Name|Digest):|Platform:'
}

cmd_push() {
  echo "==> building $PLATFORMS and pushing $REF"
  # naming the registry, because an authentication failure surfaces here as a build failure and
  # this runs inside a release, where the operator sees only this line
  docker buildx build --builder "$BUILDER" --pull -t "$REF" \
    --platform "$PLATFORMS" --push "$SCRIPT_DIR" \
    || die "build or push of $REF failed - if the registry refused it, log in to ${REGISTRY%%/*} first (the release job does that with SEQERA_PUBLIC_CR_USERNAME/SEQERA_PUBLIC_CR_PASSWORD)"

  echo "==> verifying the pushed manifest"
  verify_manifest || die "cannot inspect $REF after push"
  echo "==> OK - $REF carries every requested platform"
  echo "    point the examples at it:  agent.container = '$REF'"
}

if [[ $CMD == push ]] && (( ! FORCE )); then
  # The tag is a runtime pin, so it is immutable: an existing tag is success, not a conflict.
  # Mirrors releasePluginToRegistryIfNotExists, and makes a re-run of a failed release safe.
  # It is success only if the published manifest is COMPLETE, though - verify_manifest dies on a
  # readable but single-arch one rather than skipping it, so a re-run after a push that published
  # and then failed verification re-checks instead of inheriting the earlier run's mistake.
  # An unreadable manifest ("absent", "unreachable" and "unauthenticated" are all non-zero and
  # indistinguishable) falls through to build-and-push and fails at the push, where the error is
  # accurate. The push, not this probe, is the authority on what is published.
  if verify_manifest >/dev/null; then
    echo "==> $REF is already published with every requested platform - nothing to do (-f rebuilds and overwrites it)"
    exit 0
  fi
fi

ensure_builder
"cmd_$CMD"
