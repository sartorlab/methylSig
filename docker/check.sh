#!/usr/bin/env bash
# Run R CMD build, R CMD check, and BiocCheck on the committed HEAD in a
# Bioconductor devel container. Uncommitted changes are not checked.
#
# Caches (e.g. BiocCheck's) persist between runs in the methylsig-cache Docker
# volume.
#
# Usage: docker/check.sh [--rebuild]
#   --rebuild  Rebuild the image, pulling the latest Bioconductor devel image.
#
# Results are left in a temporary directory, printed at the end.
set -euo pipefail

repo=$(git rev-parse --show-toplevel)
image=methylsig-dev:devel

if [[ "${1:-}" == "--rebuild" ]] || ! docker image inspect "$image" > /dev/null 2>&1; then
    docker build --pull -f "$repo/docker/Dockerfile" -t "$image" "$repo"
fi

work=$(mktemp -d)
mkdir "$work/src"
git -C "$repo" archive HEAD | tar -x -C "$work/src"

echo "Checking $(git -C "$repo" rev-parse --short HEAD) in $work"

docker run --rm -v "$work:/work" -v methylsig-cache:/root/.cache -w /work \
    "$image" bash -c '
    set -e
    Rscript /opt/methylSig/install_deps.R src/DESCRIPTION
    R CMD build src
    tarball=$(ls methylSig_*.tar.gz)
    status=0
    R CMD check --no-manual "$tarball" || status=$?
    Rscript -e "BiocCheck::BiocCheck(\"$tarball\")" || status=$?
    exit $status
'

echo "Results: $work/methylSig.Rcheck and $work/methylSig.BiocCheck"
