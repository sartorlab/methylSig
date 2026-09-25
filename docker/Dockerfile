# Bioconductor devel with methylSig's dependencies preinstalled, for docker/check.sh.
# Rebuild when Bioconductor devel moves to a new version (after each release), and
# periodically within a cycle so dependencies match the Bioconductor builders:
#   docker/check.sh --rebuild
#
# The image comes from GitHub's registry; the Docker Hub copy is
# bioconductor/bioconductor_docker:devel.
ARG BIOC_IMAGE=ghcr.io/bioconductor/bioconductor:devel
FROM ${BIOC_IMAGE}

COPY docker/install_deps.R /opt/methylSig/install_deps.R
COPY DESCRIPTION /opt/methylSig/DESCRIPTION
RUN Rscript /opt/methylSig/install_deps.R /opt/methylSig/DESCRIPTION
