#!/usr/bin/env bash

set -euo pipefail

MAXWELL_USER_NAME="$1"
MAXWELL_HOST="$2"
PREFIX="$3"

BASE_REMOTE_DIR="/software/crystfel"
CURRENT_TARGET_DIR="$BASE_REMOTE_DIR/devel-$(date +%F-%H:%M:%S)"

if [ ! -d ~/.ssh ]; then \
    mkdir -p ~/.ssh; \
    chmod 700 ~/.ssh; \
fi;  \
wget https://wims.desy.de/system/ALL_afs/etc/ssh_known_hosts2 -O ~/.ssh/known_hosts

# To avoid people starting CrystFEL devel processes while those are
# being redeployed (leading to weird error messages), we copy first,
# and then atomically switch a symlink (see "mv -T" below) so we have
# minimal overlap time.

scp -r "$PREFIX/crystfel/devel" "${MAXWELL_USER_NAME}@${MAXWELL_HOST}:$CURRENT_TARGET_DIR"

ssh "${MAXWELL_USER_NAME}@${MAXWELL_HOST}" <<EOF
set -euo pipefail

rm -f "$BASE_REMOTE_DIR/crystfel-deployment-temp"
ln -s "$CURRENT_TARGET_DIR" "$BASE_REMOTE_DIR/crystfel-deployment-temp"
mv -T "$BASE_REMOTE_DIR/crystfel-deployment-temp" "$BASE_REMOTE_DIR/devel"

EOF
