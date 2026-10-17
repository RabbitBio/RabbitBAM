#!/usr/bin/env bash
set -Eeuo pipefail

# Use the already installed private glibc runtime; do not change system libraries.
CONNECT=/home/user_home/wzs/.aspera/connect
RUNTIME=/home/user_home/wzs/software/aspera-glibc-2.35/root/lib/x86_64-linux-gnu
export ASPERA_SCP_LICENSE="$(< "$CONNECT/etc/aspera-license")"
exec "$RUNTIME/ld-linux-x86-64.so.2" \
    --library-path "$RUNTIME:$CONNECT/lib:/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu" \
    "$CONNECT/bin/ascp" "$@"
