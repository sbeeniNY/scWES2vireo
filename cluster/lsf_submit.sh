#!/usr/bin/env bash
set -euo pipefail
# Wrapper: submit via bsub and print only the numeric job ID to stdout.
# cluster-generic plugin expects a single job ID on stdout.
bsub "$@" 2>&1 | sed -n 's/.*Job <\([0-9]*\)>.*/\1/p'
