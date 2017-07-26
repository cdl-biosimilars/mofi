#!/bin/sh

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
BIN_DIR="${SCRIPT_DIR}/bin"

export LD_LIBRARY_PATH="$BIN_DIR"
"${BIN_DIR}/run"
