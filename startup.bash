#!/bin/bash

if [ -n "${HOMEPATH}" ]; then
  echo "ERROR: LhARAOptcs environment appears to already be configured"
  return
fi

export HOMEPATH=$PWD
export SCRIPT_DIR=$(dirname $(realpath "${BASH_SOURCE[0]}"))
export LhARAOpticsPATH="${SCRIPT_DIR}"
export PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}${LhARAOpticsPATH}/01-Code"
export REPORTPATH="${HOMEPATH}/99-Scratch"
mkdir -p "${REPORTPATH}"

if [ "$1" == "-d" ]; then
  echo "LhARAOptcs tool run from directory: ${HOMEPATH}"
  echo "LhARAOptcs path set: ${LhARAOpticsPATH}"
  echo "Python path set: ${PYTHONPATH}"
  echo "Reports path set: ${REPORTPATH}"
fi
