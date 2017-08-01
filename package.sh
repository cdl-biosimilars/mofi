#!/bin/bash

set -e

if [[ "$(uname -a)" =~ "MINGW" ]]; then
    TARGET=windows
    VENV_BIN=venv/Scripts
else
    TARGET=linux
    VENV_BIN=venv/bin
fi

if python3 --version >/dev/null 2>&1; then
    NATIVE_PYTHON=python3
else
    NATIVE_PYTHON=python
fi

rm -rf venv
if [[ "$($NATIVE_PYTHON --version)" =~ "Continuum Analytics, Inc." ]]; then
    echo "Creating venv with Anaconda"
    $NATIVE_PYTHON -m venv --without-pip venv
    wget https://bootstrap.pypa.io/get-pip.py -O /tmp/get-pip.py
    $VENV_BIN/python /tmp/get-pip.py
else
    echo "Creating venv"
    $NATIVE_PYTHON -m venv venv
fi

$VENV_BIN/pip3 install setuptools wheel cx_Freeze --upgrade
$VENV_BIN/pip3 install .

if [[ $TARGET == windows ]]; then
    echo "Building for Windows"
    $VENV_BIN/python setup_cx.py build
    cd build
    BUILD_DIR="$(find . -type d -name "exe.win*" -print -quit)"
    if [[ -z "BUILD_DIR" ]]; then
        echo "Could not determine build directory"
        exit 1
    fi
    TARGET_DIR="mofi-windows"
    echo "Copying build files"
    rm -rf "$TARGET_DIR"
    mkdir -p "${TARGET_DIR}/bin"
    cp -r "${BUILD_DIR}"/* "${TARGET_DIR}/bin"
    cp ../package-files/win-standalone/* "${TARGET_DIR}"
    echo "Please zip build/mofi-windows now."
else
    echo "Building for Linux"
    $VENV_BIN/python setup_cx.py build
    cd build
    BUILD_DIR="$(find . -type d -name "exe.linux*" -print -quit)"
    if [[ -z "BUILD_DIR" ]]; then
        echo "Could not determine build directory"
        exit 1
    fi
    TARGET_DIR="mofi-linux"
    echo "Copying build files"
    rm -rf "$TARGET_DIR"
    mkdir -p "${TARGET_DIR}/bin"
    cp -r "${BUILD_DIR}"/* "${TARGET_DIR}/bin"
    cp ../package-files/lin-standalone/* "${TARGET_DIR}"
    ln -s bin/data "${TARGET_DIR}/data"
    ln -s bin/config "${TARGET_DIR}/config"
    rm -f mofi-linux.tar.gz
    echo "Creating archive"
    tar -czf mofi-linux.tar.gz "${TARGET_DIR}"/*
fi

