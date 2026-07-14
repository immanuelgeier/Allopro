#!/bin/bash

echo "======================================"
echo " MDGraph Workspace DECONFIGURE"
echo "======================================"
echo ""

RUN_DIR=$(pwd)

echo "This will remove MDGraph workspace setup in:"
echo "  $RUN_DIR"
echo ""
echo "It will NOT delete external data (PDB/DCD files)."
echo ""

read -p "Are you sure? (yes/no): " CONFIRM

if [ "$CONFIRM" != "yes" ]; then
    echo "Aborted."
    exit 0
fi

echo ""
echo "Cleaning workspace..."

#--------------------------------------
# remove generated runtime artifacts
#--------------------------------------

rm -f run.sh
rm -rf build/mod
rm -f run/mdgraph

# remove build symlink (IMPORTANT)
if [ -L build/src ]; then
    rm -f build/src
    echo "Removed build/src symlink"
fi

# optionally remove build dir if empty
if [ -d build ] && [ -z "$(ls -A build 2>/dev/null)" ]; then
    rmdir build
fi

#--------------------------------------
# optional project cleanup (SAFE MODE)
#--------------------------------------
echo ""
read -p "Also delete PROJECT directory? (yes/no): " DEL_PROJECT

if [ "$DEL_PROJECT" == "yes" ]; then
    if [ -d project ]; then
        rm -rf project
        echo "Removed project/"
    fi
else
    echo "Keeping project/"
fi

#--------------------------------------
# optional src_copy cleanup (MODE 3 artifact)
#--------------------------------------
if [ -d src_copy ]; then
    read -p "Delete src_copy? (yes/no): " DEL_SRC

    if [ "$DEL_SRC" == "yes" ]; then
        rm -rf src_copy
        echo "Removed src_copy/"
    fi
fi

#--------------------------------------
# optional environment hints
#--------------------------------------
echo ""
echo "Unsetting environment hints (for current shell only)"

unset MDGRAPH_PROJECT
unset MDGRAPH_SRC

echo ""
echo "======================================"
echo " DECONFIGURE COMPLETE"
echo "======================================"
echo ""
echo "Workspace reset done."
echo "You can safely run install_mdrun.sh again."
echo ""