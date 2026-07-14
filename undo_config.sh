#!/bin/bash

echo "=============================="
echo " Project structure UNDO"
echo "=============================="
echo ""

PROJECT_DIR="project"

echo "This will remove ONLY generated project structure:"
echo "  - project/state_list.dat"
echo "  - project/input/"
echo "  - project/output/"
echo ""
echo "It will NOT touch:"
echo "  - PDB files"
echo "  - DCD files"
echo "  - external directories (/media, /home, etc.)"
echo ""

read -p "Are you sure? (yes/no): " CONFIRM

if [ "$CONFIRM" != "yes" ]; then
echo "Aborted."
exit 0
fi

echo ""
echo "Removing generated project structure..."

#-----------------------------------

# safety check

#-----------------------------------
if [ ! -d "$PROJECT_DIR" ]; then
echo "Nothing to do: project/ does not exist"
exit 0
fi

#-----------------------------------

# remove generated files

#-----------------------------------
if [ -f "$PROJECT_DIR/state_list.dat" ]; then
echo "Removing state list..."
rm -f "$PROJECT_DIR/state_list.dat"
fi

if [ -d "$PROJECT_DIR/input" ]; then
echo "Removing input structure..."
rm -rf "$PROJECT_DIR/input"
fi

if [ -d "$PROJECT_DIR/output" ]; then
echo "Removing output structure..."
rm -rf "$PROJECT_DIR/output"
fi

#-----------------------------------

# cleanup project dir if empty

#-----------------------------------
if [ -d "$PROJECT_DIR" ] && [ -z "$(ls -A "$PROJECT_DIR" 2>/dev/null)" ]; then
rmdir "$PROJECT_DIR"
echo "Removed empty project/ directory"
fi

echo ""
echo "=============================="
echo " UNDO COMPLETE"
echo "=============================="
echo ""
echo "All generated project structure removed."
echo "External trajectory/data files were untouched."
echo ""
