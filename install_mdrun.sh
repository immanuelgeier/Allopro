#!/bin/bash

echo "======================================"
echo " MDGraph Workspace Setup (robust mode)"
echo "======================================"
echo ""

RUN_DIR=$(pwd)

echo "Workspace root:"
echo "  $RUN_DIR"
echo ""

mkdir -p project
mkdir -p build
mkdir -p run

PROJECT_DIR="$RUN_DIR/project"

#--------------------------------------
# MODE SELECTION
#--------------------------------------
echo "Select setup mode:"
echo "  [1] Auto-detect SRC_ROOT"
echo "  [2] Provide SRC_ROOT manually"
echo "  [3] Copy SRC_ROOT into workspace (snapshot mode)"
echo ""

read -p "Choice [1-3]: " MODE

SRC_ROOT=""

#--------------------------------------
# MODE 1: auto detect
#--------------------------------------
if [ "$MODE" == "1" ]; then

    CANDIDATES=(
        "./src"
        "../src"
        "$RUN_DIR/src"
    )

    for c in "${CANDIDATES[@]}"; do
        if [ -d "$c" ]; then
            SRC_ROOT="$c"
            break
        fi
    done

    if [ -z "$SRC_ROOT" ]; then
        echo "ERROR: no src directory found automatically"
        echo "Hint: use mode 2 and set path manually"
        exit 1
    fi

    echo "Using detected SRC_ROOT: $SRC_ROOT"

#--------------------------------------
# MODE 2: manual path
#--------------------------------------
elif [ "$MODE" == "2" ]; then

    read -p "Path to SRC_ROOT: " SRC_ROOT

    if [ ! -d "$SRC_ROOT" ]; then
        echo "ERROR: path does not exist"
        exit 1
    fi

#--------------------------------------
# MODE 3: copy into workspace
#--------------------------------------
elif [ "$MODE" == "3" ]; then

    read -p "Path to SRC_ROOT to copy: " SRC_ROOT

    if [ ! -d "$SRC_ROOT" ]; then
        echo "ERROR: invalid src path"
        exit 1
    fi

    echo "Copying SRC_ROOT into workspace..."

    rm -rf "$RUN_DIR/src_copy"
    mkdir -p "$RUN_DIR/src_copy"

    cp -r "$SRC_ROOT"/* "$RUN_DIR/src_copy/" 2>/dev/null

    SRC_ROOT="$RUN_DIR/src_copy"

else
    echo "ERROR: invalid selection"
    exit 1
fi

#--------------------------------------
# NORMALIZE PATH (important fix)
#--------------------------------------
SRC_ROOT=$(cd "$SRC_ROOT" && pwd)

echo ""
echo "Final SRC_ROOT resolved to:"
echo "  $SRC_ROOT"

#--------------------------------------
# BASIC VALIDATION
#--------------------------------------
echo ""
echo "Checking SRC_ROOT..."

F90_COUNT=$(find "$SRC_ROOT" -name "*.f90" 2>/dev/null | wc -l)

if [ "$F90_COUNT" -eq 0 ]; then
    echo "ERROR: no Fortran (.f90) files found in SRC_ROOT"
    exit 1
fi

echo "Found $F90_COUNT Fortran files"

# optional structure hint
if [ ! -d "$SRC_ROOT/core" ] && [ ! -d "$SRC_ROOT/io" ] && [ ! -d "$SRC_ROOT/features" ]; then
    echo "NOTE: no standard layout detected (core/io/features)"
    echo "      -> OK for custom project layout"
fi

# main detection
MAIN_COUNT=$(find "$SRC_ROOT" -name "main.f90" 2>/dev/null | wc -l)

if [ "$MAIN_COUNT" -gt 0 ]; then
    echo "Found $MAIN_COUNT main.f90 file(s)"
else
    echo "NOTE: no main.f90 found (library-style build)"
fi

#--------------------------------------
# LINK SRC INTO BUILD
#--------------------------------------
rm -f build/src

ln -sfn "$SRC_ROOT" build/src

echo "Linked SRC_ROOT -> build/src"

#--------------------------------------
# RUN SCRIPT
#--------------------------------------
cat > run.sh << EOF
#!/bin/bash

export MDGRAPH_PROJECT="$PROJECT_DIR"
export MDGRAPH_SRC="$SRC_ROOT"

echo "Project dir: \$MDGRAPH_PROJECT"
echo "Source dir : \$MDGRAPH_SRC"

./run/mdgraph
EOF

chmod +x run.sh

#--------------------------------------
# SUMMARY
#--------------------------------------
echo ""
echo "======================================"
echo " Workspace ready"
echo ""
echo "SRC_ROOT:    $SRC_ROOT"
echo "PROJECT_DIR: $PROJECT_DIR"
echo ""
echo "Run with: ./run.sh"
echo "======================================"