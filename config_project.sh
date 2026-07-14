#!/bin/bash

echo "=============================="
echo " Project structure builder"
echo "=============================="

PROJECT_DIR="project"
mkdir -p "$PROJECT_DIR"

#-----------------------------------
# STATES
#-----------------------------------
read -p "Number of states: " NSTATE

states=()

for ((i=1; i<=NSTATE; i++)); do
    read -p "State name $i: " s
    states+=("$s")
done

#-----------------------------------
# write state list
#-----------------------------------

#====================================================
# STATE LIST (INPUT ROOT)
#====================================================
STATE_LIST_FILE="$PROJECT_DIR/input/state_list.dat"

mkdir -p "$PROJECT_DIR/input"

printf "%s\n" "${states[@]}" > "$STATE_LIST_FILE"
echo "END" >> "$STATE_LIST_FILE"
#-----------------------------------
# LOOP STATES
#-----------------------------------
for state in "${states[@]}"; do

    echo ""
    echo "=============================="
    echo " STATE: $state"
    echo "=============================="

    #====================================================
    # INPUT / OUTPUT ROOT (NEW STRUCTURE)
    #====================================================
    INPUT_STATE_DIR="$PROJECT_DIR/input/$state"
    OUTPUT_STATE_DIR="$PROJECT_DIR/output/$state"

    mkdir -p "$INPUT_STATE_DIR"
    mkdir -p "$OUTPUT_STATE_DIR"
    #====================================================
    # TOPOLOGY FILES
    #====================================================
    echo ""
    echo "Enter reference topology for state '$state':"
    
    read -p "PDB path: " PDB_PATH
    read -p "Domain file: " DOMAIN_PATH
    
    STATE_CONFIG="$INPUT_STATE_DIR/state_config.dat"
    
    {
        echo "pdb $PDB_PATH"
        echo "domains $DOMAIN_PATH"
        echo "END"
    } > "$STATE_CONFIG"
    
    echo "Created state config:"
    echo "  $STATE_CONFIG"
    echo ""
    #-----------------------------------
    # REPLICATES
    #-----------------------------------
    echo ""
    echo "Enter replicate directory pattern."
    echo "Example:"
    echo "  /data/apo/dcd_*"
    echo ""

    read -p "Pattern: " REP_PATTERN

    REPLIST_FILE="$INPUT_STATE_DIR/dcd_dir_list.dat"
    > "$REPLIST_FILE"

    REP_DIRS=( $REP_PATTERN )

    NREP=0

    for REP_DIR in "${REP_DIRS[@]}"; do

        if [ ! -d "$REP_DIR" ]; then
            continue
        fi

        NREP=$((NREP+1))

        echo "$REP_DIR" >> "$REPLIST_FILE"

        REP_NAME="replicate_${NREP}"

        mkdir -p "$INPUT_STATE_DIR/$REP_NAME"

        #====================================================
        # OUTPUT PER REPLICATE
        #====================================================
        OUT_REP_DIR="$OUTPUT_STATE_DIR/$REP_NAME"
        mkdir -p "$OUT_REP_DIR/logs"

        FILELIST="$INPUT_STATE_DIR/$REP_NAME/dcd_file_list.dat"

        echo "Scanning: $REP_DIR"

        DCD_COUNT=$(find "$REP_DIR" -maxdepth 1 -type f -name "*.dcd" | wc -l)

        if [ "$DCD_COUNT" -eq 0 ]; then
            echo "WARNING: no dcd files found"
            > "$FILELIST"
        else
            find "$REP_DIR" -maxdepth 1 -type f -name "*.dcd" | sort > "$FILELIST"
            echo "Found $DCD_COUNT dcd files"
        fi

        echo "END" >> "$FILELIST"

        echo "Created:"
        echo "  $FILELIST"
        echo "  OUTPUT: $OUT_REP_DIR"
        echo ""

    done

    echo "END" >> "$REPLIST_FILE"

    echo ""
    echo "State '$state': found $NREP replicate directories"

done

echo ""
echo "=============================="
echo " DONE"
echo "=============================="
echo ""
echo "Created project structure:"
echo ""
echo "project/"
echo "  input/"
echo "  output/"
echo "  state_list.dat"
echo "  input/<state>/..."
echo "  output/<state>/..."
echo ""