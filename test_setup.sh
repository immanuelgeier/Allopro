#!/bin/bash

echo "================================"
echo " TEST SETUP"
echo "================================"

./install_mdrun.sh << EOF
2
/home/pbuser/Schreibtisch/mdgraph/src
EOF


# =========================================
# CONFIG PROJECT
# =========================================

PDB=/home/pbuser/Schreibtisch/MD_simulations/monomer_prod/program_files/protein.pdb
DOMAINS=/home/pbuser/Schreibtisch/3lnx_domains.txt
DCD=/media/pbuser/T7/MD_3lnx/dcd_*

./config_project.sh << EOF
1
3lnx
$PDB
$DOMAINS
$DCD
EOF

#make


