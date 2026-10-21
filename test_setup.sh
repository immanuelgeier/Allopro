#!/bin/bash

echo "================================"
echo " TEST SETUP"
echo "================================"

# Pfad zum Quellcode (Standard: aktuelles Repository)
SRC_DIR="$(pwd)"

./install_mdrun.sh << EOF
2
$SRC_DIR
EOF

echo
echo "Bitte gib die Pfade zu deinen Testdaten an."

read -e -p "PDB-Datei: " PDB
read -e -p "Domains-Datei: " DOMAINS
read -e -p "DCD-Dateien (z.B. /pfad/zu/dcd_*): " DCD

echo
echo "Projekt wird eingerichtet..."

./config_project.sh << EOF
1
test_project
$PDB
$DOMAINS
$DCD
EOF

echo
echo "Testprojekt erfolgreich eingerichtet."
#make


