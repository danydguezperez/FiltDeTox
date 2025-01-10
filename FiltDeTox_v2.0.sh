#!/bin/bash

# FiltDeTox Pipeline Script with directory navigation

# Step 1: hhmer_tx_VenomZone Module
echo "Running hhmer_tx_VenomZone Module..."
cd hhmer_tx_VenomZone/tx_VenomZone_aln || { echo "Error: Failed to navigate to hhmer_tx_VenomZone directory."; exit 1; }
bash run_hmmbuild_and_hmmsearch.sh || { echo "Error: hhmer_tx_VenomZone Module failed."; exit 1; }
echo "hhmer_tx_VenomZone module complete."
cd ../../ || { echo "Error: Failed to navigate back to the root directory."; exit 1; }

# Step 2: hhmerTxMatch Module
echo "Running hhmerTxMatch Module..."
cd hhmer_tx_VenomZone || { echo "Error: Failed to navigate to hhmer_tx_VenomZone directory."; exit 1; }
python3 hhmerTxMatch.py || { echo "Error: hhmerTxMatch Module failed."; exit 1; }
echo "hhmerTxMatch module complete."
cd .. || { echo "Error: Failed to navigate back to the root directory."; exit 1; }

# Step 3: TransDeTox Module
echo "Running TransDeTox Module..."
cd TransDeTox || { echo "Error: Failed to navigate to TransDeTox directory."; exit 1; }
python3 TransDeTox.py || { echo "Error: TransDeTox Module failed."; exit 1; }
echo "TransDeTox module complete."
cd .. || { echo "Error: Failed to navigate back to the root directory."; exit 1; }

# Step 4: ToxinKeyMatch Module
echo "Running ToxinKeyMatch Module..."
cd ToxinKeyMatch || { echo "Error: Failed to navigate to ToxinKeyMatch directory."; exit 1; }
python3 ToxinKeyMatch.py || { echo "Error: ToxinKeyMatch Module failed."; exit 1; }
echo "ToxinKeyMatch module complete."
cd .. || { echo "Error: Failed to navigate back to the root directory."; exit 1; }

# Step 5: FiltDeTox Module
echo "Running FiltDeTox Module..."
cd FiltDeTox || { echo "Error: Failed to navigate to FiltDeTox directory."; exit 1; }
Rscript ToxRecov.R || { echo "Error: FiltDeTox Module failed."; exit 1; }
echo "FiltDeTox module complete."
cd .. || { echo "Error: Failed to navigate back to the root directory."; exit 1; }

echo "Pipeline complete. All steps are finished."

