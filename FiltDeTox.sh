#!/bin/bash

# FilDeTox Pipeline Script
# This script runs the full FiltDeTox pipeline across all three modules.

# Step 1: TransDeTox Module
echo "Running TransDeTox Module..."
python3 TransDeTox/TransDeTox.py
echo "TransDeTox module complete."

# Step 2: ToxinKeyMatch Module
echo "Running ToxinKeyMatch Module..."
python3 ToxinKeyMatch/ToxinKeyMatch.py
echo "ToxinKeyMatch module complete."

# Step 3: FiltDeTox Module
echo "Running FiltDeTox Module..."
Rscript FiltDeTox/ToxRecov.R
echo "FiltDeTox module complete."

echo "Pipeline complete. All steps are finished."

