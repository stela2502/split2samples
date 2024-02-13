#!/bin/bash

# Check if RustodyFiles environment variable is set in ~/.bashrc
if ! grep -q "export RustodyFiles=" ~/.bashrc; then
    # Add RustodyFiles environment variable to ~/.bashrc
    echo 'export RustodyFiles="/opt/Rustody_files"' >> ~/.bashrc
fi

# Create directory and copy files if not already exists
mkdir -p /opt/Rustody_files
cp ./resources/CellRanger/*.gz /opt/Rustody_files
export RustodyFiles="/opt/Rustody_files"

cp target/release/split2samples /usr/bin
cp target/release/quantify_rhapsody /usr/bin
cp target/release/quantify_rhapsody_multi /usr/bin
cp target/release/bd_cell_id_2_seq /usr/bin
cp target/release/bd_get_single_cell /usr/bin
cp target/release/get_n_cell_reads /usr/bin
cp target/release/int_2_seq /usr/bin
cp target/release/create_index /usr/bin
cp target/release/create_index_te /usr/bin
cp target/release/te_analysis /usr/bin