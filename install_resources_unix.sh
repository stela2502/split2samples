#!/bin/bash

path="/opt/Rustody_files"
# Check if RustodyFiles environment variable is set in ~/.bashrc
if ! grep -q "export RustodyFiles=" ~/.bashrc; then
    # Add RustodyFiles environment variable to ~/.bashrc
    echo 'export RustodyFiles="${path}"' >> ~/.bashrc
fi

# Create directory and copy files if not already exists
mkdir -p /opt/Rustody_files
cp ./resources/CellRanger/*.gz ${path}
export RustodyFiles="${path}"

cp target/release/split2samples /usr/local/bin
cp target/release/quantify_rhapsody /usr/local/bin
cp target/release/quantify_rhapsody_multi /usr/local/bin
cp target/release/quantify_gene_mapper /usr/local/bin
cp target/release/bd_cell_id_2_seq /usr/local/bin
cp target/release/bd_get_single_cell /usr/local/bin
cp target/release/get_n_cell_reads /usr/local/bin
cp target/release/int_2_seq /usr/local/bin
cp target/release/check_read_gene_mapper /usr/local/bin
cp target/release/check_read_fast_mapper /usr/local/bin
cp target/release/create_index /usr/local/bin
cp target/release/create_index_te /usr/local/bin
cp target/release/te_analysis /usr/local/bin
cp target/release/create_gene_mapper_index /usr/local/bin
