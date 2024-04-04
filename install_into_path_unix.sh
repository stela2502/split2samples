#!/usr/bin/bash

prefix=$1

if ! -d $1; then
    echo "sorry the path ${prefix} does not exists"
fi

path="${prefix}/.Rustody_files"
# Check if RustodyFiles environment variable is set in ~/.bashrc
if ! grep -q "export RustodyFiles=" ~/.bashrc; then
    # Add RustodyFiles environment variable to ~/.bashrc
    echo 'export RustodyFiles="${path}"' >> ~/.bashrc
fi

# Create directory and copy files if not already exists

if ! -d ${path}; then
    mkdir -p ${path}
fi

cp ./resources/CellRanger/*.gz ${path}
export RustodyFiles="${path}"

if ! -d ${prefix}/bin; then
    mkdir ${prefix}/bin
fi

cp target/release/split2samples $prefix/bin
cp target/release/quantify_rhapsody $prefix/bin
cp target/release/quantify_rhapsody_multi $prefix/bin
cp target/release/quantify_gene_mapper $prefix/bin
cp target/release/bd_cell_id_2_seq $prefix/bin
cp target/release/bd_get_single_cell $prefix/bin
cp target/release/get_n_cell_reads $prefix/bin
cp target/release/int_2_seq $prefix/bin
cp target/release/check_read_gene_mapper $prefix/bin
cp target/release/check_read_fast_mapper $prefix/bin
cp target/release/create_index $prefix/bin
cp target/release/create_index_te $prefix/bin
cp target/release/te_analysis $prefix/bin
cp target/release/create_gene_mapper_index $prefix/bin

echo "finished"