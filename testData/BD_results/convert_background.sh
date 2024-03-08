#!/bin/bash

# Create a white background image

# Function to get dimensions of a file
get_dimensions() {
    identify -format "%wx%h" "$1"
}

# Loop through all PNG files in the current directory
for file in *.png; do
    # Check if the file exists and is a regular file
    if [ -f "$file" ]; then
	dimensions=$(get_dimensions "$file")
	convert -size $dimensions xc:white white_bg.png
        # Composite the PNG file onto the white background
        composite -gravity center "$file" white_bg.png "white_$file"
    fi
done

# Remove the temporary white background image
rm white_bg.png
