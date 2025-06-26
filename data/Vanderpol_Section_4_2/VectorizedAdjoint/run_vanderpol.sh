#!/bin/bash

# Create data directory if it doesn't exist
mkdir -p data

# Define the CSV file path for storing results
output_file="data/data_vanderpol_CK54.csv"

# Remove the output file if it exists to clean it
if [ -f "$output_file" ]; then
    rm "$output_file"
fi


# Write header to CSV file
#echo "tol,y,z,dydmu,dzdmu" > "$output_file"

# Path to your executable relative to script
executable="./build/vanderpol"

# Array of tolerance values to iterate over
tolerances=(1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-12 1e-13)

# Loop over each tolerance value
for tol in "${tolerances[@]}"; do
    # Run your C++ program with current tolerance value
    "$executable" "$tol"

    # Extract results from output.csv
    result_line=$(tail -n 1 output.csv)

    # Append tolerance and results to the main CSV file
    echo "$tol,$result_line" >> "$output_file"
done

echo "All simulations completed. Results stored in $output_file"

# Remove the temporary output file
if [ -f "output.csv" ]; then
    rm output.csv
    echo "Temporary output file output.csv removed."
else
    echo "No temporary output file found."
fi  
