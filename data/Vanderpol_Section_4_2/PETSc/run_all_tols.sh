#!/bin/bash

# Create data directory if it doesn't exist
mkdir -p data

# List of tolerances
tols=(1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-13)

# Number of repetitions per tolerance
TREP=3

# Path to executable
EXEC=./build/vanderpol

# Output file
OUTPUT_FILE="data/data.dat"
echo -n > "$OUTPUT_FILE"  # Clear or create the file

# Max time step file:
STEP_FILE="data/steps.dat"
echo -n > "$STEP_FILE"  # Clear or create the file

OUTPUT_DATA_FILE="data/data_vanderpol_dopri5_petsc.dat"
echo -n > "$OUTPUT_DATA_FILE"  # Clear or create the file

OUTPUT_DATA_IMEX_FILE="data/data_vanderpol_imex_petsc.dat"
echo -n > "$OUTPUT_DATA_IMEX_FILE"  # Clear or create the file


# Associative array to hold time lists
declare -A times

# Loop over tolerances and collect timing data
for tol in "${tols[@]}"; do
    echo "Running tolerance $tol..."
    
    output=$($EXEC -tol "$tol")

    read -r dtmax < <(
    echo "$output" | grep -i "^Max time step used:" | sed -E 's/^Max time step used:[[:space:]]*//')
    echo "Got max dt: $dtmax"
    echo "$dtmax " >> "$STEP_FILE"



    # Extract the five values into an array called “data”
    read -r tol x z dxdmu dzdmu < <(
    echo "$output" | grep -i "^Data is:" | sed -E 's/^Data is:[[:space:]]*//')

    # Now $d1 … $d5 hold the five doubles
    echo "Got ERK values: $tol $x $z $dxdmu $dzdmu"

    echo "$tol $x $z $dxdmu $dzdmu" >> "$OUTPUT_DATA_FILE"

    output=$($EXEC -tol "$tol" -imexform true)
    # Extract the five values into an array called “data”
    read -r tol x z dxdmu dzdmu < <(
    echo "$output" | grep -i "^Data is:" | sed -E 's/^Data is:[[:space:]]*//')
    
    # Now $d1 … $d5 hold the five doubles
    echo "Got IMEX values: $tol $x $z $dxdmu $dzdmu"
    echo "$tol $x $z $dxdmu $dzdmu" >> "$OUTPUT_DATA_IMEX_FILE"


    # Time
    for ((i=1; i<=TREP; i++)); do
        echo "  Run $i/$TREP"
        output=$($EXEC -tol "$tol")

        time_taken=$(echo "$output" | grep -i "Time taken for the event" | sed -E 's/.*: ([0-9.]+) seconds/\1/')
        times["$tol"]+="$time_taken "
    done
done

# Compute stats and write to data.dat
for tol in "${tols[@]}"; do
    time_list="${times[$tol]}"
    
    # Use awk to compute mean and stddev
    stats=$(echo "$time_list" | awk '
    {
        n=NF;
        sum=0;
        for (i=1; i<=n; i++) sum += $i;
        mean = sum / n;
        sumsq = 0;
        for (i=1; i<=n; i++) sumsq += ($i - mean)^2;
        stddev = sqrt(sumsq / n);
        printf "%.10f %.10f", mean, stddev;
    }')

    echo "$tol $stats" >> "$OUTPUT_FILE"
done

echo "Results saved in data folder"