#!/bin/bash

# Default values
DB=""
QUERY=""
MAX_ATTEMPTS=5
DELAY=10  # Delay in seconds between retries

# Function to print usage
usage() {
    echo "Usage: $0 -db <database_path> -query <SQL_query>"
    exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -db) DB="$2"; shift ;;
        -query) QUERY="$2"; shift ;;
        *) usage ;;
    esac
    shift
done

# Check if both arguments are provided
if [[ -z "$DB" || -z "$QUERY" ]]; then
    usage
fi

sqlite3 $DB <<EOF
CREATE TABLE IF NOT EXISTS participantvariants (
    participantID TEXT,
    cohortID TEXT,
    ref TEXT,
    refpath TEXT,
    alignment TEXT,
    gvcf TEXT,
    vcf TEXT,
    PRIMARY KEY (participantID, cohortID, ref, refpath) 
);
EOF

sqlite3 $DB <<EOF
CREATE TABLE IF NOT EXISTS cohortvariants (
    cohortID TEXT NOT NULL,
    ref TEXT NOT NULL,
    chr TEXT NOT NULL,
    vcf TEXT NOT NULL,
    PRIMARY KEY (cohortID, ref, chr, vcf) 
);
EOF

# Retry mechanism
attempt=0
while [ $attempt -lt $MAX_ATTEMPTS ]; do
    sqlite3 "$DB" "$QUERY" && exit 0  # Success: Exit script
    attempt=$((attempt + 1))
    echo "Attempt $attempt failed. Retrying in $DELAY seconds..." >&2
    sleep $DELAY
done

# If all attempts fail, exit with an error
echo "ERROR: Failed to execute query after $MAX_ATTEMPTS attempts." >&2
exit 1