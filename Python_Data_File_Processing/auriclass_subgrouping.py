!/usr/bin/env python3
import csv
import subprocess
import os
import argparse

# Parse command-line arguments.
parser = argparse.ArgumentParser(
    description="Assign isolates to clades using auriclass and update clade files."
)
parser.add_argument(
    "--reads",
    required=True,
    help="Path to the input reads.tab file (tab-delimited with Isolate_ID, reads1, reads2).",
)
parser.add_argument(
    "--contigs",
    required=True,
    help="Path to the input contigs.tab file (tab-delimited with Isolate_ID, contigs).",
)
parser.add_argument(
    "--clade_dir",
    required=True,
    help="Path to the directory where the current clade files exist. New entries will be appended here.",
)
args = parser.parse_args()

# Assign command-line arguments to variables.
READS_TAB = args.reads
CONTIGS_TAB = args.contigs
CLADE_DIR = args.clade_dir

# Ensure the clade directory exists.
if not os.path.exists(CLADE_DIR):
    os.makedirs(CLADE_DIR)

# Define allowed clades.
allowed_clades = ["Clade I", "Clade II", "Clade III", "Clade IV", "Clade V", "Clade VI"]

# Define the summary file name (now called Clade_designations.tab) within the clade directory.
SUMMARY_FILE = os.path.join(CLADE_DIR, "Clade_designations.tab")

# Global dictionary to track which isolates have been added to each file.
# Keys are filenames and values are sets of Isolate_IDs.
file_written_ids = {}

def append_line_if_not_exists(filename, line, isolate_id, header_line=None):
    """
    Append 'line' to file 'filename' only if 'isolate_id' is not already present.
    If header_line is provided, it will be used to skip over a header line when reading existing entries.
    Returns True if the line was appended, False if it was skipped.
    """
    # If we haven't yet loaded the file's content, load it.
    if filename not in file_written_ids:
        existing_ids = set()
        if os.path.exists(filename):
            with open(filename, "r") as f:
                for l in f:
                    if not l.strip():
                        continue  # Skip blank lines.
                    # Skip header if provided.
                    if header_line is not None and l.strip() == header_line.strip():
                        continue
                    # Assume the first column is the Isolate_ID.
                    existing_isolate = l.split("\t")[0].strip()
                    existing_ids.add(existing_isolate)
        file_written_ids[filename] = existing_ids

    if isolate_id in file_written_ids[filename]:
        print(f"{isolate_id} already exists in {filename}. Skipping.")
        return False
    else:
        with open(filename, "a") as out_f:
            out_f.write(line)
        file_written_ids[filename].add(isolate_id)
        return True

# Create the summary file with header if it doesn't exist.
if not os.path.exists(SUMMARY_FILE):
    with open(SUMMARY_FILE, "w") as summary_fh:
        summary_fh.write("Isolate_ID\tClade\n")
    file_written_ids[SUMMARY_FILE] = set()

# Load the reads file into a dictionary:
# Format: { Isolate_ID: (reads1, reads2) }
reads_dict = {}
with open(READS_TAB, "r") as fh:
    reader = csv.reader(fh, delimiter="\t")
    for row in reader:
        # Skip empty or commented lines.
        if not row or row[0].startswith("#"):
            continue
        try:
            isolate_id, reads1, reads2 = row
        except ValueError:
            print(f"Skipping malformed row in {READS_TAB}: {row}")
            continue
        reads_dict[isolate_id] = (reads1, reads2)

# Load the contigs file into a dictionary:
# Format: { Isolate_ID: contigs }
contigs_dict = {}
with open(CONTIGS_TAB, "r") as fh:
    reader = csv.reader(fh, delimiter="\t")
    for row in reader:
        if not row or row[0].startswith("#"):
            continue
        try:
            isolate_id, contigs = row
        except ValueError:
            print(f"Skipping malformed row in {CONTIGS_TAB}: {row}")
            continue
        contigs_dict[isolate_id] = contigs

# Function to run auriclass on a given reads file and return the report path.
def run_auriclass(reads1_path, isolate_id):
    # Construct the command; adjust parameters as needed.
    cmd = ["auriclass", reads1_path]
    print(f"Running auriclass for {isolate_id}: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running auriclass for {isolate_id}:\n{result.stderr}")
        return None
    # Assuming auriclass writes its output to a file named 'report.tsv' in the current directory.
    report_path = "report.tsv"
    if not os.path.exists(report_path):
        print(f"Report file {report_path} not found for {isolate_id}")
        return None
    return report_path

# Function to parse the auriclass report.
# Assumes the report is a tab-delimited file with headers,
# and that the first data row contains the clade value in the second column.
def parse_report(report_path):
    with open(report_path, "r") as rep_fh:
        reader = csv.reader(rep_fh, delimiter="\t")
        # Read header row (assuming header exists).
        header = next(reader, None)
        if header is None:
            print(f"No header found in {report_path}")
            return None
        # Read the first data row; auriclass outputs one row per run.
        row = next(reader, None)
        if row is None or len(row) < 2:
            print(f"Report file {report_path} does not contain expected columns")
            return None
        clade = row[1].strip()  # Extract clade from the second column.
        print(f"Extracted clade: {clade}")
        return clade

# Process each sample.
for isolate_id, (reads1, reads2) in reads_dict.items():
    print(f"\nProcessing sample {isolate_id} ...")
    
    # Run auriclass using the first reads file (reads1).
    report_path = run_auriclass(reads1, isolate_id)
    if report_path is None:
        print(f"Skipping {isolate_id} due to auriclass error.")
        continue

    # Parse the report to extract the clade from the second column.
    clade = parse_report(report_path)
    if clade is None:
        print(f"Skipping {isolate_id} because no clade was extracted.")
        continue

    # Check if the clade is one of the allowed clades.
    if clade not in allowed_clades:
        print(f"Warning: Clade '{clade}' for {isolate_id} is not in the allowed list. Skipping sample.")
        continue

    # Prepare output file names for the clade.
    # Replace spaces with underscores (e.g. "Clade III" becomes "Clade_III")
    safe_clade = clade.replace(" ", "_")
    clade_reads_file = os.path.join(CLADE_DIR, f"{safe_clade}.reads.tab")
    clade_contigs_file = os.path.join(CLADE_DIR, f"{safe_clade}.contigs.tab")

    # Prepare lines to append.
    reads_line = f"{isolate_id}\t{reads1}\t{reads2}\n"
    contigs_path = contigs_dict.get(isolate_id)
    contigs_line = f"{isolate_id}\t{contigs_path}\n" if contigs_path else None
    summary_line = f"{isolate_id}\t{clade}\n"

    # Append the sample's reads entry to the clade-specific reads file if not already present.
    append_line_if_not_exists(clade_reads_file, reads_line, isolate_id)

    # Append the sample's contigs entry to the clade-specific contigs file if not already present.
    if contigs_line:
        append_line_if_not_exists(clade_contigs_file, contigs_line, isolate_id)
    else:
        print(f"Warning: No contigs entry found for {isolate_id}")

    # Append the sample's result (Isolate_ID and Clade) to the summary file.
    append_line_if_not_exists(SUMMARY_FILE, summary_line, isolate_id, header_line="Isolate_ID\tClade")

    # Optionally, remove the temporary report file if not needed:
    # os.remove(report_path)

    print(f"Assigned sample {isolate_id} to {clade} and updated files:")
    print(f"  - {clade_reads_file}")
    print(f"  - {clade_contigs_file}")
    print(f"  - {SUMMARY_FILE}")

print("\nPipeline completed.")