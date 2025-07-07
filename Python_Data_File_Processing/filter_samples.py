###this script takes a core alignment and corresposding tree file and will remove samples based on a list of IDs that match both files. new tree and aln is outputted

#!/usr/bin/env python3
"""
filter_samples.py

Filter out specified samples from a FASTA alignment and Newick tree.
Usage:
    python filter_samples.py \  
        --alignment core.aln \  
        --tree treefile.nwk \  
        --remove remove_list.txt \  
        --out-aln filtered.aln \  
        --out-tree filtered.nwk
"""
import argparse
from Bio import SeqIO, Phylo

def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter samples from alignment and tree"
    )
    parser.add_argument(
        '--alignment', required=True,
        help='Input FASTA alignment file'
    )
    parser.add_argument(
        '--tree', required=True,
        help='Input Newick tree file'
    )
    parser.add_argument(
        '--remove', required=True,
        help='Text file with sample IDs to remove, one per line'
    )
    parser.add_argument(
        '--out-aln', required=True,
        help='Output filtered alignment file'
    )
    parser.add_argument(
        '--out-tree', required=True,
        help='Output filtered tree file'
    )
    return parser.parse_args()

def load_remove_list(path):
    """Read sample IDs to remove."""
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}

def filter_alignment(in_fasta, out_fasta, remove_ids):
    """Write sequences whose IDs are not in remove_ids."""
    records = (
        rec for rec in SeqIO.parse(in_fasta, 'fasta')
        if rec.id not in remove_ids
    )
    SeqIO.write(records, out_fasta, 'fasta')

def filter_tree(in_tree, out_tree, remove_ids):
    """Prune tips with names in remove_ids and write new tree."""
    tree = Phylo.read(in_tree, 'newick')
    for tid in remove_ids:
        try:
            tree.prune(target=tid)
        except ValueError:
            # tip not found; skip
            pass
    Phylo.write(tree, out_tree, 'newick')

def main():
    args = parse_args()
    remove_ids = load_remove_list(args.remove)

    filter_alignment(args.alignment, args.out_aln, remove_ids)
    filter_tree(args.tree, args.out_tree, remove_ids)

if __name__ == '__main__':  # noqa: C401
    main()
