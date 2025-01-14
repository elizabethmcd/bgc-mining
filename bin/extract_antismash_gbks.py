#!/usr/bin/env python

import argparse
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class AntismashGBKParser:
    """Parse GBK file from antiSMASH to extract relevant BGC and peptide information."""

    def __init__(self, gbk_file):
        self.gbk_file = gbk_file
        self.records = list(SeqIO.parse(self.gbk_file, 'genbank'))
        
    def parse_records(self):
        """Parse all records in the GBK file."""
        bgc_data_list = []
        peptide_data_list = []
        fasta_records = []
        for record in self.records:
            scaffold_name = record.id
            bgc_type = self.get_bgc_type(record)
            bgc_length = len(record.seq)
            bgc_completeness = self.get_bgc_completeness(record)
            leader_seqs, core_seqs = self.ex_lanthipeptide(record)
            bgc_data_list.append({
                'scaffold_name': scaffold_name,
                'bgc_file': os.path.basename(self.gbk_file),
                'bgc_type': bgc_type,
                'bgc_length': bgc_length,
                'bgc_completeness': bgc_completeness,
            })
            if leader_seqs and core_seqs:
                peptide_data_list.append({
                    'scaffold_name': scaffold_name,
                    'bgc_file': os.path.basename(self.gbk_file),
                    'leader_seq': ';'.join(leader_seqs),
                    'core_seq': ';'.join(core_seqs),
                })
                for i, seq in enumerate(core_seqs, 1):
                    record = SeqRecord(Seq(seq), 
                                       id=f"{scaffold_name}_{os.path.basename(self.gbk_file)}_{i}", 
                                       description="")
                    fasta_records.append(record)
        return bgc_data_list, peptide_data_list, fasta_records
    
    def get_bgc_type(self, record):
        """Return the type of BGC region."""
        bgc_type_list = []
        for feature in record.features:
            if feature.type == 'region':
                bgc_type_list += feature.qualifiers.get('product', [])
        return '+'.join(set(bgc_type_list))

    def get_bgc_completeness(self, record):
        """Return the completeness of the BGC."""
        for feature in record.features:
            if feature.type == 'region':
                contig_edge = feature.qualifiers.get("contig_edge", ["complete"])[0]
                if contig_edge != "complete":
                    return "incomplete"
        return "complete"

    def ex_lanthipeptide(self, record):
        """Extract leader and core sequences of lanthipeptides (RiPPs)."""
        leader_seq, core_seq = [], []
        for feature in record.features:
            if feature.type == 'CDS_motif':
                qualifiers = feature.qualifiers
                if qualifiers.get("core_sequence") and qualifiers.get("leader_sequence"):
                    leader_seq += qualifiers.get("leader_sequence")
                    core_seq += qualifiers.get("core_sequence")
        return leader_seq, core_seq

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Parse antiSMASH GBK files to summarize BGCs and extract lanthipeptides.")
    parser.add_argument("input_gbk", nargs='+', help="Paths to the input GBK files")
    parser.add_argument("output_bgc_tsv", help="Path to the output TSV file for BGC information")
    parser.add_argument("output_peptide_tsv", help="Path to the output TSV file for lanthipeptide information")
    parser.add_argument("output_peptide_fasta", help="Path to the output FASTA file for lanthipeptide sequences")
    return parser.parse_args()

def main():
    args = parse_args()

    all_bgc_data = []
    all_peptide_data = []
    all_fasta_records = []
    for gbk_file in args.input_gbk:
        parser = AntismashGBKParser(gbk_file)
        bgc_data, peptide_data, fasta_records = parser.parse_records()
        all_bgc_data.extend(bgc_data)
        all_peptide_data.extend(peptide_data)
        all_fasta_records.extend(fasta_records)

    # Write BGC data to TSV
    bgc_df = pd.DataFrame(all_bgc_data)
    bgc_df.to_csv(args.output_bgc_tsv, sep='\t', index=False)
    print(f"BGC data written to {args.output_bgc_tsv}")

    # Write peptide data to TSV if found
    if all_peptide_data:
        peptide_df = pd.DataFrame(all_peptide_data)
        peptide_df.to_csv(args.output_peptide_tsv, sep='\t', index=False)
        print(f"Lanthipeptides found and written to {args.output_peptide_tsv}.")

        # Write FASTA records
        SeqIO.write(all_fasta_records, args.output_peptide_fasta, "fasta")
        print(f"Core peptide sequences written to {args.output_peptide_fasta}.")
    else:
        print(f"No lanthipeptides found in the input GBK files. No peptide files generated.")

if __name__ == "__main__":
    main()