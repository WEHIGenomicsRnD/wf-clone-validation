#!/usr/bin/env python

"""Report component for displaying information from wf-clone-validation."""
import argparse
import json
import os

from aplanat import json_item
import numpy as np
import pandas as pd
from plannotate.annotate import annotate
from plannotate.bokeh_plot import get_bokeh
import pysam


def run_plannotate(fasta, linear=False):
    """Run annotate and create Bokeh plot."""
    with pysam.FastxFile(fasta) as fh:
        seq = next(fh).sequence
    df = annotate(
        seq, is_detailed=True, linear=linear, yaml_file="plannotate.yaml")
    plot = get_bokeh(df, linear=linear)
    plot.xgrid.grid_line_color = None
    plot.ygrid.grid_line_color = None
    old_df = df.copy(deep=True)
    clean_df = clean_results(df)
    old_df.reset_index(drop=True, inplace=True)
    return plot, clean_df, old_df


def clean_results(df):
    """Clean-up annotation dataframe for display."""
    rename = {
        'Feature': 'Feature',
        'db': 'Database',
        'pident': 'Identity', 'abs percmatch': 'Match Length',
        'Description': 'Description',
        'qstart': 'Start Location', 'qend': 'End Location', 'length': 'Length',
        'sframe': 'Strand', 'db': 'Database',
        'qlen': 'qlen'}
    numeric_columns = ['Identity', 'Match Length']
    display_columns = [
        'Feature',
        'Database',
        'Identity', 'Match Length',
        'Description',
        'Start Location', 'End Location', 'Length',
        'Strand', 'qlen']
    df = df.rename(columns=rename)[display_columns]
    df['Plasmid length'] = df.iloc[0]['qlen']
    df = df.drop(columns='qlen')
    df[numeric_columns] = np.round(df[numeric_columns], 1).astype(str) + "%"
    df.loc[df['Database'] == "infernal", 'Identity'] = "-"
    df.loc[df['Database'] == "infernal", 'Match Length'] = "-"
    df = df.set_index("Feature", drop=True).reset_index()
    return df


def bed_file(item, df):
    """Bed format for annotations."""
    display_columns = [
        'Start Location', 'End Location',
        'Feature',
        'Strand']
    df = df[display_columns]
    df["Strand"] = df["Strand"].apply(pd.to_numeric)
    df.loc[df['Strand'] == 0, 'Strand'] = "-"
    df.loc[df['Strand'] == 1, 'Strand'] = "+"
    df.insert(0, 'Name', value=str(item))
    df.to_csv(
        str(item)+'.annotations.bed', sep="\t", header=False, index=False)
    return df


def per_assembly(sample_file, item):
    """Run plannotate for a sample.

    :param database:
    :param sample_file:
    :param item: the sample
    """
    plot, annotations, clean_df = run_plannotate(sample_file)
    bed_file(item, annotations)
    with pysam.FastxFile(sample_file) as fh:
        seq_len = len(next(fh).sequence)
    tup = {
        'sample_name': item,
        'plot': plot,
        'annotations': annotations,
        'seq_len': seq_len}
    return tup, clean_df


def output_feature_table(data):
    """Build feature table text file or if no data output empty file."""
    if data:
        df = data[0]['annotations']
        sample_column = data[0]['sample_name']
        df.insert(0, 'Sample_name', sample_column)
        df.to_csv('feature_table.txt', mode='a', header=True, index=False)
        for sample in data[1:]:
            df = sample['annotations']
            sample_column = sample['sample_name']
            df.insert(0, 'Sample_name', sample_column)
            df.to_csv('feature_table.txt', mode='a', header=False, index=False)
    else:
        # If no samples passed create empty feature_table file
        feature_file = open("feature_table.txt", "w")
        feature_file.write(
"""Sample_name,Feature,Uniprot ID,Database,Identity,Match Length,Description,Start Location,End Location,Length,Strand,Plasmid length""") # noqa
        feature_file.close()


def make_yaml(database):
    """Create a yaml file for plannotate."""
    plannotate_yaml = """
Rfam:
  details:
    compressed: false
    default_type: ncRNA
    location: None
  location: {0}
  method: infernal
  priority: 3
  version: release 14.5
fpbase:
  details:
    compressed: false
    default_type: CDS
    location: Default
  location: {0}
  method: diamond
  parameters:
  - -k 0
  - --min-orf 1
  - --matrix BLOSUM90
  - --gapopen 10
  - --gapextend 1
  - --algo ctg
  - --id 75
  priority: 1
  version: downloaded 2020-09-02
snapgene:
  details:
    compressed: false
    default_type: None
    location: Default
  location: {0}
  method: blastn
  parameters:
  - -perc_identity 95
  - -max_target_seqs 20000
  - -culling_limit 25
  - -word_size 12
  priority: 1
  version: Downloaded 2021-07-23
swissprot:
  details:
    compressed: true
    default_type: CDS
    location: Default
  location: {0}
  method: diamond
  parameters:
  - -k 0
  - --min-orf 1
  - --matrix BLOSUM90
  - --gapopen 10
  - --gapextend 1
  - --algo ctg
  - --id 50
  priority: 2
  version: Release 2021_03
        """.format(database)

    with open("plannotate.yaml", "w") as text_file:
        text_file.write(plannotate_yaml)


def main():
    """Entry point to create a wf-clone-validation report."""
    parser = argparse.ArgumentParser(
        'Plannotate annotation',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)
    parser.add_argument(
        "--database", default='unknown',
        help="database to use, directory containing BLAST et. al. files.")
    parser.add_argument(
        "--sequences",
        help="sequences in directory to run plannotate on.",
        required=False)
    args = parser.parse_args()

    final_samples = []
    report_dic = {}
    plannotate_collection = {}
    make_yaml(args.database)
    json_file = open("plannotate.json", "a")
    if args.sequences:
        for filename in os.listdir(args.sequences):
            sample_file = os.path.join(args.sequences, filename)
            name = str(filename).split('.final')[0]
            tup_dic, clean_df = per_assembly(sample_file, name)
            final_samples.append(tup_dic)
            plasmid_len = tup_dic['annotations']['Plasmid length'][0]
            feature_dic = tup_dic['annotations'].drop(
                            ['Plasmid length'], axis=1)
            features = feature_dic.to_dict('records')
            output_json = json_item(tup_dic['plot'])
            plannotate_dic = {
                "reflen": float(plasmid_len),
                "features": features,
                "plot": output_json['doc']}
            plannotate_collection[name] = plannotate_dic
            report = {}
            report['sample_name'] = name
            report['plot'] = clean_df.to_json()
            report['annotations'] = tup_dic['annotations'].to_json()
            report['seq_len'] = tup_dic['seq_len']
            report_dic[name] = report

    # outputs for epi2me
    output_feature_table(final_samples)
    json_object = json.dumps(plannotate_collection, indent=4)
    json_file.write(json_object)
    json_file.close()

    # outputs for report
    json_file = open("plannotate_report.json", "a")
    json_object = json.dumps(report_dic)
    json_file.write(json_object)
    json_file.close()


if __name__ == "__main__":
    main()
