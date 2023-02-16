#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 28.11.2022
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru


import BCBio
import json
import click
from BCBio.GFF import GFFExaminer
from Bio import SeqIO
from tqdm import tqdm
from collections import defaultdict, Counter
from multiprocessing import Pool
from datetime import datetime


# function for creating feature lists for plus and minus strands for single chromosome
def creating_feat_list(record):
    print("creating feature lists")
    plus = []
    minus = []
    for i in record.features[1:]:
        if i.location.strand != 1:
            minus.append(i)
        if i.location.strand != -1:
            plus.append(i)
    plus = [x for x in plus if "match" not in x.type]
    minus = [x for x in minus if "match" not in x.type]
    plus_list = [
        [
            x.location.start,
            x.location.end,
            x.qualifiers["gene_biotype"][0]
            if "gene_biotype" in x.qualifiers
            else x.type,
        ]
        for x in plus
    ]
    minus_list = [
        [
            x.location.start,
            x.location.end,
            x.qualifiers["gene_biotype"][0]
            if "gene_biotype" in x.qualifiers
            else x.type,
        ]
        for x in minus
    ]
    return plus_list, minus_list


# forming list of features for the whole chromosome
def feature_per_nucl_whole_chrom(seq_length, gff_features):
    n_type = "not annotated"
    type_list = ["not annotated"]
    features = []
    feature_id = 0
    feature_list = []
    feature_queue = []
    length_of_features = len(gff_features)
    for i in tqdm(range(seq_length)):
        if feature_id >= length_of_features:
            if len(feature_queue) == 0:
                feature_list.append(n_type)
            elif i <= feature_queue[0][1]:
                feature_list.append(n_type)
            elif i > feature_queue[0][1]:
                t = feature_queue[0][2]
                type_list.remove(t)
                if type_list:
                    type_l = []
                    [type_l.append(x) for x in type_list if x not in type_l]
                    n_type = "/".join(type_l)
                feature_queue = feature_queue[1:]
                if not feature_queue:
                    n_type = "not annotated"
                    type_list = ["not annotated"]
                    feature_list.append(n_type)
                else:
                    feature_list.append(n_type)
        elif i < gff_features[feature_id][0]:
            if len(feature_queue) == 0:
                feature_list.append(n_type)
            elif i <= feature_queue[0][1]:
                feature_list.append(n_type)
            elif i > feature_queue[0][1]:
                t = feature_queue[0][2]
                type_list.remove(t)
                if type_list:
                    type_l = []
                    [type_l.append(x) for x in type_list if x not in type_l]
                    n_type = "/".join(type_l)
                feature_queue = feature_queue[1:]
                if not feature_queue:
                    n_type = "not annotated"
                    type_list = ["not annotated"]
                    feature_list.append(n_type)
                else:
                    feature_list.append(n_type)
        elif i == gff_features[feature_id][0]:
            k = 1
            while True:
                if feature_id + k >= length_of_features:
                    feature_id = feature_id + k - 1
                    break
                elif i == gff_features[feature_id + k][0]:
                    feature_queue.append(gff_features[feature_id + k])
                    k += 1
                else:
                    feature_id = feature_id + k - 1
                    break

            feature_queue.append(gff_features[feature_id])
            type_list = [x[2] for x in feature_queue]
            type_l = []
            if feature_id == length_of_features:
                print(feature_id)
            [type_l.append(x) for x in type_list if x not in type_l]
            n_type = "/".join(type_l)
            feature_queue = sorted(
                feature_queue, key=lambda feature_entry: feature_entry[1]
            )
            feature_id += 1
            feature_list.append(n_type)

    return feature_list


# creating  kmers dictionary (with k = range[5, 10]) with their features for one chromosome
def kmer_dictionary(record):
    plus_list, minus_list = creating_feat_list(record)
    chr_length = len(record.seq)
    print("Creating features for plus strand:")
    feat_list_plus = feature_per_nucl_whole_chrom(chr_length, plus_list)
    print("Creating features for minus strand:")
    feat_list_minus = feature_per_nucl_whole_chrom(chr_length, minus_list)

    kmers_dict = defaultdict(Counter)

    # plus strand
    for i in tqdm(range(5, 11)):
        print(f"Counting for kmer = {i}. Plus strand")
        sequence = list(record.seq[::-1])
        features_list = feat_list_plus[::-1]

        kmer = []
        features = []
        for i in range(i - 1):
            kmer.append(sequence.pop())
            features.append(features_list.pop())

        while features_list:
            kmer.append(sequence.pop())
            features.append(features_list.pop())
            kmer_string = "".join(kmer)
            result_features = []
            [result_features.append(x) for x in features if x not in result_features]
            kmers_dict[kmer_string.upper()] += Counter(["-".join(result_features)])
            kmer = kmer[1:]
            features = features[1:]

    # minus strand
    for i in tqdm(range(5, 11)):
        print(f"Counting for kmer = {i}. Minus strand")
        sequence = list(record.seq.complement())
        features_list = feat_list_minus.copy()

        kmer = []
        features = []
        for i in range(i - 1):
            kmer.append(sequence.pop())
            features.append(features_list.pop())

        while features_list:
            kmer.append(sequence.pop())
            features.append(features_list.pop())
            kmer_string = "".join(kmer)
            result_features = []
            [result_features.append(x) for x in features if x not in result_features]
            kmers_dict[kmer_string.upper()] += Counter(["-".join(result_features)])
            kmer = kmer[1:]
            features = features[1:]

    return kmers_dict


def pool_handler(n, records):
    p = Pool(n)
    results = p.map(kmer_dictionary, records)
    p.close()
    return results


def save_results(res, out):
    combined_dict = defaultdict(Counter)
    for d in res:
        for key, value in d.items():
            combined_dict[key] += Counter(value)
    json.dump(combined_dict, open(out, "w"))


@click.command()
@click.option("-f", "--fasta", type=click.Path(), required=True, help="fasta file")
@click.option("-g", "--gff", type=click.Path(), required=True, help="gff file")
@click.option("-t", "--thread", type=int, default=5, help="threads")
@click.option("-o", "--out", type=click.Path(), default="out.json", help="gff file")
def token(fasta, gff, thread, out):
    start = datetime.now()
    examiner = GFFExaminer()
    with open(fasta) as f:
        seq_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
    print("fasta is opened")
    records = []
    with open(gff) as f:
        for rec in BCBio.GFF.parse(f, base_dict=seq_dict):
            records.append(rec)
    print("gff is opened")
    results = pool_handler(thread, records)
    save_results(results, out)
    finish = datetime.now()
    time = finish - start
    print(f"spent time {time}")


if __name__ == "__main__":
    token()
