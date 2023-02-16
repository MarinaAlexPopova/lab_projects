#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 02.02.2023
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru


import pandas as pd
import os


def level_split(level, cell_types, assay_name):
    level = cell_types[f"predicted.ann_level_{level}"].unique()
    level_cells = {}
    for i in level:
        for j in list(
            cell_types.loc[
                (cell_types[f"predicted.ann_level_{level}"] == i)
                & (cell_types.assay_name == assay_name)
            ]["cell"]
        ):
            level_cells[j] = i
    return level_cells


def split_fastq_part1(assay_name, cell_types):

    level1_cells = level_split("1", cell_types, assay_name)
    level2_cells = level_split("2", cell_types, assay_name)
    level3_cells = level_split("3", cell_types, assay_name)

    barcodes = {}
    with open(f"/mnt/data/CIN/lung_cancer/german_raw_reads/{assay_name}.R2.fastq") as f:
        name = f.readline().strip()
        while True:
            if not name:
                break
            barcode = f.readline().strip()
            barcodes[name.split()[0]] = barcode
            f.readline()
            f.readline()
            name = f.readline().strip()

    # LEVEL1

    with open(f"/mnt/data/CIN/lung_cancer/german_raw_reads/{assay_name}.R1.fastq") as f:
        name = f.readline()
        while True:
            if not name:
                break
            read = f.readline()
            strand = f.readline()
            quality = f.readline()
            barcode = barcodes[name.split()[0]]
            if barcode in level1_cells:
                cell_type = level1_cells[barcode]
                with open(
                    f"/mnt/data/CIN/lung_cancer/german_raw_reads/kallisto/level1_{cell_type}_{assay_name}.fastq",
                    "a",
                ) as fw:
                    fw.write(name)
                    fw.write(read)
                    fw.write(strand)
                    fw.write(quality)
            name = f.readline()

    # LEVEL2

    with open(f"/mnt/data/CIN/lung_cancer/german_raw_reads/{assay_name}.R1.fastq") as f:
        name = f.readline()
        while True:
            if not name:
                break
            read = f.readline()
            strand = f.readline()
            quality = f.readline()
            barcode = barcodes[name.split()[0]]
            if barcode in level1_cells:
                cell_type = level2_cells[barcode]
                with open(
                    f"/mnt/data/CIN/lung_cancer/german_raw_reads/kallisto/level2_{cell_type}_{assay_name}.fastq",
                    "a",
                ) as fw:
                    fw.write(name)
                    fw.write(read)
                    fw.write(strand)
                    fw.write(quality)
            name = f.readline()

    # LEVEL3

    with open(f"/mnt/data/CIN/lung_cancer/german_raw_reads/{assay_name}.R1.fastq") as f:
        name = f.readline()
        while True:
            if not name:
                break
            read = f.readline()
            strand = f.readline()
            quality = f.readline()
            barcode = barcodes[name.split()[0]]
            if barcode in level1_cells:
                cell_type = level3_cells[barcode]
                with open(
                    f"/mnt/data/CIN/lung_cancer/german_raw_reads/kallisto/level3_{cell_type}_{assay_name}.fastq",
                    "a",
                ) as fw:
                    fw.write(name)
                    fw.write(read)
                    fw.write(strand)
                    fw.write(quality)
            name = f.readline()


def split_fastq_part2(assay_name, cell_types):
    level1_cells = level_split("1", cell_types, assay_name)
    level2_cells = level_split("2", cell_types, assay_name)
    level3_cells = level_split("3", cell_types, assay_name)

    barcodes = {}
    with open(f"/mnt/data/CIN/lung_cancer/german_raw_reads/{assay_name}_R1.fastq") as f:
        name = f.readline().strip()
        while True:
            if not name:
                break
            barcode = f.readline().strip()
            barcodes[name.split()[0]] = barcode[:16]
            f.readline()
            f.readline()
            name = f.readline().strip()

    # LEVEL1

    with open(f"/mnt/data/CIN/lung_cancer/german_raw_reads/{assay_name}_R2.fastq") as f:
        name = f.readline()
        while True:
            if not name:
                break
            read = f.readline()
            strand = f.readline()
            quality = f.readline()
            barcode = barcodes[name.split()[0]]
            if barcode in level1_cells:
                cell_type = level1_cells[barcode]
                with open(
                    f"/mnt/data/CIN/lung_cancer/german_raw_reads/kallisto/level1_{cell_type}_{assay_name}.fastq",
                    "a",
                ) as fw:
                    fw.write(name)
                    fw.write(read)
                    fw.write(strand)
                    fw.write(quality)
            name = f.readline()

    # LEVEL2

    with open(f"/mnt/data/CIN/lung_cancer/german_raw_reads/{assay_name}_R2.fastq") as f:
        name = f.readline()
        while True:
            if not name:
                break
            read = f.readline()
            strand = f.readline()
            quality = f.readline()
            barcode = barcodes[name.split()[0]]
            if barcode in level1_cells:
                cell_type = level2_cells[barcode]
                with open(
                    f"/mnt/data/CIN/lung_cancer/german_raw_reads/kallisto/level2_{cell_type}_{assay_name}.fastq",
                    "a",
                ) as fw:
                    fw.write(name)
                    fw.write(read)
                    fw.write(strand)
                    fw.write(quality)
            name = f.readline()

    # LEVEL3

    with open(f"/mnt/data/CIN/lung_cancer/german_raw_reads/{assay_name}_R2.fastq") as f:
        name = f.readline()
        while True:
            if not name:
                break
            read = f.readline()
            strand = f.readline()
            quality = f.readline()
            barcode = barcodes[name.split()[0]]
            if barcode in level1_cells:
                cell_type = level3_cells[barcode]
                with open(
                    f"/mnt/data/CIN/lung_cancer/german_raw_reads/kallisto/level3_{cell_type}_{assay_name}.fastq",
                    "a",
                ) as fw:
                    fw.write(name)
                    fw.write(read)
                    fw.write(strand)
                    fw.write(quality)
            name = f.readline()


if __name__ == "__main__":
    cell_types = pd.read_csv(
        "/mnt/data/satelome/users/mpopova/natella/azimuth_pred_with_sat.tsv", sep="\t"
    )
    cell_types = cell_types.loc[
        (cell_types["predicted.ann_level_1.score"] > 0.9)
        & (cell_types["predicted.ann_level_2.score"] > 0.8)
        & (cell_types["predicted.ann_level_3.score"] > 0.7)
    ]
    cell_types["assay_name"] = cell_types["cell"].apply(lambda x: x.split(":")[1])
    cell_types["cell"] = cell_types["cell"].apply(lambda x: x.split(":")[0])

    name_list_1 = ["1247", "BT1A", "BT1249", "BT1C", "BT2B", "BT1B", "BT2A"]
    for assay_name in name_list_1:
        command = f"wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6149/{assay_name}.R1.fastq.gz"
        os.system(command)
        command = f"wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6149/{assay_name}.R2.fastq.gz"
        os.system(command)
        command = f"gzip -d /mnt/data/CIN/lung_cancer/german_raw_reads/{assay_name}*"
        os.system(command)
        split_fastq_part1(assay_name, cell_types)
        command = f"rm /mnt/data/CIN/lung_cancer/german_raw_reads/{assay_name}*"
        os.system(command)
        print(f"{assay_name} splitted!")

    name_list_2 = [
        "BT1290",
        "BT1291",
        "BT1292",
        "BT1293",
        "BT1294",
        "BT1295",
        "BT1296",
        "BT1297",
        "BT1298",
        "BT1299",
        "BT1300",
        "BT1301",
    ]
    for assay_name in name_list_2:
        command = f"wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6149/{assay_name}_R1.fastq.gz"
        os.system(command)
        command = f"wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6149/{assay_name}_R2.fastq.gz"
        os.system(command)
        command = f"gzip -d /mnt/data/CIN/lung_cancer/german_raw_reads/{assay_name}*"
        os.system(command)
        split_fastq_part2(assay_name, cell_types)
        command = f"rm /mnt/data/CIN/lung_cancer/german_raw_reads/{assay_name}*"
        os.system(command)
        print(f"{assay_name} splitted!")
