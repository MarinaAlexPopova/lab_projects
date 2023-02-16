#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# @created: 13.02.2023
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru

import numpy as np
import os

if __name__ == "__main__":
    cell_types = [
        "level1_Endothelial",
        "level1_Epithelial",
        "level1_Immune",
        "level1_Stroma",
        "level2_Airway_epithelium",
        "level2_Alveolar_epithelium",
        "level2_Blood_vessels",
        "level2_Fibroblast_lineage",
        "level2_Lymphatic_EC",
        "level2_Lymphoid",
        "level2_Mesothelium",
        "level2_Myeloid",
        "level2_Smooth_muscle",
        "level3_AT1",
        "level3_AT2",
        "level3_B_cell_lineage",
        "level3_Basal",
        "level3_Dendritic_cells",
        "level3_EC_arterial",
        "level3_EC_capillary",
        "level3_EC_venous",
        "level3_Fibroblasts",
        "level3_Innate_lymphoid_cell_NK",
        "level3_Lymphatic_EC_mature",
        "level3_Macrophages",
        "level3_Mast_cells",
        "level3_Monocytes",
        "level3_Multiciliated_lineage",
        "level3_Myofibroblasts",
        "level3_None",
        "level3_Rare",
        "level3_Secretory",
        "level3_T_cell_lineage",
    ]

    directory = "/mnt/data/CIN/lung_cancer/german_raw_reads/kallisto"

    for i in cell_types:
        total = []
        with open(f"{directory}/{i}.fastq") as file:
            file.readline()
            seq = file.readline().strip()
            while seq:
                total.append(len(seq))
                file.readline()
                file.readline()
                file.readline()
                seq = file.readline().strip()
        mean = np.mean(total)
        std = np.std(total)
        command = f"kallisto quant -i {directory}/sat.index -l {mean} -s {std} -o {directory}/kallisto_count_results/{i} -b 1000 -t 80 --single {directory}/{i}.fastq"
        os.system(command)
