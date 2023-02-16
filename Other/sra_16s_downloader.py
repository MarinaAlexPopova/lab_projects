# !/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.04.2022
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru


import argparse
import os
from Bio import Entrez
import xml.etree.cElementTree as ET
import shutil


class NoToolException(Exception):
    pass


def is_tool(program_name: str) -> bool:
    """Check whether program_name is on PATH and executable."""
    return shutil.which(program_name) is not None


def get_ids(query: str, email: str, batch: int, start: int) -> list:
    """Finding IDs of sra that fit query"""
    Entrez.email = email
    handle = Entrez.esearch(db="sra", term=query, retmax=batch, retstart=start)
    record = Entrez.read(handle)
    ids = record["IdList"]
    return ids


def get_sra_names(ids: list) -> list:
    """Finding SRA names for IDs"""
    sra = []
    for id in ids:
        fetch_handle = Entrez.efetch(db="sra", id=id)
        data = fetch_handle.read()
        root = ET.fromstring(data)
        for rank in root.iter("RUN"):
            sra.append(rank.attrib["accession"])
    print(f"{len(sra)} fastq files will be downloaded")
    return sra


def downloading_commands(sra: list, output_folder: os.path) -> os.path:
    """Making a file with commands for xargs"""
    commands = []
    for sra_file in sra:
        command = f"fastq-dump --split-3 {sra_file}"
        commands.append(command)
    file_name = os.path.join(output_folder, "commands.txt")
    with open(file_name, "w") as file:
        file.write("\n".join(commands))
    return file_name


def download_sra(commands_list: path, threads: int):
    """Downloading data"""
    if not is_tool("fastq-dump"):
        exception = "Please, install 'mamba install -c bioconda sra-tools'"
        raise NoToolException(exception)
    command = f"less {commands_list} | xargs -I {{}} -P {threads} -n 1 sh -c {{}}"
    print("Starting download data. Please wait.")
    os.system(command)
    print("Done with this batch!")
    pass


def main(settings):
    """ """
    query = settings["query"]
    output = settings["output"]
    threads = settings["threads"]
    batch = settings["batch"]
    email = settings["email"]
    start = 0
    while True:
        ids = get_ids(query, email, batch, start)
        if len(ids) == 0:
            break
        sra = get_sra_names(ids)
        commands_list = downloading_commands(sra, output)
        download_sra(commands_list, threads)
        os.remove(commands_list)
        start += batch
    print("All 16S data was downloaded!")
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Downloading raw data")
    parser.add_argument(
        "-q",
        "--query",
        help="Query that you want to find",
        required=False,
        default='16s[All Fields] AND (cluster_public[prop] AND "biomol dna"[Properties])',
    )
    parser.add_argument(
        "-o", "--output", help="Ouput folder", required=False, default="./"
    )
    parser.add_argument(
        "-b",
        "--batch",
        help="Batches that the whole raw data will be divided and downloaded",
        required=False,
        default=100,
    )
    parser.add_argument("-t", "--threads", help="Threads", required=False, default=50)
    parser.add_argument(
        "-e",
        "--email",
        help="Threads that will be used",
        required=False,
        default="a@gmail.com",
    )
    args = vars(parser.parse_args())
    query = args["query"]
    output = args["output"]
    threads = args["threads"]
    batch = args["batch"]
    e_mail = args["email"]
    settings = {
        "query": query,
        "output": output,
        "batch": batch,
        "threads": threads,
        "email": e_mail,
    }
    main(settings)
