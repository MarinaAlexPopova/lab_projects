# !/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.04.2022
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru


import argparse
from importlib.resources import path
import os
from Bio import Entrez
import xml.etree.cElementTree as ET
import shutil


class NoToolException(Exception):
    pass


def is_tool(program_name: str) -> bool:
    """Check whether program_name is on PATH and executable."""
    return shutil.which(program_name) is not None


def get_ids(query: str, email: str, max_ids: int) -> list:
    """Finding IDs of sra that fit query"""
    Entrez.email = email
    handle = Entrez.esearch(db="sra", term=query, retmax=max_ids)
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


def downloading_commands(sra: list, output_folder: path) -> path:
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
    print("Done!")
    pass


def main(settings):
    """ """
    query = settings["query"]
    output = settings["output"]
    threads = settings["threads"]
    max_ids = settings["max_ids"]
    email = settings["email"]
    ids = get_ids(query, email, max_ids)
    sra = get_sra_names(ids)
    commands_list = downloading_commands(sra, output)
    download_sra(commands_list, threads)
    os.remove(commands_list)
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Downloading raw data")
    parser.add_argument(
        "-q", "--query", help="Query that you want to find", required=True
    )
    parser.add_argument(
        "-o", "--output", help="Ouput folder", required=False, default="./"
    )
    parser.add_argument(
        "-m",
        "--max_ids",
        help="Maximum of raw data that will be downloaded",
        required=False,
        default=2,
    )
    parser.add_argument("-t", "--threads", help="Threads", required=False, default=1)
    parser.add_argument(
        "-e", "--email", help="Threads", required=False, default="a@gmail.com"
    )
    args = vars(parser.parse_args())
    query = args["query"]
    output = args["output"]
    threads = args["threads"]
    max_ids = args["max_ids"]
    e_mail = args["email"]
    settings = {
        "query": query,
        "output": output,
        "max_ids": max_ids,
        "threads": threads,
        "email": e_mail,
    }
    main(settings)
