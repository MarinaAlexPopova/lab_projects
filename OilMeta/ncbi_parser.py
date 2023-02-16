#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 02.02.2023
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru

import re
import requests
import click
import pandas as pd
from Bio import Entrez
from collections import defaultdict


Entrez.email = "user@gmail.com"


def find_biosamples_id(term):
    handle = Entrez.esearch(db="bioproject", term=term, idtype="acc")
    record = Entrez.read(handle)
    handle.close()
    id_bioproject = record["IdList"][0]
    handle = Entrez.elink(
        dbfrom="bioproject", id=id_bioproject, linkname="bioproject_biosample"
    )
    record = Entrez.read(handle)
    handle.close()
    id_biosamples = []
    for i in record[0]["LinkSetDb"][0]["Link"]:
        id_biosamples.append(i["Id"])
    return id_biosamples


def find_sample_name(string):
    query = re.compile(r"(?<=Sample name: ).+(?=;)")
    result = query.search(string)
    if result:
        return result.group()
    else:
        return "None"


def find_srs(string):
    query = re.compile(r"(?<=SRA: ).+")
    result = query.search(string)
    if result:
        return result.group()
    else:
        return "None"


def find_location(string):
    query = re.compile(r'(?<="geographic location">).+?(?=</Attribute>)')
    result = query.search(string)
    if result:
        return result.group()
    else:
        return "None"


def find_lat_lon(string):
    query = re.compile(r'(?<="latitude and longitude">).+?(?=</Attribute>)')
    result = query.search(string)
    if result:
        return result.group()
    else:
        return "None"


def find_collection_date(string):
    query = re.compile(r'(?<="collection date">).+?(?=</Attribute>)')
    result = query.search(string)
    if result:
        return result.group()
    else:
        return "None"


def find_total_size(string):
    query = re.compile(r'(?<=total_size=")\d+?(?=")')
    result = query.search(string)
    if result:
        return result.group()
    else:
        return "None"


def find_instrument(string):
    query = re.compile(r'(?<=Instrument).+?(?="/>)')
    result = query.search(string)
    if result:
        return result.group().split('="')[1]
    else:
        return "None"


def find_library_name(string):
    query = re.compile(r"(?<=<LIBRARY_NAME>).+?(?=</LIBRARY_NAME>)")
    result = query.search(string)
    if result:
        return result.group()
    else:
        return "None"


def find_library_layout(string):
    query = re.compile(r"PAIRED|SINGLE")
    result = query.search(string)
    if result:
        return result.group()
    else:
        return "None"


def find_srr(string):
    query = re.compile(r'(?<=Run acc=").+?(?=")')
    result = query.search(string)
    if result:
        return result.group()
    else:
        return "None"


def find_spots(string):
    query = re.compile(r'(?<=total_spots=").+?(?=")')
    result = query.search(string)
    if result:
        return result.group()
    else:
        return "None"


def library_data(srs):
    handle = Entrez.esearch(db="sra", term=srs, idtype="acc")
    record = Entrez.read(handle)
    handle.close()
    library = []
    for i in record["IdList"]:
        handle = Entrez.esummary(db="sra", id=i, retmode="text")
        record = Entrez.read(handle)
        handle.close()
        total_size = find_total_size(record[0]["ExpXml"])
        instrument = find_instrument(record[0]["ExpXml"])
        library_name = find_library_name(record[0]["ExpXml"])
        library_layout = find_library_layout(record[0]["ExpXml"])
        srr = find_srr(record[0]["Runs"])
        total_spots = find_spots(record[0]["Runs"])
        library.append(
            [total_size, instrument, library_name, library_layout, srr, total_spots]
        )
    return library


def find_data(term):
    id_biosamples = find_biosamples_id(term)
    final_data = {term: defaultdict(dict)}
    for i in id_biosamples:
        handle = Entrez.esummary(db="biosample", id=i)
        record = Entrez.read(handle)
        handle.close()
        final_data[term].update({i: defaultdict(dict)})
        sample_name = find_sample_name(
            record["DocumentSummarySet"]["DocumentSummary"][0]["Identifiers"]
        )
        srs = find_srs(
            record["DocumentSummarySet"]["DocumentSummary"][0]["Identifiers"]
        )
        location = find_location(
            record["DocumentSummarySet"]["DocumentSummary"][0]["SampleData"]
        )
        lat_lon = find_lat_lon(
            record["DocumentSummarySet"]["DocumentSummary"][0]["SampleData"]
        )
        collection_date = find_collection_date(
            record["DocumentSummarySet"]["DocumentSummary"][0]["SampleData"]
        )
        final_data[term][i] = {
            "Sample name": sample_name,
            "archive sample accession": {srs: library_data(srs)},
            "location": location,
            "lat_lon": lat_lon,
            "collection date": collection_date,
        }
    return final_data


def find_lat(lat_lon):
    query = re.compile(r"[\d|.]+(?= [N|S])")
    result = query.search(lat_lon)
    if result:
        return result.group()
    else:
        return "None"


def find_lon(lat_lon):
    query = re.compile(r"[\d|.]+(?= [E|W])")
    result = query.search(lat_lon)
    if result:
        return result.group()
    else:
        return "None"


def find_md5(srr):
    r = requests.get(
        f"https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?filetype=run&location-type=forced&location=s3.us-east-1&accept-charges=aws&acc={srr}",
        auth=("user", "pass"),
    )
    query = re.compile(r'(?<="md5": ").+?(?=",)')
    result = query.search(r.text)
    if result:
        return result.group()
    else:
        return "None"


def find_link(srr):
    r = requests.get(
        f"https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?filetype=run&location-type=forced&location=s3.us-east-1&accept-charges=aws&acc={srr}",
        auth=("user", "pass"),
    )
    query = re.compile(r'(?<="link": ").+?(?=",)')
    result = query.search(r.text)
    if result:
        return result.group()
    else:
        return "None"


@click.command()
@click.option(
    "-id", "--bioproject", type=str, help="BioProject accession id", required=True
)
@click.option(
    "-o",
    "--output_file",
    type=click.Path(),
    help="output file name",
    default="table.tsv",
)
def make_table(bioproject, output_file):
    data = find_data(bioproject)
    rows = []
    for i in data[bioproject]:
        row = [bioproject]
        for j in data[bioproject][i]:
            if j == "archive sample accession":
                for z in data[bioproject][i][j]:
                    row.append(z)
                    row.extend(data[bioproject][i][j][z][0])
            else:
                row.append(data[bioproject][i][j])
        rows.append(row)
        row = [bioproject]
    table = pd.DataFrame(
        rows,
        columns=[
            "Bioproject",
            "Sample name",
            "Archive sample accession",
            "Total size",
            "Instrument",
            "Library name",
            "Library layout",
            "Archive data accession",
            "Read count",
            "Place",
            "Lat_lon",
            "Collection date",
        ],
    )
    table["Geo loc name"] = table["Place"].apply(
        lambda x: x.split(":")[0] if x != "N/A" else "N/A"
    )
    table["Site name"] = table["Place"].apply(
        lambda x: x.split(":")[1] if x != "N/A" else "N/A"
    )
    table["Collection date"] = table["Collection date"].apply(
        lambda x: x.split("-")[0] if x != "N/A" else "N/A"
    )
    table["Latitude"] = table["Lat_lon"].apply(lambda x: find_lat(x))
    table["Longitude"] = table["Lat_lon"].apply(lambda x: find_lon(x))
    table["md5"] = table["Archive data accession"].apply(lambda x: find_md5(x))
    table["link"] = table["Archive data accession"].apply(lambda x: find_link(x))
    table.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    make_table()
