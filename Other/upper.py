#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.04.2022
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru

import click


def read_file(file_name: str) -> list:
    """
    reading and transforming into upper case fna file
    """
    new_file = []
    with open(file_name) as f:
        line = f.readline()
        while line:
            if line.startswith(">"):
                new_file.append(line)
            else:
                new_file.append(line.upper())
            line = f.readline()
    return new_file


def writing_file(new_file: list, file_name: str):
    name = file_name[:-4] + "_upper.fna"
    with open(name, "w") as f:
        f.write("".join(new_file))


@click.command()
@click.argument("file_name", type=str)
def main(file_name):
    new_file = read_file(file_name)
    writing_file(new_file, file_name)


if __name__ == "__main__":
    main()
