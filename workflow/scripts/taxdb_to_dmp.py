#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert a KrakenUniq ``taxDB`` into NCBI ``names.dmp`` + ``nodes.dmp``.

ngsLCA needs NCBI-style taxonomy dumps, but a KrakenUniq database ships only a
compact ``taxDB`` (tab-separated: taxid, parent_taxid, name, rank). This emits the
minimal dmp columns ngsLCA reads:

  nodes.dmp : tax_id | parent_tax_id | rank | ...
  names.dmp : tax_id | name_txt | | scientific name |

Fields are separated by the NCBI ``\\t|\\t`` delimiter and each line ends ``\\t|``.

Usage: taxdb_to_dmp.py TAXDB NAMES_DMP NODES_DMP
"""
import sys


def main() -> None:
    if len(sys.argv) != 4:
        sys.exit("usage: taxdb_to_dmp.py TAXDB NAMES_DMP NODES_DMP")
    taxdb, names_out, nodes_out = sys.argv[1:4]

    n = 0
    with open(taxdb) as fh, open(names_out, "w") as names, open(nodes_out, "w") as nodes:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            taxid, parent, name, rank = parts[0], parts[1], parts[2], parts[3]
            nodes.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\n")
            names.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")
            n += 1
    print(f"Converted {n} taxa from {taxdb}", file=sys.stderr)


if __name__ == "__main__":
    main()
