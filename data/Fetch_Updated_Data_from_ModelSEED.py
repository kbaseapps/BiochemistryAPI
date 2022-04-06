#!/usr/bin/env python
from urllib.request import urlopen
import json

MSD_git_url = "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/"
MSD_commit = "master"

#
# Compounds
#
file = urlopen(MSD_git_url + MSD_commit + "/Biochemistry/compounds.json")
compounds = json.load(file)

compound_headers = [
    "id",
    "abbreviation",
    "name",
    "formula",
    "mass",
    "source",
    "structure",
    "charge",
    "is_core",
    "is_obsolete",
    "is_cofactor",
    "deltaG",
    "deltaGErr",
    "pka",
    "pkb",
    "aliases",
]

cpd_header_translation = {
    "deltaG": "deltag",
    "deltaGErr": "deltagerr",
    "structure": "inchikey",
}

with open("compounds.tsv", "w") as cpd_fh:
    cpd_fh.write("\t".join(compound_headers) + "\n")

    for cpd in compounds:
        line_values = list()
        for field in compound_headers:
            true_field = field
            if field not in cpd and field in cpd_header_translation:
                true_field = cpd_header_translation[field]

            line_values.append(str(cpd[true_field]))

        cpd_fh.write("\t".join(line_values) + "\n")
cpd_fh.close()

#
# Reactions
#
file = urlopen(MSD_git_url + MSD_commit + "/Biochemistry/reactions.json")
reactions = json.load(file)

reaction_headers = [
    "id",
    "abbreviation",
    "name",
    "code",
    "stoichiometry",
    "is_transport",
    "equation",
    "definition",
    "reversibility",
    "direction",
    "abstract_reaction",
    "pathways",
    "aliases",
    "ec_numbers",
    "deltag",
    "deltagerr",
    "compound_ids",
    "status",
    "is_obsolete",
    "linked_reaction",
]

rxn_header_translation = {}

with open("reactions.tsv", "w") as rxn_fh:
    rxn_fh.write("\t".join(reaction_headers) + "\n")

    for rxn in reactions:
        line_values = list()
        for field in reaction_headers:
            true_field = field
            if field not in rxn and field in rxn_header_translation:
                true_field = rxn_header_translation[field]

            line_values.append(str(rxn[true_field]))

        rxn_fh.write("\t".join(line_values) + "\n")
rxn_fh.close()

#
# Aliases
#
alias_headers = ["MS ID", "Old MS ID", "External ID", "Source"]

with open("Compounds_Aliases.tsv", "w") as cpda_fh:
    cpda_fh.write("\t".join(alias_headers) + "\n")

    file = urlopen(
        MSD_git_url
        + MSD_commit
        + "/Biochemistry/Aliases/Unique_ModelSEED_Compound_Aliases.txt"
    )
    header = 1
    for line in file.readlines():
        if header == 1:
            header -= 1
            continue

        line = line.decode("utf-8").strip()
        array = line.split("\t")

        cpda_fh.write("\t".join([array[0], "", array[1], array[2]]) + "\n")

with open("Reactions_Aliases.tsv", "w") as rxna_fh:
    rxna_fh.write("\t".join(alias_headers) + "\n")

    file = urlopen(
        MSD_git_url
        + MSD_commit
        + "/Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt"
    )
    header = 1
    for line in file.readlines():
        if header == 1:
            header -= 1
            continue

        line = line.decode("utf-8").strip()
        array = line.split("\t")

        rxna_fh.write("\t".join([array[0], "", array[1], array[2]]) + "\n")

with open("Enzyme_Class_Reactions_Aliases.tsv", "w") as rxnec_fh:
    rxnec_fh.write("\t".join(alias_headers) + "\n")

    file = urlopen(
        MSD_git_url
        + MSD_commit
        + "/Biochemistry/Aliases/Unique_ModelSEED_Reaction_ECs.txt"
    )
    header = 1
    for line in file.readlines():
        if header == 1:
            header -= 1
            continue

        line = line.decode("utf-8").strip()
        array = line.split("\t")

        rxnec_fh.write("\t".join([array[0], "", array[1], array[2]]) + "\n")
