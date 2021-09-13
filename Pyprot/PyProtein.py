from urllib.request import urlopen
from urllib.parse import urlencode
from pandas import DataFrame
from Bio import SeqIO
import os


def search(text, output_format="tab", sort="score", organism="", columns=(),
           isoform=False, compress=False, offset=0, limit=0):
    """
    Perform a query over the UniProt API.
    :param text: text to be searched
    :param output_format: html, tab, xls, fasta, gff, txt, xml, rdf, list or rss
           Currently only tab format is supported and converted to DataFrame
    :param sort: by which parameter the results are sorted
    :param organism: specify the organism (optional)
    :param columns: list/tuple of strings from the following options:
           citation, clusters, comments, domains, domain, ec, id, entry name,
           existence, families, features, genes, go, go-id, interactor,7
           keywords, last-modified, length, organism, organism-id, pathway,
           protein names, reviewed, sequence, 3d, version, virus hosts
    :param isoform: should isoforms be considered
    :param compress: should the file be compressed
    :param offset:
    :param limit:
    :return: DataFrame or urlopen response
    More at: https://www.uniprot.org/help/api_queries
    """
    cgi = "https://www.uniprot.org/uniprot/?"
    variables = {"query": text,
                 "format": output_format,
                 "sort": sort,
                 "offset": str(offset)}
    if organism:
        variables["organism"] = organism
    if columns:
        variables["columns"] = ",".join(columns)
    if isoform:
        variables["isoform"] = "Yes"
    if compress:
        variables["compress"] = "Yes"
    if limit:
        variables["limit"] = str(limit)

    fullcgi = "".join((cgi, urlencode(variables)))
    if output_format == "tab":
        response = urlopen(fullcgi)
        contents = [x.split("\t") for x in response.read().decode(response.headers.get_content_charset()).split("\n")]
        return DataFrame(contents[1:], columns=contents[0])
    else:
        return urlopen(fullcgi)


def seq_download(name, organism="Homo sapiens", gaba=False):
    """
    Downloads the sequence from UniProt. Returns Seq object or -1 if sequence not found.
    :param name: Name of the searched protein
    :param organism: Organism the protein comes from.
    :param gaba: Is it a reference gaba subunit (the inner subunit naming reference is then used). If
    False, then the name itself is searched.
    :return: Seq object containing all the info from fasta download of the sequence. If sequence not found,
    then -1 is returned.
    """

    subunits = {
        "Alpha-1": "Gabra1",
        "Alpha-2": "Gabra2",
        "Alpha-3": "Gabra3",
        "Alpha-4": "Gabra4",
        "Alpha-5": "Gabra5",
        "Alpha-6": "Gabra6",
        "Beta-1": "Gabrb1",
        "Beta-2": "Gabrb2",
        "Beta-3": "Gabrb3",
        "Gamma-1": "Gabrg1",
        "Gamma-2": "Gabrg2",
        "Gamma-3": "Gabrg3",
        "Delta": "Gabrd",
        "Pi": "Gabrp",
        "Rho-1": "Gabrr1",
        "Rho-2": "Gabrr2",
        "Rho-3": "Gabrr3",
        "Epsilon": "Gabre",
        "Theta": "Gabrq"
    }
    if gaba:
        results = search(subunits[name])
    else:
        results = search(name)
    results = results[results["Organism"].str.contains(organism, na=False)]
    if len(results):
        if gaba:
            target = results[results["Gene names"].str.contains(subunits[name].upper())]["Entry"].max()
        else:
            target = results[results["Gene names"].str.contains(name)]["Entry"].max()
        response = urlopen(f"https://www.uniprot.org/uniprot/{target}.fasta").read().decode("utf-8")
        with open("Temp_seq.fasta", "w") as file:
            file.write(response)
        seq = SeqIO.read("Temp_seq.fasta", "fasta")
        os.remove("Temp_seq.fasta")

        return seq

    else:
        return -1


if __name__ == "__main__":
    _subunit = "Alpha-6"
    _organism = "Homo sapiens"
    seq_download(_subunit, _organism)
