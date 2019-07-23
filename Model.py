import modeller
import Bio.AlignIO as AlignIO
import Bio.PDB.PDBParser as PDBParser
import Bio.PDB.Polypeptide as Polypeptide
import Bio.pairwise2 as pairwise2
import Bio.SubsMat.MatrixInfo as MatrixInfo
import Bio.Seq
import Bio.PDB as PDB


class _Select(PDB.Select):
    """
    Class to overwrite the selection class in Bio.PDB. Used for selecting the included chains
    during the cleaning of the modeller input structure.
    """
    def __init__(self, selected_chains):
        super(_Select, self).__init__()
        self.chains = selected_chains

    def accept_chain(self, chain):
        """
        Overwritten method for returning the selection.
        :param chain: Chain to be checked.
        :return: True if chain is in the selected list, else False
        """
        if chain in self.chains:
            return True
        return False


def align(seq1, seq2, seq1_name, seq2_name, breaks):
    """
    Aligns the sequences and creates the chain brakes at the appropriate places.
    The flanking extra amino acids of seq2 are also removed.
    :param seq1: First sequence (as string); this sequence is intended as a reference
    from a structure
    :param seq2: Second sequence (as string)
    :param breaks: A list of indices of brakes, adjusted to the seq1
    :return: Tuple of two sequences as strings
    """
    # Splitting the empty indices into coherent blocks and removing the first and
    # last block - flanking regions.
    break_blocks = [[breaks[0]]]
    for index in breaks[1:]:
        if index == break_blocks[-1][-1] + 1:
            break_blocks[-1].append(index)
        else:
            break_blocks.append([index])
    # Extracting the sequences of the blocks
    break_sequences = [seq1[sequence[0]: sequence[-1]+1] for sequence in break_blocks]
    if breaks[0] == 0:
        break_sequences = break_sequences[1:]
    if seq1.endswith(break_sequences[-1]):
        break_sequences = break_sequences[:-1]

    # Getting the aligned sequences
    seq1_out = str([x for x in sequences if (seq1_name in x.name and seq2_name[0].lower() == x.name[-2])][0].seq)
    seq2_out = str([x for x in sequences if seq2_name in x.name][0].seq)
    # Removing gaps present in both sequences
    gaps = [i for i, (s1, s2) in enumerate(zip(seq1_out, seq2_out)) if s1 == s2 == "-"]
    for gap in gaps[::-1]:
        seq1_out = seq1_out[:gap] + seq1_out[gap+1:]
        seq2_out = seq2_out[:gap] + seq2_out[gap+1:]

    # Inserting "/" as a chain break sign at the end of the break blocks (this is why the chains
    # are operated on in the reversed state.
    seq1_out, seq2_out = seq1_out[::-1], seq2_out[::-1]
    for break_sequence in break_sequences:
        i = seq1_out.index(break_sequence[::-1])
        seq1_out = seq1_out[:i] + "/" + seq1_out[i:]
        seq2_out = seq2_out[:i] + "/" + seq2_out[i:]
    # Removal of the flanking sequences, which are not present in the structure
    while seq1_out.startswith("-"):
        seq1_out, seq2_out = seq1_out[1:], seq2_out[1:]
    while seq1_out.endswith("-"):
        seq1_out, seq2_out = seq1_out[:-1], seq2_out[:-1]

    return seq1_out[::-1], seq2_out[::-1]


def pir_write(data, filename):
    """
    Writes the PIR file, needed as an input to the modeller.
    :param data: A list of two dictionaries with sequences.
    First should be the template, the second the target structure.
    :param filename: Name of the file, to which it saves.
    :return: True if successful else False
    """
    try:
        # Sorted for possible better viewing of the generated file by humans
        out = []
        for chain in sorted(data[0].keys()):
            out.append("\n".join([f">P1;{used_structure_name}",
                                  f"structure:{used_structure_name}:{data[0][chain]['start']}:{chain}:{data[0][chain]['end']}:{chain}:{used_structure_name}:::",
                                  *[data[0][chain]["sequence"][i:i+51] for i in range(0, len(data[0][chain]["sequence"]), 50)],
                                  "*"]))

            out.append("\n".join([f">P1;{target_name}",
                                  f"sequence:{target_name}:1:{chain}::{chain}:{target_name}:::",
                                  *[data[1][chain]["sequence"][i:i+51] for i in range(0, len(data[1][chain]["sequence"]), 50)],
                                  "*"]))

        with open(filename, "w") as file:
            file.write("\n".join(out))

        return True
    except Exception as e:
        print(e)
        return False


if __name__ == "__main__":
    # Variables for the input

    # PIR file (.ali) contains the sequences of the target protein in the format:
    # >P1;name
    # sequence:model_name
    # AA SEQUENCE*
    pir_input = "C:/Users/Jure/Documents/Jure files/GitHub/ligand_dock_project/MasterAli.pir"
    structure_loc = "C:/Users/Jure/Documents/Jure files/GitHub/ligand_dock_project/Structures/"
    used_structure_name = "6a96"
    used_structure_loc = f"{structure_loc}{used_structure_name}.pdb"
    target = ("Alpha-6", "Beta-3")
    target_name = "a6b3"

    structure_parser = PDBParser(PERMISSIVE=True)
    structure = structure_parser.get_structure(used_structure_name, used_structure_loc)
    compounds = structure.header["compound"]

    sequences = AlignIO.read(pir_input, "pir")
    sequences_names = [x.name for x in sequences]

    # Extraction of the chain names
    chain_names = [(compound["molecule"].split(",")[0].split(" ")[-1].capitalize(),
                    [chain.capitalize() for chain in compound["chain"].split(", ")]) for compound in compounds.values() if "gamma" in compound["molecule"]]

    # Aligning the chains to the reference sequences and getting the gap/linker indices
    chain_reference = {}
    builder = Polypeptide.CaPPBuilder()
    for sub_name, chains in chain_names:
        for chain in chains:
            asequence = builder.build_peptides(structure[0][chain])
            asequence = asequence[0].get_sequence()
            ref_seq = sequences[sequences_names.index(sub_name)]
            alignment = pairwise2.align.globalds(asequence.ungap("-"),
                                                 ref_seq.seq.ungap("-"),
                                                 MatrixInfo.blosum62, -10, -0.5,
                                                 one_alignment_only=True)
            chain_reference[chain] = {"name": sub_name,
                                      "gaps": [i for i, aa in enumerate(alignment[0][1]) if aa == "-"],
                                      "sequence": alignment[0][0]}

    # Aligning the target sequences to the chains
    # Breaks need to be inserted (by insertion of "/")
    target_reference = {}
    for chain, data in chain_reference.items():
        target_subunit = [x for x in target if data["name"].split("-")[0] in x][0]
        chain_reference[chain]["sequence"], target_ali = align(data["sequence"],
                                                               [str(x.seq) for x in sequences if x.name == target_subunit][0],
                                                               f"{used_structure_name.upper()}",
                                                               target_subunit,
                                                               data["gaps"])
        target_reference[chain] = {"name": target_subunit,
                                   "sequence": target_ali}

    # Getting the start and end AA number in the template structure (for now a very brute force approach)
    for chain in chain_reference.keys():
        aa = list(structure[0][chain].get_residues())
        chain_reference[chain]["start"] = aa[0].full_id[-1][1]
        chain_reference[chain]["end"] = aa[-1].full_id[-1][1]

    # Writing of the PIR file; format:
    # >P1;name
    # structure/sequence:name:beginning_residue:beginning_chain:ending_residue:ending_chain:name:::
    # SEQ
    # *
    pir_write([chain_reference, target_reference], f"{structure_loc}{target_name}.pir")

    # Structure cleaning - removal of unused chains (extra molecules etc.)
    # coupled with structure saving
    io = PDB.PDBIO()
    io.set_structure(structure)
    selected_chains = _Select([x for x in structure[0].get_list() if x.id in chain_reference.keys()])
    io.save(f"{structure_loc}{used_structure_name}_cleaned.pdb", selected_chains)

    # modeller.log.verbose()
    # environment = modeller.environ()

    # Database of sequence alignments reading/writing
    # Modeller can create a binary DB, which works faster  than reading PIR file.
    # We should opt for this once we create the alignment master-PIR.

    # Alignment object creation
    # alignment = modeller.alignment(environment)
