import os
import pathlib

import Bio.AlignIO as AlignIO
import Bio.PDB as PDB
import Bio.SeqUtils as SeqUtils
import Bio.SubsMat.MatrixInfo as MatrixInfo
import Bio.pairwise2 as pairwise2
import modeller
import modeller.automodel as automodel

_verbose = False


def verbose():
    global _verbose
    _verbose = True
    modeller.log.verbose()


class _Select(PDB.Select):
    """
    Class to overwrite the selection class in Bio.PDB. Used for selecting the included chains
    during the cleaning of the modeller input structure.
    """
    def __init__(self):
        super(_Select, self).__init__()
        self.chains = []
        self.residues = []

    def extend(self, extension, ext_type):
        """
        Extends the list of the selected chains or residues.
        :param extension: List of the objects to include.
        :param ext_type: Type of the objects in the list:
        c == chain
        r == residue
        :return: None
        """
        if ext_type == "c":
            self.chains += extension
        elif ext_type == "r":
            self.residues += extension

    def accept_chain(self, chain):
        """
        Overwritten method for returning the selection.
        :param chain: Chain to be checked.
        :return: True if chain is in the selected list, else False
        """
        if chain in self.chains:
            return True
        return False

    def accept_residue(self, residue):
        """
        Overwritten method for returning the selection.
        :param residue: Residue to be checked.
        :return: True if residue iis in the selected list, else False
        """
        if residue in self.residues:
            return True
        return False


class Input:
    def __init__(self, structure_loc, target_map, structure_name="", target_name="X", master_ali=""):
        """
        Builds the input files correctly for running a model.
        :param structure_loc: Location of the template structure. If not a pathlib.Path, it gets converted to one.
        :param target_map: Tuple of used subunits; format:
        [(template_sub, target_sub), (template_sub2, target_sub2)...]
        :param structure_name: Optional name for the template structure
        (else a name is generated from structure_loc)
        :param target_name: Optional name for the target structure (else it is substituted with X)
        :param master_ali: Location of the MasterAli.pir file.
        """
        self.structure_loc = pathlib.Path(structure_loc)
        self.input_loc = self.structure_loc / "Model"
        if structure_name:
            self.structure_name = structure_name
        else:
            self.structure_name = self.structure_loc.parts[-1]
        self.target = target_map
        self.target_name = target_name

        # Checking if the structure location points to the folder or the file itself
        if self.structure_loc.is_file():
            self.structure = PDB.PDBParser(PERMISSIVE=True,
                                           QUIET=not _verbose).get_structure(self.structure_name,
                                                                             self.structure_loc)
        else:
            structure = self.structure_loc / f"{self.structure_name}.pdb"
            self.structure = PDB.PDBParser(PERMISSIVE=True,
                                           QUIET=not _verbose).get_structure(self.structure_name,
                                                                             structure)
        self.sequences = AlignIO.read(master_ali, "pir")

        # Extraction of information about chains
        self.structure_chains = {}
        self.structure_chains_reference = {}
        self.chains_extract()
        # Order of chains, used in writing the alignments, to conserve the order
        # used in the template structure file
        self.chains_ordered = [ch.id for ch in self.structure.get_chains() if ch.id in self.structure_chains.keys()]

        # Building the information about target chains in relation to template chains
        self.target_chains = {}
        self.chains_build()

        # Creating output directory
        try:
            os.mkdir(self.input_loc)
        except FileExistsError:
            pass

        # Cleaning and preparing the PDB of the template
        # Getting the start and end AA number in the cleaned template structure
        self.input_pdb_loc = self.input_loc / f"{self.structure_name.upper()}.pdb"
        self.pdb_clean()

        # Writing of the PIR/ALI file
        self.input_ali_loc = self.input_loc / f"{self.target_name}.ali"
        self.ali_write()

    def chains_extract(self):
        """
        Extracts the chain information (name, gaps and sequence) from the structure.
        :return:
        """
        builder = PDB.Polypeptide.CaPPBuilder()
        chain_names = []
        for compound in self.structure.header["compound"].values():
            if "gamma" in compound["molecule"]:
                name = compound["molecule"].split(",")[0].split(" ")[-1].capitalize()
                chain_names += [(chain, name) for chain in compound["chain"].upper().split(", ")]

        for chain, sub_name in chain_names:
            asequence = builder.build_peptides(self.structure[0][chain])
            asequence = asequence[0].get_sequence()

            ref_seq = [seq for seq in self.sequences if sub_name in seq.name][0]
            alignment = pairwise2.align.globalds(asequence.ungap("-"),
                                                 ref_seq.seq.ungap("-"),
                                                 MatrixInfo.blosum62, -10, -0.5,
                                                 one_alignment_only=True)
            self.structure_chains[chain] = {"name": sub_name,
                                            "gaps": [i for i, aa in enumerate(alignment[0][1]) if aa == "-"],
                                            "sequence": alignment[0][0]}

            self.structure_chains_reference[chain] = {"name": sub_name,
                                                      "gaps": [i for i,
                                                               aa in enumerate(alignment[0][1]) if aa == "-"],
                                                      "sequence": alignment[0][0]}

    def chains_build(self):
        """
        Building/preparing the target sequence chains.
        :return: None
        """
        # Transforming the mapping of the input to the mapping in the template file
        input_mapping = [x[0] for x in self.target]
        input_mapping_target = [x[1] for x in self.target]
        template = [self.structure_chains[chain]["name"] for chain in self.chains_ordered]
        normal_orientation = [template] + [template[i:]+template[:i] for i in range(1, len(template))]
        # Checking the orientation of the input
        if input_mapping not in normal_orientation:
            input_mapping = input_mapping[::-1]
            input_mapping_target = input_mapping_target[::-1]
        # Checking the phase of the input
        while input_mapping != template:
            input_mapping = input_mapping[1:] + input_mapping[:1]
            input_mapping_target = input_mapping_target[1:] + input_mapping_target[:1]

        for target_subunit, chain in zip(input_mapping_target, self.chains_ordered):
            data = self.structure_chains[chain]
            self.structure_chains[chain]["sequence"], target_ali = self.align(data["sequence"],
                                                                              f"{self.structure_name.upper()}",
                                                                              target_subunit,
                                                                              data["gaps"])
            self.target_chains[chain] = {"name": target_subunit,
                                         "sequence": target_ali}

    def align(self, seq1, seq1_name, seq2_name, breaks):
        """
        Aligns the sequences and creates the chain brakes at the appropriate places.
        The flanking extra amino acids of seq2 are also removed.
        :param seq1: First sequence (as string); this sequence is intended as a reference
        from a structure
        :param seq1_name: Name of the first sequence
        :param seq2_name: Name of the second sequence
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
        break_sequences = [seq1[sequence[0]: sequence[-1] + 1] for sequence in break_blocks]
        if breaks[0] == 0:
            break_sequences = break_sequences[1:]
        if seq1.endswith(break_sequences[-1]):
            break_sequences = break_sequences[:-1]

        # Getting the aligned sequences
        seq1_out = str([x for x in self.sequences if
                        (seq1_name in x.name and seq2_name[0].lower() == x.name[-2])][0].seq)
        seq2_out = str([x for x in self.sequences if seq2_name in x.name][0].seq)
        # Removing gaps present in both sequences
        gaps = [i for i, (s1, s2) in enumerate(zip(seq1_out, seq2_out)) if s1 == s2 == "-"]
        for gap in gaps[::-1]:
            seq1_out = seq1_out[:gap] + seq1_out[gap + 1:]
            seq2_out = seq2_out[:gap] + seq2_out[gap + 1:]

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

    def ali_write(self):
        """
        Writes the ALI file, needed as an input to the modeller.
        >P1;name
        structure/sequence:name:beginning_residue:beginning_chain:ending_residue:ending_chain:name:::
        :return: None
        """
        # Sorted because the modeller reads it so.
        seq1_out = []
        seq2_out = []
        for chain in self.chains_ordered:
            seq1_out.append("\n".join([self.structure_chains[chain]["sequence"][i:i+50] for
                                       i in range(0, len(self.structure_chains[chain]["sequence"]), 50)]))
            seq2_out.append("\n".join([self.target_chains[chain]["sequence"][i:i+50] for
                                       i in range(0, len(self.target_chains[chain]["sequence"]), 50)]))
        seq1_out = "\n".join([i+"/" if len(i) % 50 != 0 else i+"\n/" for i in seq1_out])
        seq2_out = "\n".join([i+"/" if len(i) % 50 != 0 else i+"\n/" for i in seq2_out])

        start = self.structure_chains[self.chains_ordered[0]]['start']
        end = self.structure_chains[self.chains_ordered[-1]]['end']
        out = "\n".join([f">P1;{self.structure_name.upper()}",
                         "structure:{0}:{1}:{2}:{3}:{4}:{0}:::".format(self.structure_name.upper(),
                                                                       start,
                                                                       self.chains_ordered[0],
                                                                       end,
                                                                       self.chains_ordered[-1]),
                         seq1_out,
                         "*",
                         "",
                         f">P1;{self.target_name}",
                         f"sequence:{self.target_name}:1::::{self.target_name}:::",
                         seq2_out,
                         "*"])
        with open(self.input_ali_loc, "w") as file:
            file.write(out)

    def pdb_clean(self):
        """
        Prepares the template pdb file for the modeling.
        Removes the unnecessary chains and extra molecules and saves the structure.
        :return: None
        """
        io = PDB.PDBIO()
        io.set_structure(self.structure)
        selected = _Select()
        selected.extend([x for x in self.structure[0].get_list() if x.id in self.structure_chains.keys()], "c")

        residues = {ch.id: list(ch.get_residues()) for ch in self.structure[0].get_list()
                    if ch.id in self.structure_chains.keys()}
        accepted_res = []
        for chain, sequence in residues.items():
            while SeqUtils.seq1("".join(aa.get_resname() for
                                        aa in sequence[:5])) != self.structure_chains[chain]["sequence"][:5]:
                sequence = sequence[1:]
            while SeqUtils.seq1("".join(aa.get_resname() for
                                        aa in sequence[-5:])) != self.structure_chains[chain]["sequence"][-5:]:
                sequence = sequence[:-1]
            accepted_res = accepted_res + sequence
            residues[chain] = sequence
        selected.extend(accepted_res, "r")

        io.save(str(self.input_pdb_loc.absolute()), selected)

        # Adding the beginning and ending AAs to the list
        for chain in self.structure_chains.keys():
            self.structure_chains[chain]["start"] = residues[chain][0].full_id[-1][1]
            self.structure_chains[chain]["end"] = residues[chain][-1].full_id[-1][1]


class Model:
    def __init__(self, settings, starting_model=5, ending_model=5):
        """
        The model creation object.
        :param settings: Input object created for the input settings to the model.
        """
        self.env = modeller.environ()

        # Changing the working directory
        if settings.input_loc.is_file():
            os.chdir(settings.input_loc.parent)
        else:
            os.chdir(settings.input_loc)
        self.env.io.atom_files_directory = ["./"]

        self.alignment = modeller.automodel.automodel(self.env,
                                                      alnfile=settings.input_ali_loc.parts[-1],
                                                      knowns=settings.structure_name.upper(),
                                                      sequence=settings.target_name)

        self.alignment.starting_model = starting_model
        self.alignment.ending_model = ending_model

        self.alignment.make()

        # Rearranging and renaming of the chains
        subchains_count = {ch: len(settings.target_chains[ch]["sequence"].split("/")) for ch in settings.chains_ordered}
        self.created_files = [x for x in os.listdir(".") if x.startswith(settings.target_name) and x.endswith(".pdb")]
        for file in self.created_files:
            structure = PDB.PDBParser(PERMISSIVE=True,
                                      QUIET=not _verbose).get_structure(settings.target_name,
                                                                        file)
            chains = list(structure[0].get_chains())
            # Rename all chains so they can be named correctly in the end
            for i, chain in enumerate(chains):
                chain.id = f"CHAIN{i}"
            # Extend the chains with the initial AAs with AAs from the same subunits
            i = 0
            for chain in settings.chains_ordered:
                base = chains[i]
                base.id = chain
                i += 1
                for _ in range(subchains_count[chain]-1):
                    for residue in chains[i]:
                        base.add(residue)
                    i += 1
            # Remove extra chains
            for chain in chains:
                if chain.id.startswith("CHAIN"):
                    del structure[0][chain.id]

            for chain in structure[0].get_chains():
                # sequence_target = settings.target_chains[chain.id]["sequence"]
                sequence_target = settings.target_chains[chain.id]["sequence"].replace("/", "")
                # sequence_template = settings.structure_chains[chain.id]["sequence"]
                sequence_template = settings.structure_chains[chain.id]["sequence"].replace("/", "")
                chain_target = list(structure[0][chain.id].get_residues())
                chain_template = list(settings.structure[0][chain.id].get_residues())

                # Changing the numbering so to not introduce unwanted clashes in the assignment of the real numbers
                j = 0
                chain_target_dict = {}
                for i, residue_target in enumerate(sequence_target):
                    if residue_target != "-":
                        while SeqUtils.seq1(chain_target[j].get_resname()) != residue_target:
                            j += 1
                        chain_target[j].id = (" ", 10000+i, " ")
                        chain_target_dict[10_000+i] = chain_target[j]
                        j += 1

                # Introducing the numbering from the template
                j = 0
                for i, (residue_target, residue_template) in enumerate(zip(sequence_target, sequence_template)):
                    if residue_target != "-":
                        while SeqUtils.seq1(chain_template[j].get_resname()) != residue_template:
                            j += 1
                        chain_target_dict[10_000+i].id = chain_template[j].id
                        j += 1

            # Save back into the file
            io = PDB.PDBIO()
            io.set_structure(structure)
            io.save(file)


if __name__ == "__main__":
    # Variables for the input

    # PIR file (.ali) contains the sequences of the target protein in the format:
    # >P1;name
    # sequence:model_name
    # AA SEQUENCE*
    _pir_input = pathlib.Path(".") / "MasterAli.pir"
    _structure_loc = pathlib.Path(".") / "Structures"
    _used_structure_name = "6a96"
    _used_structure_loc = _structure_loc / f"{_used_structure_name}.pdb"
    _target_name = "a6b3"
    _target = [("Alpha-5", "Alpha-6"),
               ("Beta-3", "Beta-3"),
               ("Beta-3", "Beta-3"),
               ("Beta-3", "Beta-3"),
               ("Beta-3", "Beta-3")]

    _model_input = Input(_structure_loc,
                         _target,
                         _used_structure_name,
                         _target_name, _pir_input)

    _model = Model(_model_input)
