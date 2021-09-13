import os
import pathlib

import Bio.AlignIO as AlignIO
import Bio.PDB as PDB
import Bio.SeqUtils as SeqUtils
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
        :return: True if residue is in the selected list, else False
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
        # Setting a new directory as a location
        i = 1
        while os.path.exists(self.structure_loc / "Model_{}_{:0>2d}".format(target_name, i)):
            i += 1
        self.input_loc = self.structure_loc / "Model_{}_{:0>2d}".format(target_name, i)
        # Creating output directory
        try:
            os.mkdir(self.input_loc)
        except FileExistsError:
            pass

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

        self.input_pdb_loc = self.input_loc / f"{self.structure_name.upper()}.pdb"  # Where the cleaned input structure will be saved
        self.input_ali_loc = self.input_loc / f"{self.target_name}.ali" # The path to the model-speciffic alignment

        # Extraction of information about chains
        self.chains_ordered = []
        self.structure_chains = {}
        self.structure_chains_reference = {}
        self.target_chains = {}

    def chain_preparation(self):
        """
        Extracts the  chains from the structure and builds the reference sequences.
        Cleans the sequences where necessary (gaps in both sequences and extra AAs in
        the beggining and end of the target sequence).
        Inserts '/' where the chain brake is located.
        :return:
        """
        # Extracting the chain names from the structure header
        chain_names = {}
        for compound in self.structure.header["compound"].values():
            if "gamma" in compound["molecule"]:
                name = compound["molecule"].split(",")[0].split(" ")[-1].capitalize()
                for chain in compound["chain"].upper().split(", "):
                    chain_names[chain] = name

        # getting the order if the chains
        self.chains_ordered = [ch.id for ch in self.structure.get_chains() if ch.id in chain_names.keys()]

        # Shifting the input mapping to achieve the correct substitutions
        input_mapping = [x[0] for x in self.target]
        input_mapping_target = [x[1] for x in self.target]
        template = [chain_names[ch] for ch in self.chains_ordered]
        normal_orientation = [template] + [template[i:] + template[:i] for i in range(1, len(template))]
        # Checking the orientation of the input
        if input_mapping not in normal_orientation:
            input_mapping = input_mapping[::-1]
            input_mapping_target = input_mapping_target[::-1]
        # Checking the phase of the input
        while input_mapping != template:
            input_mapping = input_mapping[1:] + input_mapping[:1]
            input_mapping_target = input_mapping_target[1:] + input_mapping_target[:1]

        # Fetching the sequences
        # template_seq = {"subunit": Seq()}
        template_seq = {seq.name.split("_")[0].split(".")[-1]: seq for
                        seq in self.sequences if seq.name.startswith(self.structure_name.upper())}
        target_seq = {seq.name: seq for seq in self.sequences if seq.name in input_mapping_target}
        target_seq = {ch: target_seq[name] for ch, name in zip(self.chains_ordered, input_mapping_target)}

        # Modifying all the gaps
        for ch in self.chains_ordered:
            # Cleaning the extra gaps in both sequences
            seq_1 = template_seq[ch]
            seq_2 = target_seq[ch]
            extra_i = [i for i, (aa1, aa2) in enumerate(zip(seq_1, seq_2)) if aa1 == "-" and aa2 == "-"]
            for i in extra_i[::-1]:
                seq_1 = seq_1[:i] + seq_1[i+1:]
                seq_2 = seq_2[:i] + seq_2[i+1:]

            # Correct the beginning and ending of the sequence to the structure (the MasterAli.pir already has
            # the correct AAs for the structures, which means we need to remove extra AAs on both sides of the
            # target sequence only)
            while seq_1[0] == "-" and seq_2[0] != "-":
                seq_1 = seq_1[1:]
                seq_2 = seq_2[1:]
            while seq_1[-1] == "-" and seq_2[-1] != "-":
                seq_1 = seq_1[:-1]
                seq_2 = seq_2[:-1]

            # Add extra periods to sequences for non-AA residues at the end of the chains in the structure
            for res in [SeqUtils.seq1(r.get_resname()) for r in self.structure[0][ch].get_residues()][::-1]:
                if res == "X":
                    seq_1 += "."
                    seq_2 += "."
                else:
                    break

            # Getting the chain breaks from the structure
            # The indices of gaps are the indices on which to break with python syntax (ie. Seq[:i] + "/" + Seq[i:])
            gaps_temporary = []
            for aa1, aa2, aa3, aa4, aa5 in zip(list(self.structure[0][ch].get_residues())[0:],
                                               list(self.structure[0][ch].get_residues())[1:],
                                               list(self.structure[0][ch].get_residues())[2:],
                                               list(self.structure[0][ch].get_residues())[3:],
                                               list(self.structure[0][ch].get_residues())[4:]):
                # Check if there is a jump in residue numbering and if the next residue is an AA
                # !! To implement: add the non-AA residues to the sequences and remove the AA checking !!
                if aa1.id[1] != aa2.id[1]-1 and SeqUtils.seq1(aa2.get_resname()) != "X":
                    # Find the index of the post-gap sequence
                    i = seq_1.seq.find("".join([SeqUtils.seq1(aa.get_resname()) for aa in (aa2, aa3, aa4, aa5)]))
                    # Shift the index back if there are multiple gaps on the place of
                    while seq_1[i-1] == "-":
                        i -= 1
                    gaps_temporary.append(i)
            # Inserting the gaps into the sequences (in reverse order to not change the index of the subsequent
            # gap indices)
            for i in gaps_temporary[::-1]:
                seq_1 = seq_1[:i] + "/" + seq_1[i:]
                seq_2 = seq_2[:i] + "/" + seq_2[i:]

            self.structure_chains[ch] = {"name": input_mapping[self.chains_ordered.index(ch)],
                                         "gaps": gaps_temporary,
                                         "sequence": str(seq_1.seq)}

            self.target_chains[ch] = {"name": input_mapping[self.chains_ordered.index(ch)],
                                      "sequence": str(seq_2.seq)}

    def ali_write(self):
        """
        Writes the ALI file, needed as an input to the modeller.
        >P1;name
        structure/sequence:name:beginning_residue:beginning_chain:ending_residue:ending_chain:name:::
        :return: None
        """
        # Ordered because the modeller reads it so.
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
        selected.extend([res for ch in residues.values() for res in ch], "r")

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
        self.settings = settings

        # Changing the working directory
        if settings.input_loc.is_file():
            os.chdir(settings.input_loc.parent)
        else:
            os.chdir(settings.input_loc)
        self.env.io.atom_files_directory = ["./"]

        # Check if there are non-protein residues in the alignment and let modeller know he should consider them
        for chain in settings.structure_chains.values():
            if chain["sequence"].endswith("."):
                self.env.io.hetatm = True
                break

        self.alignment = modeller.automodel.automodel(self.env,
                                                      alnfile=settings.input_ali_loc.parts[-1],
                                                      knowns=settings.structure_name.upper(),
                                                      sequence=settings.target_name)

        self.alignment.starting_model = starting_model
        self.alignment.ending_model = ending_model

    def run_model(self):
        """
        Method that runs the modelling.
        :param: None
        :return: None
        """
        self.alignment.make()

    def chain_cleanup(self):
        """
        Cleans up the created models and saves the cleaned versions into the read files.
        Cleaning up consists of renaming the chains the same way as in the template and renumbering the residues to the template numbers.
        :param: None
        :return: None
        """
        # Rearranging and renaming of the chains
        subchains_count = {ch: len(self.settings.target_chains[ch]["sequence"].split("/")) for ch in self.settings.chains_ordered}
        self.created_files = [x for x in os.listdir(".") if x.startswith(self.settings.target_name) and x.endswith(".pdb")]
        for file in self.created_files:
            structure = PDB.PDBParser(PERMISSIVE=True,
                                      QUIET=not _verbose).get_structure(self.settings.target_name,
                                                                        file)
            chains = list(structure[0].get_chains())
            # Rename all chains so they can be named correctly in the end
            for i, chain in enumerate(chains):
                chain.id = f"CHAIN{i}"
            # Extend the chains with the initial AAs with AAs from the same subunits
            i = 0
            for chain in self.settings.chains_ordered:
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
                # Fetch all resources needed for renumbering
                chain_target = list(structure[0][chain.id].get_residues())
                chain_template = list(self.settings.structure[0][chain.id].get_residues())

                # Construct aligned lists of chain residues
                chain_target_aligned = []
                i = 0
                for seq in self.settings.target_chains[chain.id]["sequence"].replace("/", ""):
                    if seq == "-":
                        chain_target_aligned.append(seq)
                    else:
                        chain_target_aligned.append(chain_target[i])
                        i += 1
                chain_template_aligned = []
                i = 0
                for seq in self.settings.structure_chains[chain.id]["sequence"].replace("/", ""):
                    if seq == "-":
                        chain_template_aligned.append(seq)
                    else:
                        chain_template_aligned.append(chain_template[i])
                        i += 1

                # Change the numbering so to not introduce unwanted clashes in the assignment of the real numbers
                for residue in chain_target_aligned:
                    if not type(residue) == str:
                        # residue is not "-"; the number is increased by 10.000
                        residue.id = (" ", 10_000+residue.id[1], " ")
                # Introduce the correct numbering
                for residue_template, residue_target in zip(chain_template_aligned, chain_target_aligned):
                    if not (type(residue_target) == str or type(residue_template) == str):
                        residue_target.id = residue_template.id
                # Checking which AAs in the target don't have the correct numbering and renumber them
                for i, residue in enumerate(chain_target_aligned):
                    if not type(residue) == str:
                        if residue.id[1] > 999:
                            # Find the number before and after that is not high
                            x = 1
                            while type(chain_target_aligned[i-x]) == str or chain_target_aligned[i-x].id[1] > 999:
                                x += 1
                            y = 1
                            # Accounting for the AA being the last on the list
                            try:
                                chain_target_aligned[i+y]
                            except IndexError:
                                # There is no residue after the current one, so the number is just one higher than
                                # the previous residues
                                residue.id = (" ", chain_target_aligned[i-x].id[1]+1, " ")
                            else:
                                # There are AAs after the one being renumbered
                                while type(chain_target_aligned[i+y]) == str or chain_target_aligned[i+y].id[1] > 999:
                                    if i+y == len(chain_target_aligned)-1:
                                        break
                                    y += 1
                                if chain_target_aligned[i-x].id[1] - chain_target_aligned[i+y].id[1] == 1:
                                    # Adding 10.000 to the number before  and adding x (distance in the alignment), to
                                    # avoid multiple consecutive AAs having the same number
                                    residue.id = (" ", chain_target_aligned[i-x].id[1]+10_000+x, " ")
                                else:
                                    # Adding 1 to the previous AA
                                    # Works if there is a gap in numbering in the template, or if at the end of
                                    # the chain no AA has a number assigned
                                    residue.id = (" ", chain_target_aligned[i-x].id[1]+1, " ")

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
    _used_structure_name = "6hup"
    _used_structure_loc = _structure_loc / f"{_used_structure_name}.pdb"
    _target_name = "a6b2g2_6hup"
    _target = [("Alpha-1", "Alpha-6"),
               ("Beta-3", "Beta-3"),
               ("Alpha-1", "Alpha-6"),
               ("Beta-3", "Beta-3"),
               ("Gamma-2", "Gamma-2")]

    _model_input = Input(_structure_loc,
                         _target,
                         _used_structure_name,
                         _target_name, _pir_input)

    _model_input.chain_preparation()
    # Cleaning and preparing the PDB of the template
    # Getting the start and end AA number in the cleaned template structure
    _model_input.pdb_clean()
    # Writing of the PIR/ALI file
    _model_input.ali_write()

    _model = Model(_model_input)
    _model.run_model()
    _model.chain_cleanup()
