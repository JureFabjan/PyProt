import modeller
import Bio.AlignIO as AlignIO


def seq_clean(seq1, seq2):
    """
    Cleans the gaps from the sequences.
    :param seq1: First SeqRecord object
    :param seq2: Second SeqRecord object
    :return: (SeqRecord, SeqRecord, array) representing two cleaned sequences and an array of coordinates
    of removed spots.
    """
    # Get the coordinates where there are gaps.
    indices = sorted(list(set(pos for pos, aa in enumerate(seq1) if aa == "-") | set(pos for pos, aa in enumerate(seq2) if aa == "-")), reverse=True)
    for index in indices:
        seq1 = seq1[:index] + seq1[index+1:]
        seq2 = seq2[:index] + seq2[index+1:]
    return seq1, seq2, indices[::-1]


if __name__ == "__main__":
    # Variables for the input

    # PIR file (.ali) contains the sequences of the target protein in the format:
    # >P1;name
    # sequence:model_name
    # AA SEQUENCE*
    pir_input = "C:/Users/Jure/Documents/Jure files/GitHub/ligand_dock_project/MasterAli.pir"
    sequences = AlignIO.read(pir_input, "pir")

    # modeller.log.verbose()
    # environment = modeller.environ()

    # Database of sequence alignments reading/writing
    # Modeller can create a binary DB, which works faster  than reading PIR file.
    # We should opt for this once we create the alignment master-PIR.

    # Alignment object creation
    # alignment = modeller.alignment(environment)
