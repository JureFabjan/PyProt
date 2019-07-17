import modeller


if __name__ == "__main__":
    # Variables for the input

    # PIR file (.ali) contains the sequences of the target protein in the format:
    # >P1;name
    # sequence:model_name
    # AA SEQUENCE*
    pir_input = ""

    modeller.log.verbose()
    environment = modeller.environ()

    # Database of sequence alignments reading/writing
    # Modeller can create a binary DB, which works faster  than reading PIR file.
    # We should opt for this once we create the alignment master-PIR.

    # Alignment object creation
    alignment = modeller.alignment(environment)
    