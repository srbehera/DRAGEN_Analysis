import argparse
import re

header = """##FILTER=<ID=HOMREF,Description="Homozygous reference call">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=EVENTTYPE,Number=A,Type=String,Description="Type of associated event">
##EVENTTYPE=<ID=STR,Description="Short tandem repeat">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the longest variant described in this record">
##INFO=<ID=REFRUC,Number=1,Type=Float,Description="Reference copy number of the tandem repeat">
##INFO=<ID=CN,Number=A,Type=Float,Description="Copy number of allele">
##INFO=<ID=RUS,Number=1,Type=String,Description="Repeat unit sequence of the corresponding repeat sequence">
##INFO=<ID=RUL,Number=1,Type=Integer,Description="Repeat unit length of the corresponding repeat sequence">
##INFO=<ID=RUC,Number=A,Type=Float,Description="Repeat unit count of the corresponding repeat sequence">
##INFO=<ID=RUCCHANGE,Number=A,Type=Float,Description="Change in repeat unit count of the corresponding repeat sequence">
##INFO=<ID=CNVTRLEN,Number=A,Type=Float,Description="Change in total length of the corresponding repeat sequence">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=CN,Number=1,Type=Float,Description="Sum of copy numbers">
##ALT=<ID=CNV:TR,Description="Tandem repeat determined based on DNA abundance">
"""


def translateLine(line: str) -> str:
    lineSplit = line.strip("\n").split("\t")

    contigsAndAlleles = lineSplit[:-3]
    contigsAndAlleles = "\t".join(contigsAndAlleles)
    contigsAndAlleles = re.sub(r"<STR\d+>", "<CNV:TR>", contigsAndAlleles)

    infoField = lineSplit[-3]
    infoDict = re.findall(r"([A-Z]+)=([\+\-A-Za-z\.\d_:]+);", infoField)
    assert len(infoDict) != 0
    infoDict = dict(infoDict)

    formatKeys = lineSplit[-2].split(":")
    formatValues = lineSplit[-1].split(":")
    formatDict = dict(zip(formatKeys, formatValues))

    contigsAndAlleles = re.sub(r"\.", infoDict["VARID"], contigsAndAlleles, count=1)

    rul = len(infoDict["RU"])
    refruc = float(infoDict["RL"]) / rul

    if formatDict["GT"] in ["./.", "."]:
        eventtype = "."
        svlen = infoDict["RL"]
        cn = "."
        ruc = "."
        rucchange = "."
        cnvtrlen = "."
        cnformat = "."
    else:
        if formatDict["GT"] == "1/2":
            eventtype = ["STR", "STR"]
            svlen = [infoDict["RL"]] * 2
            ruc = [float(x) for x in formatDict["REPCN"].split("/")]
            cn = [x / refruc for x in ruc]
            rucchange = [x - refruc for x in ruc]
            cnvtrlen = [x * rul for x in rucchange]
        else:
            eventtype = ["STR"]
            svlen = [infoDict["RL"]]
            ruc = float(formatDict["REPCN"].split("/")[-1])
            cn = [ruc / refruc]
            rucchange = [ruc - refruc]
            cnvtrlen = [rucchange[0] * rul]
            ruc = [ruc]
        eventtype = ",".join(eventtype)
        cnformat = f"{sum(cn):.6f}"
        svlen = ",".join(svlen)
        cn = ",".join([f"{x:.2f}" for x in cn])
        ruc = ",".join([f"{x:.2f}" for x in ruc])
        rucchange = ",".join([f"{x:.6f}" for x in rucchange])
        cnvtrlen = ",".join([f"{x:.6f}" for x in cnvtrlen])

    info = ["SVTYPE=CNV"]
    info.append(f"EVENTTYPE={eventtype}")
    info.append(f"SVLEN={svlen}")
    info.append(f"END={infoDict['END']}")
    info.append(f"REFRUC={refruc:.2f}")
    info.append(f"RUS={infoDict['RU']}")
    # info.append(f"RUL={len(infoDict['RU'])}") Gives issues with Wittyer (RulMismatch) not sure why
    info.append(f"CN={cn}")
    info.append(f"RUC={ruc}")
    info.append(f"RUCCHANGE={rucchange}")
    info.append(f"CNVTRLEN={cnvtrlen}")

    info = ";".join(info)

    formatFields = ":".join(["GT", "CN"])

    gt = formatDict["GT"]

    format = ":".join([gt, cnformat])

    """
    HOMREF calls have to be "masked" when using Witty.er, otherwise a HOMREF call that is not
    present in the truth will be counted as a query FP, decreasing precision for no reason.
    It's sufficient to change the FILTER from PASS or . to something else (like HOMREF)
    """
    if gt in ["0/0", "0"]:
        contigsAndAlleles = re.sub("PASS", "HomRef", contigsAndAlleles)

    return "\t".join([contigsAndAlleles, info, formatFields, format]) + "\n"


def convertVcf(inputFile: str, outputFile: str) -> None:
    inputHandle = open(inputFile, "r")
    outputHandle = open(outputFile, "w")

    for line in inputHandle:
        if re.search("^#[^#]", line):
            outputHandle.write(header)
            outputHandle.write(line)
            continue

        if re.search("^##(INFO)|(FORMAT)|(ALT)", line):
            continue

        if re.search("^##", line):
            outputHandle.write(line)
            continue

        towrite = translateLine(line)
        outputHandle.write(towrite)

    inputHandle.close()
    outputHandle.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert the repeat genotyping module VCF output in a VCFv4.4 compliant format (similar to the VNTR module output)")
    parser.add_argument("-i", "--input", help="input VCF", type=str, required=True)
    parser.add_argument("-o", "--output", help="output VCF", type=str, required=True)
    args = parser.parse_args()

    convertVcf(args.input, args.output)
