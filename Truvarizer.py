import argparse
import re
import math


def addSequenceCalls(line: str) -> str:
    lineSplit = line.strip("\n").split("\t")

    contigsAndAlleles = lineSplit[:-3]
    contigsAndAlleles = "\t".join(contigsAndAlleles)

    infoField = lineSplit[-3]
    infoDict = re.findall("([A-Z]+)=([\+\-A-Za-z\.\d_:]+)", infoField)
    assert len(infoDict) != 0
    infoDict = dict(infoDict)

    if infoDict["CNVTRLEN"] == "." or infoDict["CNVTRLEN"] == "0":
        return ""

    motif = infoDict["RUS"]
    svlen = int(infoDict["SVLEN"])
    altlen = int(infoDict["CNVTRLEN"])
    repeat = motif * abs(max(svlen, svlen + altlen))
    ref = repeat[:svlen]
    alt = repeat[: (svlen + altlen)]
    svtype = "REF"
    end = infoDict["END"]
    lineSplit[4] = lineSplit[3] + alt
    lineSplit[3] = lineSplit[3] + ref

    if altlen == 0:
        lineSplit[3] = "."
        lineSplit[4] = "."
    elif altlen < 0:
        svtype = "DEL"
    else:
        svtype = "INS"

    return (
        "\t".join(
            lineSplit[0:7]
            + [
                f"SVTYPE={svtype};SVLEN={altlen};RUS={motif}"
            ]
            + lineSplit[-2:]
        )
        + "\n"
    )


def convertVcf(inputFile: str, outputFile: str) -> None:
    inputHandle = open(inputFile, "r")
    outputHandle = open(outputFile, "w")

    for line in inputHandle:
        if re.search("^#", line):
            outputHandle.write(line)
            continue

        towrite = addSequenceCalls(line)
        if towrite != "":
            outputHandle.write(towrite)

    inputHandle.close()
    outputHandle.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert the repeat genotyping module VCF output in a VCFv4.4 compliant format (similar to the VNTR module output)"
    )
    parser.add_argument(
        "-i", "--input", help="input VCF", type=str, required=True
    )
    parser.add_argument(
        "-o", "--output", help="output VCF", type=str, required=True
    )
    args = parser.parse_args()

    convertVcf(args.input, args.output)
