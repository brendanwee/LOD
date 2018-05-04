

# this script adds additional attribute information to a GTF file from a reference GTF file.
# Start and stop codon information and attributes like biotype if available


from sys import argv


if len(argv) != 4:
    print """USAGE: python preprocessGTF.py <experimental_GTF> <reference_GTF> <outfile_name>
    This script adds additional attribute information to a GTF file from a reference GTF file.
    Start and stop codon information and attributes like biotype if available"""
    exit(0)


s, exp, ref, out = argv

""" DEBUGGING ARTIFACT. PLZ IGNORE or utilize if you find this has a bug when you use it.
exp = "/Users/Brendan/cmu/dataPracticum/project/newestAnnotations/RNAseq.transcriptome.gtf"
ref  = "/Users/Brendan/cmu/dataPracticum/project/newestAnnotations/ref.gtf"
out = "/Users/Brendan/cmu/dataPracticum/project/newestAnnotations/out.gtf"
"""



# Makes a dictionary mapping Transcript ID's of exons to their index in data
def MakeReferenceDictionary(data):
    ref_dict = {}
    for i in range(0,len(data)):
        columns = data[i].split("\t")
        if columns[0][0] == "#":
            continue
        if columns[2] == "exon":
            attr = columns[8].split(";")
            transcript_id = FindTranscriptID(attr)
            if transcript_id not in ref_dict:
                ref_dict[transcript_id] = i
    return ref_dict

def FindTranscriptID(attr):
    id = None
    for a in attr:
        if "transcript_id" in a:
            id = a.split(" ")[2]
            break
    return id



# iterates through the reference data until the entry matching the ref_id is found
# Assuming the entries are for a given transcript are grouped together in hierachical order,
# the exon line, start_codon line, and stop_codon line are returned if found.
# if not found, then None is returned
# this is the computational bottleneck - The search is linear and iterates through the entire reference file.

def FindReferenceEntry(ref_id, exon_num, data, ref_dict, exp_attr):
    ref_entry = None
    startline = None
    stopline = None
    CDS = None
    try:
        i = ref_dict[ref_id]
    except KeyError:
        print ref_id + " not found in the reference file. Please ensure you are using the correct reference file and that" \
                       "your experimental GTF was created with this reference as a guide"
        exit(1)
    entry = data[i]
    columns = entry.split("\t")
    attr = columns[8].split(";")
    this_id = FindTranscriptID(attr)
    while ref_id == this_id:
        if columns[2] == "exon":
            if exon_num:
                for item in attr:
                    if "exon_number" in item and exon_num == item.split(" ")[2]:
                            ref_entry = entry
                            break
            else:
                ref_entry = entry
        if columns[2] == "CDS":
            if exon_num:
                for item in attr:
                    if "exon_number" in item and exon_num == item.split(" ")[2]:
                        CDS = "\t".join(columns[:8]) + "\t" + ";".join(exp_attr) + "\n" # edit attributes to match gene_id and transcript_id of experimental gtf - exp_attr
        if columns[2] == "start_codon":
            if exon_num:
                for item in attr:
                    if "exon_number" in item and exon_num == item.split(" ")[2]:
                        startline = "\t".join(columns[:8]) + "\t" + ";".join(exp_attr) + "\n" # edit attributes to match gene_id and transcript_id of experimental gtf - exp_attr
        if columns[2] == "stop_codon":
            if exon_num:
                for item in attr:
                    if "exon_number" in item and exon_num == item.split(" ")[2]:
                        stopline = "\t".join(columns[:8]) + "\t" + ";".join(exp_attr) + "\n"
        i += 1
        entry = data[i]
        columns = entry.split("\t")
        attr = columns[8].split(";")
        this_id = FindTranscriptID(attr)
    return ref_entry, startline, stopline, CDS


# given the attribute column of the GTF as a string, finds the exon number if available
def FindExonNum(attr):
    for item in attr:
        if "exon_number" in item:
            return item.split(" ")[2]
    return False


# Extract additional attribute data from the reference entry matching your exon. Attempts to avoid redundant data
def ExtractAttrItems(ref_entry):
    attr = ref_entry.split("\t")
    items = attr[8].split(";")
    i = 0
    while True:
        if "gene_name" in items[i]:
            break
        i += 1

    return ";".join(items[i:])


# opens the three files provided as args
# writes the data from the experimental GTF into the outfile
# upon finding an entry with reference data matching, the extra data is extracted from the reference GTF
# this data is written to the outfile maintaining Hierachical structure
def main():
    outfile = open(out, 'a+')
    outfile.seek(0,0)
    test = outfile.read().strip()
    if test != "":
        print "Outfile already exists"
        outfile.close()
        exit(1)
    infile = open(exp,'r')
    infile.seek(0,0)
    datafile = open(ref,"r")
    datafile.seek(0,0)
    data = datafile.readlines()
    datafile.close()
    ref_dict = MakeReferenceDictionary(data)
    line = infile.readline()
    while line[0] == "#":
        outfile.write(line)
        line = infile.readline()
    while line != "":
        line_written = False
        l = line.strip().split("\t")
        if l[6] == ".": # ribocode requires that every read have a strand sense. The bam files are each mapped to transcipts here, so we cant remove any transcripts. + is arbitrarily chosen
            l[6] = "+"
            line = "\t".join(l) + "\n"
            outfile.write(line)
            line_written = True
            pass
        elif l[2] == "exon":
            attr = l[8].split(";")
            for item in attr:
                if "reference_id" in item:
                    exon_num = FindExonNum(attr)
                    #get reference id from attr
                    ref_id = item.split(" ")[2]
                    ref_entry, startline, stopline, CDS = FindReferenceEntry(ref_id, exon_num, data, ref_dict, attr)
                    ref_attr_add_ons = ExtractAttrItems(ref_entry) # joined string to add at end of line
                    line = line[:-1] + ref_attr_add_ons
                    outfile.write(line)
                    line_written = True
                    if CDS:
                        outfile.write(CDS)
                    if startline:
                        outfile.write(startline)
                    if stopline:
                        outfile.write(stopline)
                    break
        if not line_written:
            outfile.write(line)
        line = infile.readline()

    outfile.close()
    infile.close()

if __name__ == "__main__":
    main()