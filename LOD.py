#!/usr/bin/env python

import os, sys, argparse, subprocess, inspect

def parse_arguments():

    parser = argparse.ArgumentParser(description='Argument Parser')
    parser.add_argument('--refFa', dest='refFa', type=str, required=True)
    parser.add_argument('--refGTF', dest='refGTF', type=str, required=True)
    parser.add_argument('--RNASeq', dest='RNASeq', type=str, required=True)
    parser.add_argument('--RiboSeq', dest='RiboSeq', type=str, required=True)
    parser.add_argument('--f', dest='f', type=int, default=0.9)
    parser.add_argument('--n', dest='n', type=int, default=25)
    parser.add_argument('--i', dest='i', type=int, default=1)

    return parser.parse_args()
    
def get_bin():
    script_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    if os.path.isdir(script_dir+"/bin/"):
        return script_dir+"/bin/"
    else: # assume it is a symlink
        # ln -s /path/to/LOD.py /usr/local/bin/
        # any dir typ LOD tab -> LOD.py
        # call it
        # return path to project folder bin
        p = subprocess.Popen(["readlink `which LOD.py`"], stdout=subprocess.PIPE, shell=True)
        o, e = p.communicate()
        p.wait()
        #smthg/LOD.py
        folders = o.split("/")
        project_dir = "/".join(folders[:-1])
        return project_dir+"/bin/"



def exe_commands(args):

    # Create directory to store intermediate files
    bin = get_bin()

    os.popen("mkdir intermediate").read()
    # HISAT2


    if(args.i==1):
        # HISAT2
        print("Running HISAT2 ...")
        print "extracting splice sites"
        os.popen(bin + "hisat2_extract_splice_sites.py " + args.refGTF + " > ./intermediate/reference.ss").read()
        print "extracting exons"
        os.popen(bin + "hisat2_extract_exons.py " + args.refGTF + " > ./intermediate/reference.exon").read()
        print "building index"
        os.popen(bin + "hisat2-build -p " + str(args.n) + " --ss ./intermediate/reference.ss --exon ./intermediate/reference.exon " + args.refFa + " ./intermediate/reference_ht2_index").read()
        print "aligning reads"
        os.popen(bin + "hisat2 -p 15 -x ./intermediate/reference_ht2_index -U " + args.RNASeq + "-S ./intermediate/RNAseq.alignedto.reference.sam").read()
        print("Finished!")
    else:
        print("Index Provided: Running HISAT2 ...")
        os.popen(
            bin + "hisat2 -p 15 -x "+ args.i + " -U " + args.RNASeq + " -S ./intermediate/RNAseq.alignedto.reference.sam").read()
        print("Finished!")


    # Samtools
    print("Running Samtools ...")
    os.popen(bin + "samtools view -bS ./intermediate/RNAseq.alignedto.reference.sam > ./intermediate/RNAseq.alignedto.reference.bam")
    os.popen(bin + "samtools sort ./intermediate/RNAseq.alignedto.reference.bam > ./intermediate/RNAseq.sorted.bam").read()
    print("Finished!")

    # StringTie
    print("Running StringTie ...")
    os.popen(bin + "stringtie ./intermediate/RNAseq.sorted.bam -f " + str(args.f) + "-p 16 -G " + args.refGTF + "-o ./intermediate/RNAseq.transcriptome.gtf").read()
    print("Finished!")

    # preprocess
    print("Preprocessing ...")
    os.popen("python ./bin/preprocessGTF.py ./intermediate/RNAseq.transcriptome.gtf "+ args.refGTF + "./intermediate/RNAseq.processed.transcriptome.gtf").read()
    os.popen("python ./bin/ChangeName.py ./intermediate/RNAseq.processed.transcriptome.gtf ./intermediate/RNAseq.newnames.transcriptome.gtf").read()
    print("Finished!")
    
    # STAR
    print("Running STAR ...")
    os.popen(bin + "STAR --runThreadN " + str(args.n) + "--genomeDir ./intermediate/experimental_genome/ --genomeFastaFiles " + args.refFa + " --sjdbGTFfile ./intermediate/RNAseq.newnames.transcriptome.gtf").read()
    os.popen(bin + "STAR --outFilterType BySJout --runThreadN " + str(args.n) + " --outFilterMismatchNmax 2 --genomeDir ./intermediate/experimental_genome/ --readFilesIn "+ args.RiboSeq +" --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFilterMultmapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd --outFileNamePrefix ./intermediate/Riboseq").read()
    
    print("Finished!")
    
    # RiboCode
    print("Running RiboCode ...")
    os.popen("GTFupdate ./intermediate/RNAseq.newnames.transcriptome.gtf > ./intermediate/RNAseq.processed.updated.gtf").read()
    os.popen("prepare_transcripts -g ./intermediate/RNAseq.processed.updated.gtf -o ./intermediate/RiboCode_annot -f " + args.refFa).read()
    os.popen("metaplots -a ./intermediate/RiboCode_annot -r ./intermediate/Riboseq.AlignedTo.Transcriptome.Bam").read()
    os.popen("RiboCode -a ./intermediate/RiboCode_annot -c metaplots_pre_config.txt -l no -g -o RiboCode_ORFs_result").read()
    print("Finished!")
    
def main(args):

    args = parse_arguments()
    exe_commands(args)
    
if __name__ == '__main__':
    main(sys.argv)
