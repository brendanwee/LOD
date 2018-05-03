import os, sys, argparse, subprocess

def parse_arguments():

    parser = argparse.ArgumentParser(description='Argument Parser')
    parser.add_argument('--refFa', dest='refFa', type=str, required=True)
    parser.add_argument('--refGTF', dest='refGTF', type=str, required=True)
    parser.add_argument('--RNASeq', dest='RNASeq', type=str, required=True)
    parser.add_argument('--RiboSeq', dest='RiboSeq', type=str, required=True)
    parser.add_argument('--f', dest='f', type=int, default=0.9)

    return parser.parse_args()
    

def exe_commands(args):

    # Create directory to store intermediate files

    os.popen("mkdir intermediate").read()
    # STAR
    print("Running HISAT2 ...")

    os.popen("STAR --runThreadN 7 --runMode genomeGenerate --genomeDir ./outputs/reference_genome/ --genomeFastaFiles "+args.refFa+"--sjdbGTFfile "+args.refGTF+"--sjdb 29").read()
    os.popen("STAR --runThreadN 7 --genomeDir ./outputs/reference_genome/ --readFilesIn "+args.RNASeq+"--outFileNamePrefix RNAseq").read()
    print("Finished!")

    # Samtools
    print("Running Samtools ...")
    os.popen("Samtools view -bS RNAseq.Aligned.out.sam | samtools sort > ./outputs/RNAseq.sorted.bam").read()
    print("Finished!")

    # StringTie
    print("Running StringTie ...")
    os.popen("stringtie ./outputs/RNAseq.sorted.bam -f "+args.f+"-p 16 -G "+args.refGTF+"-o ./outputs/RNAseq.transcriptome.gtf").read()
    print("Finished!")

    # STAR
    print("Running STAR ...")
    os.popen("STAR --runThreadN 7 --genomeDir ./outputs/experimental_genome/ --genomeFastaFiles "+args.refFa+"--sjdbGTFfile ./outputs/RNAseq.transcriptome.gtf").read()
    os.popen("STAR --outFilterType BySJout --runThreadN 7 --outFilterMismatchNmax 2 --genomeDir ./outputs/experimental_genome/ --readFilesIn "+args.RiboSeq+"--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFilterMultmapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd").read()
    os.popen("python preprocessGTF.py ./outputs/RNAseq.transcriptome.gtf "+args.refGTF+"RNAseq.processed.transcriptome.gtf").read()
    print("Finished!")

    # RiboCode
    print("Running RiboCode ...")
    os.popen("prepare_transcripts -g RNAseq.processed.transcriptome.gtf -o ./outputs/RiboCode_annot -f "+args.RiboSeq).read()
    os.popen("metaplots -a ./outputs/RiboCode_annot -r Riboseq.AlignedTo.Transcriptome.Bam").read()
    os.popen("RiboCode -a ./outputs/RiboCode_annot -c config.txt -l no -g -o RiboCode_ORFs_result").read()
    print("Finished!")

def main(args):

    #args = parse_arguments()
    exe_commands(args)
    
if __name__ == '__main__':
    main(sys.argv)