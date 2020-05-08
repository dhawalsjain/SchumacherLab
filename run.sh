### RNA fold -23
PATH=$PATH:/n/data1/hms/dbmi/park/Dhawal/Softwares/ViennaRNA-2.4.9/SW/bin
sbatch -t 0-11:59 -n 4 --mem-per-cpu=4G -p short --wrap "RNAfold -p -d2 --noLP DRE_Clustered_RefSeqs_Genomic.fa > DRE_Clustered_RefSeqs_Genomic_OUT.fa"
sbatch -t 0-11:59 -n 4 --mem-per-cpu=4G -p short --wrap " RNAfold -p -d2 --noLP DRE_Clustered_AltSeqs_Genomic.fa > DRE_Clustered_AltSeqs_Genomic_OUT.fa"

module load ghostscript/9.24 
for f in *_ss.ps
do
f0=$(echo "$f" | sed -e 's/_ss.ps//g')
relplot.pl ${f} ${f0}_dp.ps > ${f0}_rel.ps
gs -sDEVICE=jpeg -r300 -sPAPERSIZE=a4 -dBATCH -dNOPAUSE -sOutputFile=${f0}_rel.jpg ${f0}_rel.ps
done


### RNA fold -33
PATH=$PATH:/n/data1/hms/dbmi/park/Dhawal/Softwares/ViennaRNA-2.4.9/SW/bin
sbatch -t 0-11:59 -n 4 --mem-per-cpu=4G -p short --wrap "RNAfold -p -d2 --noLP 33_DRE_Clustered_RefSeqs_Genomic.fa > 33_DRE_Clustered_RefSeqs_Genomic_OUT.fa"
sbatch -t 0-11:59 -n 4 --mem-per-cpu=4G -p short --wrap " RNAfold -p -d2 --noLP 33_DRE_Clustered_AltSeqs_Genomic.fa > 33_DRE_Clustered_AltSeqs_Genomic_OUT.fa"

module load ghostscript/9.24 
for f in *_ss.ps
do
f0=$(echo "$f" | sed -e 's/_ss.ps//g')
relplot.pl ${f} ${f0}_dp.ps > ${f0}_rel.ps
gs -sDEVICE=jpeg -r300 -sPAPERSIZE=a4 -dBATCH -dNOPAUSE -sOutputFile=${f0}_rel.jpg ${f0}_rel.ps
done


### RNA fold -33hc
PATH=$PATH:/n/data1/hms/dbmi/park/Dhawal/Softwares/ViennaRNA-2.4.9/SW/bin
sbatch -t 0-11:59 -n 4 --mem-per-cpu=4G -p short --wrap "RNAfold -p -d2 --noLP 33hc_DRE_Clustered_RefSeqs_Genomic.fa > 33hc_DRE_Clustered_RefSeqs_Genomic_OUT.fa"
sbatch -t 0-11:59 -n 4 --mem-per-cpu=4G -p short --wrap " RNAfold -p -d2 --noLP 33hc_DRE_Clustered_AltSeqs_Genomic.fa > 33hc_DRE_Clustered_AltSeqs_Genomic_OUT.fa"

module load ghostscript/9.24 
for f in *_ss.ps
do
f0=$(echo "$f" | sed -e 's/_ss.ps//g')
relplot.pl ${f} ${f0}_dp.ps > ${f0}_rel.ps
gs -sDEVICE=jpeg -r300 -sPAPERSIZE=a4 -dBATCH -dNOPAUSE -sOutputFile=${f0}_rel.jpg ${f0}_rel.ps
done



for f in *.fastq
do
sbatch -t 0-11:59 -n 2 --mem-per-cpu=2G -p short --wrap "gzip -c $f > ${f}.gz"
done




## meme
PATH=$PATH:/home/dj70/meme/bin
meme Edit_25bpContext.fa -oc all -dna -maxsize 100000000 -revcomp 
meme Edit_25bpContext_Control.fa -oc control -dna -maxsize 100000000 -revcomp -nmotifs 3
meme Edit_25bpContext_D7.fa -oc D7 -dna -maxsize 100000000 -revcomp -nmotifs 3
meme Edit_25bpContext_D14.fa -oc D14 -dna -maxsize 100000000 -revcomp -nmotifs 3

fimo --thresh 1e-2 --oc control_fimo edits.meme Edit_25bpContext_Control.fa
fimo --thresh 1e-2 --oc D7_fimo edits.meme Edit_25bpContext_D7.fa
fimo --thresh 1e-2 --oc D14_fimo edits.meme Edit_25bpContext_D14.fa


## meme _33
PATH=$PATH:/home/dj70/meme/bin
sbatch -t 0-11:59 -n 2 --mem-per-cpu=2G -p short --wrap "meme 33_Edit_25bpContext_all.fa -oc all -dna -maxsize 100000000 -revcomp" 
sbatch -t 0-11:59 -n 2 --mem-per-cpu=2G -p short --wrap "meme 33_Edit_25bpContext_Control.fa -oc control -dna -maxsize 100000000 -revcomp -nmotifs 3"
sbatch -t 0-11:59 -n 2 --mem-per-cpu=2G -p short --wrap "meme 33_Edit_25bpContext_D7.fa -oc D7 -dna -maxsize 100000000 -revcomp -nmotifs 3"
sbatch -t 0-11:59 -n 2 --mem-per-cpu=2G -p short --wrap "meme 33_Edit_25bpContext_D14.fa -oc D14 -dna -maxsize 100000000 -revcomp -nmotifs 3"

fimo --thresh 1e-2 --oc control_fimo 33_edits.meme 33_Edit_25bpContext_Control.fa
fimo --thresh 1e-2 --oc D7_fimo 33_edits.meme 33_Edit_25bpContext_D7.fa
fimo --thresh 1e-2 --oc D14_fimo 33_edits.meme 33_Edit_25bpContext_D14.fa






PATH=$PATH:/home/dj70/meme/bin
sbatch -t 0-11:59 -n 2 --mem-per-cpu=2G -p short --wrap "meme 33_Edit_25bpContext_3_prime_UTR.fa -oc utr3 -dna -maxsize 100000000 -revcomp; \ 
meme 33_Edit_25bpContext_intergenic_region.fa -oc intergenic -dna -maxsize 100000000 -revcomp; \ 
meme 33_Edit_25bpContext_intron.fa -oc intron -dna -maxsize 100000000 -revcomp; \ 
meme 33_Edit_25bpContext_missense.fa -oc missense -dna -maxsize 100000000 -revcomp; \ 
meme 33_Edit_25bpContext_splice.fa -oc splice -dna -maxsize 100000000 -revcomp; \ 
meme 33_Edit_25bpContext_synonymous.fa -oc synonym -dna -maxsize 100000000 -revcomp; \
" 




## meme
PATH=$PATH:/home/dj70/meme/bin
sbatch -t 0-11:59 -n 2 --mem-per-cpu=2G -p short --wrap "meme 33hc_Edit_25bpContext_3_prime_UTR.fa -oc utr3 -dna -maxsize 100000000 -revcomp; \ 
meme 33hc_Edit_25bpContext_all.fa -oc all -dna -maxsize 100000000 -revcomp; \ 
meme 33hc_Edit_25bpContext_intergenic_region.fa -oc intergenic -dna -maxsize 100000000 -revcomp; \ 
meme 33hc_Edit_25bpContext_intron.fa -oc intron -dna -maxsize 100000000 -revcomp; \ 
meme 33hc_Edit_25bpContext_missense.fa -oc missense -dna -maxsize 100000000 -revcomp; \ 
meme 33hc_Edit_25bpContext_splice.fa -oc splice -dna -maxsize 100000000 -revcomp; \ 
meme 33hc_Edit_25bpContext_synonymous.fa -oc synonym -dna -maxsize 100000000 -revcomp; \
" 









#################################################################################################3
########### hyperedited reads
PATH=$PATH:/n/data1/hms/dbmi/park/Dhawal/Softwares/mpileup2readcounts
PATH=$PATH:/n/data1/hms/dbmi/park/Dhawal/Softwares/bwa_0.7.17/

sbatch -t 0-11:59 -n 4 --mem-per-cpu=8G -p short --wrap "\
sed -i 's/A/G/g' mm10_A2G.fa; \
bwa index mm10_A2G.fa; \
"

for f in /n/data1/hms/dbmi/park/DATA/Kreidberg/Sept2019_VSchumacker/trim_fastq/rnaEditor/*_R1_trim/*_R1_trim.bam
do
f0=$(echo "$f" | sed -e 's/_R1_trim.bam//g')
f0=$(echo "$f0" | sed -e 's/.*\///g')
echo $f0
sbatch -t 0-11:59 -n 4 --mem-per-cpu=8G -p short --wrap "\
perl extractUnmap.pl -bam $f; \
bwa mem -t 8 /n/data1/hms/dbmi/park/Dhawal/Genomes/mm10/bwa_A2G/mm10_A2G.fa ${f0}_R1_trim.fq.gz | perl EditReport.pl -bam - -hd header.sam -out ${f0};
samtools sort -@ 8 -o ${f0}.sort.vis.bam ${f0}.vis.bam; \
samtools index ${f0}.sort.vis.bam; \ 
"  
done


for f in /n/data1/hms/dbmi/park/DATA/Kreidberg/Sept2019_VSchumacker/trim_fastq/rnaEditor/*_R1_trim/*_R1_trim.bam
do
f0=$(echo "$f" | sed -e 's/_R1_trim.bam//g')
f0=$(echo "$f0" | sed -e 's/.*\///g')
echo $f0
sbatch -t 0-31:59 -n 4 --mem-per-cpu=8G -p medium --wrap "\
bwa mem -t 8 /n/data1/hms/dbmi/park/Dhawal/Genomes/mm10/bwa_A2G/mm10_A2G.fa ${f0}_R1_trim.fq.gz | perl EditReport.pl -bam - -hd header.sam -out ${f0}_AG -type AG;
samtools sort -@ 8 -o ${f0}_AG.sort.vis.bam ${f0}_AG.vis.bam; \
samtools index ${f0}_AG.sort.vis.bam; \ 
"  
done

for f in /n/data1/hms/dbmi/park/DATA/Kreidberg/Sept2019_VSchumacker/trim_fastq/rnaEditor/*_R1_trim/*_R1_trim.bam
do
f0=$(echo "$f" | sed -e 's/_R1_trim.bam//g')
f0=$(echo "$f0" | sed -e 's/.*\///g')
echo $f0
sbatch -t 0-31:59 -n 4 --mem-per-cpu=8G -p medium --wrap "\
bwa mem -t 8 /n/data1/hms/dbmi/park/Dhawal/Genomes/mm10/bwa_A2G/mm10_A2G.fa ${f0}_R1_trim.fq.gz | perl EditReport.pl -bam - -hd header.sam -out ${f0}_TC -type TC;
samtools sort -@ 8 -o ${f0}_TC.sort.vis.bam ${f0}_TC.vis.bam; \
samtools index ${f0}_TC.sort.vis.bam; \ 
"  
done

mkdir tmp
for f in *.report.txt.gz
do
f0=$(echo "$f" | sed -e 's/.report.txt.gz/.report.uq.txt.gz/g')
echo $f0
sbatch -t 0-11:59 -n 4 --mem-per-cpu=2G -p short --wrap "zcat $f | sort -k1,1 -k2,2n -k5,5 --parallel=4 -T tmp -u -S 1G - |gzip -c - > ${f0}"
done    
