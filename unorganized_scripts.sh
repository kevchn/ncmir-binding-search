# Download 3' and 5' (control) UTR regions from UTRdb
3UTRaspic.Hum.fasta
5UTRapsic.Hum.fasta

# Remove line breaks from fasta
awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' 3UTRaspic.Hum.fasta > 3UTR.fasta
awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' 5UTRaspic.Hum.fasta > 5UTR.fasta

# Find sequences with 7nt motif
grep -i -B 1 'GCGGAAC' 3UTR.fasta > 3UTR.GCGGAAC.fasta
grep -i -B 1 'GCGGAAC' 5UTR.fasta > 5UTR.GCGGAAC.fasta

# Find sequences with 7nt motif and AA
grep -i -B 1 'AAGCGGAAC' 3UTR.fasta > 3UTR.AAGCGGAAC.fasta
grep -i -B 1 'AAGCGGAAC' 5UTR.fasta > 5UTR.AAGCGGAAC.fasta

grep -i -B 1 'CCGCGGAAC' 3UTR.fasta > 3UTR.CCGCGGAAC.fasta
grep -i -B 1 'CCGCGGAAC' 5UTR.fasta > 5UTR.CCGCGGAAC.fasta

### Analyze file sizes / # of sequences
ls -lah
wc -l

232K 3UTR.AAGCGGAAC.fasta (0.0014) (0.0473 e=0.0625)
4.9M 3UTR.GCGGAAC.fasta   (0.0295)
166M 3UTR.fasta

51K 5UTR.AAGCGGAAC.fasta  (0.0014) (0.0561 e=0.0625)
908K 5UTR.GCGGAAC.fasta   (0.0245)
37M 5UTR.fasta

102 3UTR.AAGCGGAAC.fasta  (0.00026)
2393 3UTR.GCGGAAC.fasta   (0.00615)
389006 3UTR.fasta

124 5UTR.AAGCGGAAC.fasta  (0.00049)
2261 5UTR.GCGGAAC.fasta   (0.00909)
248630 5UTR.fasta

No significantly greater amount of 3p UTR regions that can bind to miRNA motif than 5p UTR regions. (Actually slightly less)
No significantly greater amount of AA in 3p UTR regions that can bind to miRNA motif than 5p UTR regions. (Actually slightly less)
No significantly greater amount of AA in any UTR regions than any other nt combinations. (Actually slightly less)

All 3UTR AA sequences have only 1 binding site
ALL 3UTR no AA sequences have only 1 binding site

### Analyze AA content
grep -o ..gcggaac 3UTR.GCGGAAC.fasta > 3UTR.GCGAAC.cut.fasta
grep -o ..gcggaac 5UTR.GCGGAAC.fasta > 5UTR.GCGAAC.cut.fasta
sed -i -e 's/gcggaac//g' 3UTR.GCGAAC.cut.fasta 
sed -i -e 's/gcggaac//g' 5UTR.GCGAAC.cut.fasta 
sort 3UTR.GCGAAC.cut.fasta | uniq -c > 3UTR.GCGAAC.nt_distribution
sort 5UTR.GCGAAC.cut.fasta | uniq -c > 5UTR.GCGAAC.nt_distribution

nt_distribution_plot.png 
ca ct ta are enriched, cc (validated, no) and aa (validated, yes) are not

### Find 3' UTR gene names
grep '>' 3UTR.AAGCGGAAC.fasta > 3UTR.AAGCGGAAC.fasta.names
grep '>' 3UTR.CCGCGGAAC.fasta > 3UTR.CCGCGGAAC.fasta.names

grep 'AC\|DE' 5UTRaspic.Hum.dat > 5UTR.genenames

sed 's/;//g' 3UTR.genenames > temp
cat temp > 3UTR.genenames
rm temp

sed -e 's/\(>\).*\( \)/\1\2/' 3UTR.genenames | cut -b 5- > temp
cat temp > 3UTR.genenames
rm temp

while read line; do echo "$line|"; done < 7nt/3UTR.GCGGAAC.fasta.names | tr -d '\n' > temp
cat temp > 7nt/3UTR.GCGGAAC.fasta.names
rm temp

RM ALL NEW LINES
sed -i ':a;N;$!ba;s/\n/\t/g' 3UTR.genenames
gsed 's/CA......./\n&/g' < 3UTR.genenames > temp
cat temp > 3UTR.genenames
rm temp
#
python insert_genenames.py -f isoform/5UTR.AAGCGGAAC.fasta 5UTR.AAGCGGAAC.genenames

# Move (string) to beginning of line
gsed -r 's/(.*)\((.*)\)/>\2\t\1/' 3UTR.fasta_with_names
