sed '1,1d' edgeR_saida/salmon.gene.counts.matrix.sp_log_vs_sp_plat.edgeR.DE_results| awk '{ if ($7 <= 0.05) print $1;}' > DE_genes.txt

cat trinity_saida/Trinity.fasta | tr -d '\n' > temp1.txt

sed 's/>/\n>/g' temp1.txt > temp2.txt

sed 's/]A/]\nA/g' temp2.txt > temp3.txt

sed 's/]T/]\nT/g' temp3.txt > temp4.txt

sed 's/]C/]\nC/g' temp4.txt > temp5.txt

sed 's/]G/]\nG/g' temp5.txt > temp6.txt

grep -A 1 -f DE_genes.txt temp6.txt --no-group-separator > DE_genes.fasta

rm temp1.txt temp2.txt temp3.txt temp4.txt temp5.txt temp6.txt
