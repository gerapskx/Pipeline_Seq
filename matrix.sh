#RemoveGeneID_fromfiles
ls -1 *.txt | parallel 'cat {} | sed 1d | cut -f2 {} > {/.}_clean.txt'

#ExtractGeneID
ls -1 *.txt | head -1 | xargs cut -f1 > genes.txt

#PasteGeneIDwithmatchingcount
paste genes.txt *clean.txt > Dros_matrix.txt
