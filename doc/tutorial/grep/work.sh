../../../bin/ngfftools grep --infile template.ngff --outfile grep1.ngff --grep 'seqId !~ /ChrM/'
../../../bin/ngfftools grep --infile template.ngff --outfile grep2.ngff --grep 'level == "locus"'
../../../bin/ngfftools grep --infile template.ngff --outfile grep3.ngff --grep 'level == "processed" && product_biotype == "ORF"'
../../../bin/ngfftools grep --infile template.ngff --outfile grep4.ngff --grep 'level == "processed" && parentId != "L01.pt1"'
