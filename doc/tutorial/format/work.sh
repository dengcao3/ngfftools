../../../bin/ngfftools format --infile template.ngff --informat ngff --outfile out.gff3 &>ngff_to_gff3.error
cat template.ngff | ../../../bin/ngfftools format --infile - --informat ngff --outfile OUT.gtf &>ngff_to_gtf.error
../../../bin/ngfftools format --infile template.gtf --informat gtf --outfile gtf2.gff3 &>gtf_to_gff3.error 
../../../bin/ngfftools format --infile template.gtf --informat gtf --outfile gtf2.ngff &>gtf_to_ngff.error
../../../bin/ngfftools format --infile template.gff3 --informat gff3 --outformat ngff --prefix template_ >GFF3.ngff 2>gff3_to_ngff.error
../../../bin/ngfftools format --infile sub.GCF_000001405.38_GRCh38.p12_genomic.gff --informat gff3 --outformat ngff >sub.GCF_000001405.38_GRCh38.p12_genomic.ngff 2>sub.GCF_000001405.38_GRCh38.p12_genomic.gff3_to_ngff.error
../../../bin/ngfftools format --name_type ncbi --infile sub.GCF_000001405.38_GRCh38.p12_genomic.gff --informat gff3 --outformat ngff >sub.GCF_000001405.38_GRCh38.p12_genomic.protein_id.ngff 2>sub.GCF_000001405.38_GRCh38.p12_genomic.gff3_to_ngff.protein_id.error
../../../bin/ngfftools format --infile fusion.ngff --informat ngff --outformat ngff >fusion.new.ngff 2>fusion.error
