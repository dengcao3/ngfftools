cat template.ngff |../../../bin/ngfftools abs2rel --abs2rel 'product=relative' --infile - --informat ngff --outfile product_relative.ngff &>abs2rel.productRel.error
../../../bin/ngfftools abs2rel --abs2rel 'processed=relative;product=relative' --infile template.ngff --outfile proc_prod_relative.ngff &>abs2rel.productRel_processRel.error
../../../bin/ngfftools abs2rel --abs2rel 'product=relative' --infile fusion.ngff --outfile fusion.proc_prod_relative.ngff &>abs2rel.productRel_processRel.fusion.error
../../../bin/ngfftools abs2rel --abs2rel 'product=relative' --infile template.2.ngff --informat ngff --outfile product_relative.2.ngff
