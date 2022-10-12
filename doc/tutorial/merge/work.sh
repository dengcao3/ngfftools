../../../bin/ngfftools merge --list merge_gff3.list --outfile merge_gff3.ngff --overlap_percent 0.8 &>merge_gff3.error
../../../bin/ngfftools merge --list merge_gtf.list --outfile merge_gtf.gtf --overlap_percent 0.8 --header forMerge.gtf &>merge_gtf.error
../../../bin/ngfftools merge --list mergeTwo.list --outfile mergeTwo.ngff --overlap_percent 0.8 &>mergeTwo.error
