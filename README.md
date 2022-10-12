# The next-generation genomic feature format and ngfftools

The ngfftools is a suit of tools for manipulating files with next-generation genomic feature format (NGFF), and contains various functional modules including formats conversion among NGFF and other formats (GFF3, GTF etc.), retrieval subsets (such as grep, merge, extract), retrieval sequences and coordinate transformation.


## Installing and Building ngfftools

The installing and building of ngfftools is very simple:

    git clone https://github.com/dengcao3/ngfftools.git    
    cd ngfftools
    bash install.sh

Then, You can copy the resulting ngfftools executable into somewhere to your $PATH, or add /path/to/ngfftools/bin to your $PATH, or run it where it is.

## Running the tutorial

To test the installation, you can run the following codes:

    cd /path/to/ngfftools/doc/tutorial
    bash tutorial.sh

If the commands donn't throw any exceptions, the ngfftools is installed correctedly.

You also can run each command separately in the tutorial.sh file. You can find the manual file at /path/to/ngfftools/doc/ngfftools_manual_v202210.pdf

## NGFF format

You can find the format specification file at /path/to/ngfftools/doc/ngff_format_specifications_v202210.pdf

## Citing

Please cite this paper when using NGFF format and/or ngfftools for your publications:

    Cao Deng, et al. The next-generation genomic feature format and ngfftools.
