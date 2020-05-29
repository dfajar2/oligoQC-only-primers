# oligoQC-only-primers

This script will align the primer sequences on the reference genome and 
predict, based on distance and orientation, potential product amplifications. 

By default, the maximum amplicon size is 2kb (which can be modified).


USAGE:
    $ /path/to/script.sh OPTIONS

        [ -p Primers fasta file ]
        [ -r Reference file Blast DB ]
        [ -o Output directory  ]
        [ -m Primers distance apart to merge blast hits
             as potential amplicons. Default: 2000 ]
        [ -v verbose ]
