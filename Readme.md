#


## Scripts

### Check_Restriction_site

### Create_Haplotype_Structure

### Filter_PCR_duplicates

### Filter_Reads

>Usage: FilterReads.py [options] input.fasta [output.fasta]  

>Options:  
>  -h, --help            show this help message and exit  
>  -q, --quiet           do not print status messages to the screen  
>  -u, --fastq           input file is fastq  
>  -z, --gzip            input file is gzipped  
>  -l X, --min-length=X  write only sequence with lengths at least X  
>  -i X, --id-list=X     write only sequence with an ID from this list. List  
>                        can be comma separated string of IDs or a path to a  
>                        file with a line separated list of IDs  
>  -r X, --random=X      randomly sample X sequence from input file  
>  -e, --regexp          use regular expression instead of exact matching for  
>                        IDs  
>  -a, --ignore-at       ignore the first letter of the query IDs if it is an @
>                        (this is for more convenient filter list creation from
>                        fastq files)  
>  -n, --negative        do exactly the opposite of what would normally be done  

Filter fasta/fastq files in different way.

1. Filter by minimum length (`-l`)  
Only write sequences of a certain length to the output.
2. Filter by list of IDs (`-i`)  
Only write sequences with an ID from the list given. The list can either be given as a comma separated list of IDs or as in a file with one ID per line. With the `-e` option the given "IDs" will be used as regular expression instead of exact matching (using python regex).
3. "Filter" randomly (`-r`)  
Write a random subset of the sequences of the given size.

The `-n` option will switch to negative mode. Meaning the script will do exactly the opposite of what it normally does.

The `-a` option is will make the script ignore @-signs at the begining of IDs in the ID list.
The main use case for this is with two fastq files (A.fq and B.fq) and all your ID lines start with @M01271 (because M01271 is the serial number of your sequencer). If we want to keep in B only the sequences that are also in A, we can run the following:
>grep "^@M01271" A.fq > id_list.txt  
>python filterFasta.py -i id_list.txt -a B.fq > B_and_A.fq

Input data can be provided as a file (first argument) or be piped in.