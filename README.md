# Bioinfo-tools

## interlace.py

Tool used to merge data matrix by column.

Example:

| RowName | Nb | Occ |
|---------|----|-----|
| Foo     | 4  | 4   |
| Bar     | 3  | 3   |
| Foobar  | 3  | 3   |

and 

| RowName | Nb | Occ |
|---------|----|-----|
| Foo     | 8  | 1   |
| Bar     | 2  | 4   |
| Foobar  | 2  | 7   |

Will generate:

| RowName | Nb1 | Occ1 | Nb2 | Occ2 |
|---------|-----|------|-----|------|
| Foo     | 4   | 4    | 8   | 1    |
| Bar     | 3   | 3    | 2   | 4    |
| Foobar  | 3   | 3    | 2   | 7    |

Warning: data matrix must be row-name sorted. The program won't index rownames and it will assume that they always come in the same order.
Warning2: this program isn't optimized ... it will consume a LOT of memory if you are playing with a lot of files

## pyVCF.py

Tool used in order to synchronize variants from two or more VCF files. Adapted to plug before Sciclone or QuantumClone.
Reads VCF files from any Tumor or Normal sample, looks up in BAM files for matching variants in any other sample, and outputs a result for every sample.

## Annotate_custom_rows.py

Tool used to add genomic information stored in a flat file database bgziped and tabix indexed to a generic other file containing genomic information.
Both files must be TSV-type files, though other separators may be really easy to implement.
Takes in input two files (the one that you want to annotate and the database), indexes (the indexes of the genomic regions for the file that you want to annotate and the index of the columns that you want to
add from your database). The program may remove trailing "chr" in order to ensure compatibility for all databases (in the same way, adding chr would be trivial, but I had no use for that).

Ex: python Annotate_custom_rows.py -i svdetect.txt -o - -c 4-5,7-8 --db ./genes.gz --remove-trailing-chr --db-index 2,3
=> This will take as genomic region input columns 4-5 and columns 7-8. The program will make it out for most situations, like chr3 in a first column, and 1:2 in a second column.


## remove_annot_*.py

Those scripts helps filtering snpEff annotated (either EFF format or new ANN format) vcf files.
They support either manual curation or batch removal of certain annotations.

## try_mutload.sh

Weird stuff that is supposed to help us filtering Integragen's variants... For now, we've tried like 10 different ways of filtering those variants and we still *NEVER* find the same number they find...
Usage: bash try_mutload.sh <in.snp> <in.indels> <out.somatics>

## collapse_files.py

This program will allow you to collapse N files by columns, without interlacing.


Example:

| RowName | Nb | Occ |
|---------|----|-----|
| Foo     | 4  | 4   |
| Bar     | 3  | 3   |
| Foobar  | 3  | 3   |

and

| RowName | Nb | Occ |
|---------|----|-----|
| Foo     | 8  | 1   |
| Bar     | 2  | 4   |
| Foobar  | 2  | 7   |

Will generate:

| RowName | Nb  | Occ  | RowName | Nb  | Occ  |
|---------|-----|------|---------|-----|------|
| Foo     | 4   | 4    | Foo     | 8   | 1    |
| Bar     | 3   | 3    | Bar     | 2   | 4    |
| Foobar  | 3   | 3    | Foobar  | 2   | 7    |

I just now get that it's looking awfully like interlace.py ... 
Well, this one is very memory efficient (reads only one line per file at a time).
Warning: All files must have same number of lines
