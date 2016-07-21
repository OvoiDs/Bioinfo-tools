# Bioinfo-tools

# interlace.py

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
