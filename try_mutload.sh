inSNP=$1
inIndel=$2
out=$3

sed "s/\t\t/\t0\t/g" $inSNP | sed "s/\tNA\t/\t0\t/g" | sed 's/chr//' | grep -P 'splice|missense|nonsense' | awk -F "\t" 'NR>0 && $3 >4 && $4="Somatic"  && $9 <= 0.5 && $11 <= 0.5 && $13 <=1 && $14<=1 && $57 >= 0.06 && $56 <=0.04 && $59 >=3{OFS="\t"; print "SNP: ", $1, $2, $15,$60, $16}' | sort -k 2,2n -k 3,3n | uniq > $out
nbSNP=$(sed "s/\t\t/\t0\t/g" $inSNP | sed "s/\tNA\t/\t0\t/g" | sed 's/chr//' | grep -P 'splice|missense|nonsense' | awk -F "\t" 'NR>0 && $3 >4 && $4="Somatic"  && $9 <= 0.5 && $11 <= 0.5 && $13 <=1 && $14<=1 && $57 >= 0.06 && $56 <=0.04 && $59 >=3{OFS="\t"; print "SNP :", $1, $2, $15,$60, $16}' | sort -k 2,2n -k 3,3n | uniq | wc -l )
sed "s/\t\t/\t0\t/g" $inIndel | sed "s/\tNA\t/\t0\t/g" | sed 's/chr//' | awk -F "\t" 'NR>0 && $3 >=5 && $10 <= 0.5 && $11 <=1 && $12<=1 && $52 >= 0.05 && $51 < 0.04 && $54 >=3 {OFS="\t"; print "Indel: ", $2, $1}' | sort -k 2,2n -k 3,3n | uniq >> $out
nbIndel=$(sed "s/\t\t/\t0\t/g" $inIndel | sed "s/\tNA\t/\t0\t/g" | sed 's/chr//' |  awk -F "\t" 'NR>0 && $3 >=5 && $10 <= 0.5 && $11 <=1 && $12<=1 && $52 >= 0.05 && $51 < 0.04 && $54 >=3 {OFS="\t"; print "Indel: ", $2, $1}' | sort -k 2,2n -k 3,3n | uniq | wc -l)
echo "SNP: $nbSNP"
echo "Indel: $nbIndel"
echo $(( $nbSNP + $nbIndel ))

