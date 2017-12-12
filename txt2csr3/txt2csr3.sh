#!/bin/sh

if [ -f "$1" ]; then 
  FIN=$1
  file=`basename $FIN .dat`
  FSED=${file}_key.sed
  FREN=Renumerated_${file}
  cut -f1 $FIN | sort -g | uniq > tmp1  #cut -d' '
  cut -f2 $FIN | sort -g | uniq > tmp2
  cat tmp1 tmp2 | sort -g | uniq | awk '{print "s/\\b"$0"\\b/"NR"/g"}' > $FSED
  sed -f $FSED < $FIN > $FREN
  echo "Renumerated matrix written to $FREN, sed commands for renumeration to $FSED"
  rm tmp1 tmp2 
  FOUT=`basename $FREN`_CSR3.dat
  # Преобразование матрицы в формат CSR3
  awk 'BEGIN {C = 0; max = 0; ORS=" "} {if($1 > C) {print 0 > "values"; print $1 > "columns"; print $3 > "values"; print $2 > "columns"; C = $1} else {print $3 > "values"; print $2 > "columns"} if($1>max) max=$1; if($2>max) max=$2; N[$1]++} END {C++; if(C<max) for(; C<=max; C++) {print 0 > "values"; print C > "columns"}; print "\n" > "values"; print "\n" > "columns"; rowIndex[1] = 1; for(i=2; i<=max; i++) {if(i-1 in N) rowIndex[i] = N[i-1] + rowIndex[i-1]; else rowIndex[i] = rowIndex[i-1];} for(i=1; i<=max; i++) print rowIndex[i] > "rowIndex"; print NR+1 > "rowIndex"; print "\n" > "rowIndex"}' $FREN
  cat values columns rowIndex > $FOUT
  echo "Matrix in CSR3 format written to $FOUT"
  rm values columns rowIndex
fi

