#!/bin/bash

BeginFlag="\[Atoms\]"
EndFlag="\[STO\]"

declare -i Bnum
declare -i Enum
declare -i nums

Bnum=$(grep -n "$BeginFlag" ./molden.dat | cut -d: -f1)
Enum=$(grep -n "$EndFlag" ./molden.dat | cut -d: -f1)
nums=$(($Enum-$Bnum-1))


file_out=$1'.gjf'

echo '%nprocs=28' >> $file_out
echo '%mem=30GB' >> $file_out
echo '# b3lyp def2svp opt freq' >> $file_out
echo -e >> $file_out
echo $1 >> $file_out
echo -e >> $file_out
echo '0 1 ' >> $file_out
grep -A $nums "$BeginFlag" ./molden.dat >> $file_out
sed -i '/\[Atoms\]/d' ${file_out}
sed -i 's/    [0-9]    [0-9]    //g' ${file_out}
sed -i 's/   [0-9][0-9]    [0-9]    //g' ${file_out}
echo -e >> $file_out
echo -e >> $file_out





