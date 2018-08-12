#!/bin/bash
a=$(wc -l < *.xyz)                 # 'wc -l' = line count. Finds the number of lines in your main .xyz file and outputs the value.
b=$(expr $a - 2)                   # Used to not account for the first 2 lines.  
c=$(expr $b / 3)                   # Divides the number of atoms in your main .xyz file by 3 to give the loop range.
for q in $(seq 1 $c); do           # "seq" is used to bypass not being able to use dollar sign in a loop, e.g. {1..$a}.

  cd water_$q
   ../../src/a.out < water_input.dat > aout &
   cd ../
done
