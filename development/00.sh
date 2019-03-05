#!/bin/bash
a=$(wc -l < tip4p-prism.xyz)   # Finds the number of lines in your main .xyz file and outputs it to a variable.
b=$(expr $a - 2)
c=$(expr $b / 3)         # Divides the number of atoms in your main .xyz file by 3 to give the loop range.
i=3
k=5
for q in $(seq 1 $c); do # "seq" is used to bypass not being able to use dollar sign in a loop, e.g. {1..$a}.
                         # "sed -n ''$i','$k'p" ../prism.xyz" copies lines i - k in the .xyz file and output to new file.
  mkdir water_$q
  cd water_$q
  sed -n '1,2p' ../tip4p-prism.xyz >> water_$q.xyz
  sed -n ''$i','$k'p' ../tip4p-prism.xyz >> water_$q.xyz
  i=$(expr $i + 3)
  k=$(expr $k + 3)
  cd ../
done
