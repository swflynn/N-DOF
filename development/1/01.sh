#!/bin/bash
a=$(wc -l < tip4p-prism.xyz)       # 'wc -l' = line count. Finds the number of lines in your main .xyz file and outputs the value.
b=$(expr $a - 2)                   # To not account for the first 2 lines.  
c=$(expr $b / 3)                   # Divides the number of atoms in your main .xyz file by 3 to give the loop range.
i=3                                # Start at line 3.
k=5                                # Start at line 5.
for q in $(seq 1 $c); do           # "seq" is used to bypass not being able to use dollar sign in a loop, e.g. {1..$a}.

# >> is used in bash to append an output to a file.
# head -2 echos the first 2 lines of a file.
# "sed -n ''$i','$k'p" ../prism.xyz" copies lines i - k in the .xyz file. 
# tail -18 echos the last 18 lines of a file.
# "sed -i ''$i','$k'd' water_$q.xyz" deletes lines i - k in the .xyz file.
  mkdir water_$q
  cd water_$q
  head -2 ../tip4p-prism.xyz >> water_$q.xyz           
  sed -n ''$i','$k'p' ../tip4p-prism.xyz >> water_$q.xyz
  tail -18 ../tip4p-prism.xyz >> water_$q.xyz 
  i=$(expr $i + 3)
  k=$(expr $k + 3)
  sed -i ''$i','$k'd' water_$q.xyz 
  cd ../
done
