#!/bin/bash
a=$(wc -l < *.xyz)                 # 'wc -l' = line count. Finds the number of lines in your main .xyz file and outputs the value.
b=$(expr $a - 2)                   # Used to not account for the first 2 lines.  
c=$(expr $b / 3)                   # Divides the number of atoms in your main .xyz file by 3 to give the loop range.
i=3                                # Start at line 3.
k=5                                # Start at aline 5.
for q in $(seq 1 $c); do           # "seq" is used to bypass not being able to use dollar sign in a loop, e.g. {1..$a}.

# >> is used in bash to append an output to the end of a file.
# head -2 echos the first 2 lines of a file.
# "sed -n ''$i','$k'p" ../prism.xyz" copies lines i - k in the .xyz file. 
# tail -18 echos the last 18 lines of a file.
# "sed -i ''$i','$k'd' water_$q.xyz" deletes lines i - k in the .xyz file.
# "sed -n '/.xyz/!p' ../input.dat" deletes any line in the file that contains '.xyz'. 
# sed -i automatically edits and saves the changes done to a file.
# sed -n echos the resulting file and does not save the changes. 
  mkdir water_$q
  cd water_$q
  sed -n '/.xyz/!p' ../input.dat > water_input.dat 
  echo 'water_'$q'.xyz' >> water_input.dat
  head -2 ../*.xyz>> water_$q.xyz           
  sed -n ''$i','$k'p' ../*.xyz >> water_$q.xyz
  tail -18 ../*.xyz >> water_$q.xyz 
  i=$(expr $i + 3)
  k=$(expr $k + 3)
  sed -i ''$i','$k'd' water_$q.xyz 
  cd ../
done
