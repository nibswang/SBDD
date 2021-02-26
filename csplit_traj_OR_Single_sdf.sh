####traj.pdb is trajtory PDBs from MD, 'my.' is prefix of generated PDB files 
a=`grep ENDMDL traj.pdb | wc -l`
b=`expr $a - 2`
csplit -k -s -n 3 -f my. traj.pdb '/^ENDMDL/+1' '{'$b'}'

##rename
list=`ls my.*`
for i in $list
 do echo $i
 mv $i $i.pdb
 done


####
csplit Enamine_PriAmines_Aliphatic_Policyclic_979cmpds_20200331.sdf '/^\$\$/+1' -n3 -s {*} -f prefix
