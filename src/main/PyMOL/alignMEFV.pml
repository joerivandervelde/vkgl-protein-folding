# repair mutated structure first?
# cd /Users/joeri/git/vkgl-protein-folding/data/MEFV/tmp/NA766H
# /Applications/FoldX/foldx5MacStd/foldx_20231231 --command=RepairPDB --pdb=AF-O15553-F1-model_v4_Repair_1.pdb

#Pymol:

reinitialize
load /Users/joeri/git/vkgl-protein-folding/data/MEFV/AF-O15553-F1-model_v4.pdb, MEFV
load /Users/joeri/git/vkgl-protein-folding/data/MEFV/tmp/NA766H/AF-O15553-F1-model_v4_Repair_1_Repair.pdb, MEFV_N766H
align MEFV, MEFV_N766H
