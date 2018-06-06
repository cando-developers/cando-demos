loadAmberParams "gaff.dat"
x = loadMol2 "trimemwide-water-carbonyl.mol2"
loadamberparams "trimemwide.frcmod"
saveAmberParm x "trimemwide-water-carbonyl.prmtop" "trimemwide-water-carbonyl.inpcrd"
quit
