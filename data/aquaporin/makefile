
all:
	make charges-gas
	make amber-input

charges-gas:
	antechamber -i trimemwide.mol2 -fi mol2 -o trimemwide-gas.mol2 -fo mol2 -c gas

amber-input:
	tleap -f leap.cmd
