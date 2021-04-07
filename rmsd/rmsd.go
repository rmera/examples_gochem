package main

import (
	"fmt"
	"os"
	"strings"

	chem "github.com/rmera/gochem"
	"github.com/rmera/scu"
)

//Usage:

func main() {
	if len(os.Args) < 3 {
		fmt.Println("Usage: ", os.Args[0], "test.xyz template.xyz  \"indexes for test\"  \"indexes for template\" \n The indexes, given between quotations, are a list of integers separated by spaces which contains the atoms that will be used from each molecule to calculate the rotation matrix. If one or both indexes lists is empty (\"\"), the whole molecule will be used.")
		panic(":-)") //Not very subtle but oh well.
	}
	testname := os.Args[1]
	templaname := os.Args[2]
	var test, templa *chem.Molecule
	if strings.Contains(testname, ".pdb") {
		test, _ = chem.PDBFileRead(testname, false)
		templa, _ = chem.PDBFileRead(templaname, false)
	} else {
		test, _ = chem.XYZFileRead(testname)
		templa, _ = chem.XYZFileRead(testname)
	}
	testlist, _ := scu.IndexStringParse(os.Args[3])
	templalist, _ := scu.IndexStringParse(os.Args[4])
	if os.Args[len(os.Args)-1] == "-CA" {
		for i := 0; i < templa.Len(); i++ {
			if templa.Atom(i).Name == "CA" {
				templalist = append(templalist, i)
			}
		}

		for i := 0; i < test.Len(); i++ {
			if test.Atom(i).Name == "CA" {
				testlist = append(testlist, i)
			}
		}
	}
	test2, err := chem.Super(test.Coords[0], templa.Coords[0], testlist, templalist)
	if err != nil {
		panic(err)
	}
	chem.XYZFileWrite("superimposed.xyz", test2, test)
	chem.PDBFileWrite("superimposed.pdb", test2, test, nil)
	rmsd, err := chem.RMSD(test2, templa.Coords[0])
	fmt.Println(rmsd)
}
