package main

import (
"fmt"
"os"
"github.com/rmera/gochem"
"github.com/rmera/scu"
)

//Usage:


func main(){
	if len(os.Args)<3{
		fmt.Println("Usage: ",os.Args[0], "test.xyz template.xyz  \"indexes for test\"  \"indexes for template\" \n The indexes, given between quotations, are a list of integers separated by spaces which contains the atoms that will be used from each molecule to calculate the rotation matrix. If one or both indexes lists is empty (\"\"), the whole molecule will be used." )
		panic(":-)") //Not very subtle but oh well.
	}
	test,_:=chem.XYZRead(os.Args[1])
	templa,_:=chem.XYZRead(os.Args[2])
	testlist,_:=scu.IndexStringParse(os.Args[3])
	templalist,_:=scu.IndexStringParse(os.Args[4])
	test2,err:=chem.Super(test.Coords[0],templa.Coords[0],testlist,templalist)
	if err!=nil{
		panic(err)
	}
	chem.XYZWrite("superimposed.xyz",test2,test)
	rmsd,err:=chem.RMSD(test2,templa.Coords[0])
	fmt.Println(rmsd)
}
