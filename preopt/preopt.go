
package main

import (
	"fmt"
	"flag"
	"github.com/rmera/gochem"
	"github.com/rmera/gochem/qm"
	"strings"
	"os"
)
type qmdef interface{
	qm.Handle
	SetDefaults()
}

func main() {
	//flag parsing
	//Lots of options, but the defaults are sane enough that you shouldnt need to use more than the filename and qm program.
	charge:=flag.Int("charge",0,"The charge of the system.")
	multi:=flag.Int("multi",1,"The multiplicity of the system.")
	filename:=flag.String("file","file.xyz","The XYZ file containing the coordinates for the system.")
	functional:=flag.String("func", "BP86", "The density functional used. TPSS and BP86 activate RI when possible.")
	program:=flag.String("program","nwchem", "The QM program used: qcmine, nwchem or orca.")
	basis:=flag.String("basis", "def2-SVP","the basis set to use. Use Karlsruhe basis.")
	dielectric:=flag.Float64("epsilon", -1, "The dielectric constant. -1 indicates no dielectric used.")
	optimize:=flag.Bool("opt",true,"Wether to optimize or run an SP calculation.")
	flag.Parse()
	//We set the calculation to the values in the flags. Not big deal.
	mol, err := chem.XYZRead(*filename)
	if err != nil {
		panic(err.Error())
	}
	mol.SetMulti(*multi)
	mol.SetCharge(*charge)
	calc := new(qm.Calc)
	if strings.Contains("TPSS,BP86,PBE",*functional){
		calc.RI=true
	}
	calc.Memory = 1000
	calc.Dielectric = *dielectric
	calc.Grid = 3
	if !strings.HasPrefix(*basis,"def2"){
		fmt.Println("Told ya to use Karlsruhe basis! RI cannot be activated. Happy now?")
		calc.RI=false
	}
	calc.Basis = *basis
	calc.Optimize=*optimize
	if calc.Optimize{
		calc.SCFTightness=1
	}
	calc.Disperssion="D3"
	//The preeliminar optimization will have a different name.
	namep:=strings.Replace(*filename,".xyz","OPT",-1)
	name:=strings.Replace(*filename,".xyz","",-1)
	pre:=qm.NewMopacHandle()
	pre.SetName(namep)
	pre.BuildInput(mol.Coords[0],mol,calc)
	pre.Run(true)
	ncoords,err:=pre.OptimizedGeometry(mol)
	if err!=nil{
		panic(err.Error())
	}
	//MOPAC will overwrite the method, as it doesnt support DFT. This is why we set it AFTER the preeliminar calculations
	calc.Method=*functional
	var QM qmdef
	switch *program{
		default:
			QM=qmdef(qm.NewNWChemHandle())
		case "orca":
			QM=qmdef(qm.NewOrcaHandle())
//		default:
//			//There is a hicup with ChemShell since I have not implemented a few methods.
//			//I should have an interface in goChem that does not require Energy() and OptimizedGeometry()
//			CS:=qm.NewCSHandle()
//			CS.SetDefaults()
//			CS.SetName(name)
//			CS.BuildInput(ncoords,mol,calc)
	}
	if QM!=nil{
		QM.SetDefaults()
		QM.SetName(name)
		QM.BuildInput(ncoords,mol,calc)
	}
	fmt.Fprintln(os.Stderr,"Apparently, we made it :-)")
}

