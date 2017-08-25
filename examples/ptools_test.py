from ptools import *

pdb = Rigidbody("1cgi.pdb")
chainE = pdb.SelectChainId("E")
chainI = pdb.SelectChainId("I")
protE = chainE.CreateRigid()
protI = chainI.CreateRigid()
WritePDB(protE, "receptor.pdb")
WritePDB(protI, "ligand.pdb")