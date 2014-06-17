import sys
sys.path.append('../src/')
import xltable

### data labels
ddir='data/'
odir='output/'
chain_to_names_map={"A":"Rpb1",
                    "B":"Rpb2",
                    "C":"Rpb3",
                    "D":"Rpb4",
                    "E":"Rpb5",
                    "F":"Rpb6",
                    "G":"Rpb7",
                    "H":"Rpb8",
                    "I":"Rpb9",
                    "J":"Rpb10",
                    "K":"Rpb11",
                    "L":"Rpb12"}

prot_list=["Rpb1","Rpb2","Rpb3","Rpb4","Rpb5","Rpb6","Rpb7","Rpb8","Rpb9","Rpb10","Rpb11","Rpb12"]

### loading data
xlt=xltable.XLTable()
xlt.load_pdb_coordinates(ddir+"1WCM.pdb",chain_to_names_map)
for chainname in chain_to_names_map:
    xlt.load_sequence_from_fasta_file(fasta_file=ddir+"1WCM.fasta.txt",
                                      id_in_fasta_file="1WCM:"+chainname+"|PDBID|CHAIN|SEQUENCE",
                                      protein_name=chain_to_names_map[chainname])

### creating contact map
xlt.setup_contact_map(protein_list=prot_list,upperbound=20)
xlt.plot_table(prot_listx=prot_list,
               prot_listy=prot_list,
               alphablend=0.4,
               scale_symbol_size=1.0,
               gap_between_components=100,
               filename=odir+"contacts",
               contactmap=True,
               display_residue_pairs=False)
