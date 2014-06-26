import sys
sys.path.append('../src/')
import xltable

### data labels
ddir='data/'
odir='./'
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

field_map={"prot1":"pep1.accession",
           "prot2":"pep2.accession",
           "res1":"pep1.xlinked_aa",
           "res2":"pep2.xlinked_aa",
           "score":"SVM dval"}

prot_list=["Rpb1","Rpb2","Rpb3","Rpb4","Rpb5","Rpb6","Rpb7","Rpb8","Rpb9","Rpb10","Rpb11","Rpb12"]

### loading data
xlt=xltable.XLTable(contact_threshold=20)
for chainname in chain_to_names_map:
        xlt.load_sequence_from_fasta_file(fasta_file=ddir+"1WCM.fasta.txt",
                                      id_in_fasta_file="1WCM:"+chainname+"|PDBID|CHAIN|SEQUENCE",
                                      protein_name=chain_to_names_map[chainname])
xlt.load_pdb_coordinates(ddir+"1WCM.pdb",chain_to_names_map)
xlt.load_crosslinks(ddir+"polii_xlinks.csv",field_map)
#xlt.set_residue_pairs_to_display(("K","K"))

### creating contact map
xlt.setup_contact_map()

### plotting
xlt.plot_table(prot_listx=prot_list,
               prot_listy=prot_list,
               alphablend=0.4,
               scale_symbol_size=1.0,
               gap_between_components=100,
               filename=odir+"model_with_xl",
               contactmap=True,
               display_residue_pairs=False)
