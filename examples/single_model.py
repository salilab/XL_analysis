import sys,os
sys.path.append('../src/')
import xltable

### data labels
ddir='data/'
odir='./'
chain_to_names_map={"A":"Spc97",
                    "B":"Spc98",
                    "G":"Spc97A"}
prot_list=["Spc97","Spc98","Spc97A"]

### reading in sequences and coordinates
xlt=xltable.XLTable(contact_threshold=20)
xlt.load_sequence_from_fasta_file(ddir+'yGCP2_full.fasta',
                                  'GCP2_YEAST',
                                  'Spc97')
xlt.load_sequence_from_fasta_file(ddir+'yGCP3_full.fasta',
                                  'GCP3_YEAST',
                                  'Spc98')
xlt.load_sequence_from_fasta_file(ddir+'yGCP2_full.fasta',
                                  'GCP2_YEAST',
                                  'Spc97A')
xlt.load_pdb_coordinates(ddir+"tusc/tusc_flex0.pdb",chain_to_names_map)

### set up contact map and plot
xlt.setup_contact_map()
xlt.plot_table(prot_listx=prot_list,
               prot_listy=prot_list,
               alphablend=0.4,
               scale_symbol_size=1.0,
               gap_between_components=100,
               filename=odir+"single_model",
               contactmap=True,
               display_residue_pairs=False)
