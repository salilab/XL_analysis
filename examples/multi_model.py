import sys,os
sys.path.append('../src/')
import xltable

for i in range(1,20):
    os.system('cp /Users/cgreen/projects/gtusc/promals4/model_sse2/models/%i/tusc_flex.pdb data/tusc/tusc_flex%i.pdb'%(i,i))

### data labels
ddir='data/'
odir='output/'
chain_to_names_map={"A":"Spc97",
                    "B":"Spc98",
                    "G":"Spc97A"}

prot_list=["Spc97","Spc98","Spc97A"]


xlt=xltable.XLTable()
xlt.load_sequence_from_fasta_file(ddir+'yGCP2_full.fasta',
                                  'GCP2_YEAST',
                                  'Spc97')
xlt.load_sequence_from_fasta_file(ddir+'yGCP3_full.fasta',
                                  'GCP3_YEAST',
                                  'Spc98')
xlt.load_sequence_from_fasta_file(ddir+'yGCP2_full.fasta',
                                  'GCP2_YEAST',
                                  'Spc97A')
for i in range(1,20):
    xlt.load_pdb_coordinates(ddir+"tusc/tusc_flex%i.pdb"%i,chain_to_names_map)

xlt.setup_contact_map(upperbound=20)
xlt.plot_table(prot_listx=prot_list,
               prot_listy=prot_list,
               alphablend=0.4,
               scale_symbol_size=1.0,
               gap_between_components=100,
               filename=odir+"gtusc",
               contactmap=True,
               display_residue_pairs=False)
