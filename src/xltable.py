from math import sqrt
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
import numpy as np
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import utilities
from collections import defaultdict

class XLTable():
    ''' class to read, analyze, and plot xlink data on contact maps
    Canonical way to read the data:
    1) load sequences and name them
    2) load coordinates for those sequences from PDB file
    3) add crosslinks
    4) create contact map
    5) plot
    '''

    def __init__(self):
        self.sequence_dict={}
        self.field_map={}
        self.cross_link_db=[]
        self.residue_pair_list=[]          # list of special residue pairs to display
        self.distance_maps=[]              # distance map for each copy of the complex
        self.index_dict=defaultdict(list)  # location in the dmap of each residue

        # internal things
        self._first=True
    def _colormap(self, dist, threshold=35, tolerance=0):
        if dist < threshold - tolerance:
            return "Green"
        elif dist >= threshold + tolerance:
            return "Orange"
        else:
            return "Orange"

    def _get_distance(self,r1,c1,r2,c2):
        print '1:',c1,r1,len(self.index_dict[c1]),'2:',c2,r2,len(self.index_dict[c2])

        idx1=self.index_dict[c1][r1]
        idx2=self.index_dict[c2][r2]
        return self.av_dist_map[idx1,idx2]

    def load_sequence_from_fasta_file(self,fasta_file,id_in_fasta_file,protein_name):
        ''' read sequence. structures are displayed in the same order as sequences are read.
        fasta_file:        file to read
        id_in_fasta_file:  id of desired sequence
        protein_name:      identifier for this sequence (use same name when handling coordinates)
        can provide the fasta name (for retrieval) and the protein name (for storage) '''
        handle = open(fasta_file, "rU")
        record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        handle.close()
        if id_in_fasta_file is None:
            id_in_fasta_file = name
        try:
            length = len(record_dict[id_in_fasta_file].seq)
        except KeyError:
            print "add_component_sequence: id %s not found in fasta file" % id_in_fasta_file
            exit()
        self.sequence_dict[protein_name] = str(record_dict[id_in_fasta_file].seq).replace("*", "")

    def load_pdb_coordinates(self,pdbfile,chain_to_name_map):
        ''' read coordinates from a pdb file.
        pdbfile:             file for reading coords
        chain_to_name_map:   correspond chain ID with protein name (will ONLY read these chains)
        This function returns an error if the sequence for each chain has NOT been read

        '''
        pdbparser = PDBParser()
        structure = pdbparser.get_structure(pdbfile,pdbfile)
        total_len = sum(len(self.sequence_dict[s]) for s in self.sequence_dict)
        coords = np.ones((total_len,3)) * 1e5 #default to coords "very far away"
        prev_stop=0
        for n,model in enumerate(structure):
            for cid in chain_to_name_map:
                cname=chain_to_name_map[cid]
                if cname not in self.sequence_dict:
                    print "ERROR: chain",cname,'has not been read or has a naming mismatch'
                    return
                if self._first:
                    self.index_dict[cname]=range(prev_stop,prev_stop+len(self.sequence_dict[cname]))
                for residue in model[cid]:
                    if "CA" in residue:
                        ca=residue["CA"]
                        rnum=residue.id[1]
                        coords[rnum+prev_stop-1,:]=ca.get_coord()
                    #else:
                    #    print residue
                prev_stop+=len(self.sequence_dict[cname])
        dists = cdist(coords, coords)
        self.distance_maps.append(dists)
        if self._first:
            self._first=False

    def load_crosslinks(self,crosslinkfile,field_map):
        ''' read crosslinks from a CSV file.
        provide a dictionary to explain the columns
        (must contain prot1,prot2,res1,res2,score)'''
        if len(set(field_map.keys()) & set(("prot1","prot2","res1","res2","score")))<5:
            return
            print "ERROR: your field_map dictionary does not contain required fields"
        self.cross_link_db=utilities.get_db_from_csv(crosslinkfile)
        self.field_map=field_map

    def set_residue_pairs_to_display(self,residue_type_pair):
        ''' select the atom names of residue pairs to plot on the contact map
        list of residues types must be single letter code
        e.g. residue_type_pair=("K","K")
        '''
        rtp=sorted(residue_type_pair)
        for prot1 in self.sequence_dict:
            seq1=self.sequence_dict[prot1]
            for nres1,res1 in enumerate(seq1):
                for prot2 in self.sequence_dict:
                    seq2=self.sequence_dict[prot2]
                    for nres2,res2 in enumerate(seq2):
                        if sorted((res1,res2))==rtp:
                           self.residue_pair_list.append((nres1+1,prot1,nres2+1,prot2))

    def setup_contact_map(self,upperbound=20):
        ''' loop through each distance map and get frequency of contacts
        upperbound:   maximum distance to be marked
        '''

        # filter each distance and get the frequency of contact
        print 'filtering distance maps'
        all_dists = np.dstack(self.distance_maps)
        self.av_dist_map=1.0/len(self.distance_maps) * np.sum(all_dists,axis=2)

        binary_dists = np.where((all_dists <= upperbound) & (all_dists >= 1.0), 1.0, 0.0)
        if all_dists.shape[2]>1:
            self.contact_map = 1.0/len(self.distance_maps) * np.sum(binary_dists,axis=2)
        else:
            self.contact_map = binary_dists[:,:,0]

    def plot_table(self, prot_listx=None,
                   prot_listy=None,
                   no_dist_info=False,
                   confidence_info=False,
                   filter=None,
                   display_residue_pairs=False,
                   contactmap=False,
                   filename=None,
                   confidence_classes=None,
                   alphablend=0.1,
                   scale_symbol_size=1.0,
                   gap_between_components=0):
        ''' plot the xlink table with optional contact map.
        prot_listx:             list of protein names on the x-axis
        prot_listy:             list of protein names on the y-axis
        no_dist_info:           plot only the cross-links as grey spots
        confidence_info:
        filter:                 list of tuples to filter on. each one contains:
                                    keyword in the database to be filtered on
                                    relationship ">","==","<"
                                    a value
                                example ("ID_Score",">",40)
        display_residue_pairs:  display all pairs defined in self.residue_pair_list
        contactmap:             display the contact map
        filename:               save to file (adds .pdf extension)
        confidence_classes:
        alphablend:
        scale_symbol_size:      rescale the symbol for the crosslink
        gap_between_components:
        '''

        # prepare figure
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)
        ax.set_xticks([])
        ax.set_yticks([])

        # set the list of proteins on the x axis
        if prot_listx is None:
            prot_listx = self.sequence_dict.keys()
            prot_listx.sort()
        nresx = gap_between_components + \
            sum([len(self.sequence_dict[name])
                + gap_between_components for name in prot_listx])

        # set the list of proteins on the y axis
        if prot_listy is None:
            prot_listy = self.sequence_dict.keys()
            prot_listy.sort()
        nresy = gap_between_components + \
            sum([len(self.sequence_dict[name])
                + gap_between_components for name in prot_listy])

        # this is the residue offset for each protein
        resoffsetx = {}
        resendx = {}
        res = gap_between_components
        for prot in prot_listx:
            resoffsetx[prot] = res
            res += len(self.sequence_dict[prot])
            resendx[prot] = res
            res += gap_between_components

        resoffsety = {}
        resendy = {}
        res = gap_between_components
        for prot in prot_listy:
            resoffsety[prot] = res
            res += len(self.sequence_dict[prot])
            resendy[prot] = res
            res += gap_between_components

        resoffsetdiagonal = {}
        res = gap_between_components
        for prot in utilities.OrderedSet(prot_listx + prot_listy):
            resoffsetdiagonal[prot] = res
            res += len(self.sequence_dict[prot])
            res += gap_between_components

        # plot protein boundaries
        xticks = []
        xlabels = []
        for n, prot in enumerate(prot_listx):
            res = resoffsetx[prot]
            end = resendx[prot]
            for proty in prot_listy:
                resy = resoffsety[proty]
                endy = resendy[proty]
                ax.plot([res, res], [resy, endy], 'k-', lw=0.4)
                ax.plot([end, end], [resy, endy], 'k-', lw=0.4)
            xticks.append((float(res) + float(end)) / 2)
            xlabels.append(prot)

        yticks = []
        ylabels = []
        for n, prot in enumerate(prot_listy):
            res = resoffsety[prot]
            end = resendy[prot]
            for protx in prot_listx:
                resx = resoffsetx[protx]
                endx = resendx[protx]
                ax.plot([resx, endx], [res, res], 'k-', lw=0.4)
                ax.plot([resx, endx], [end, end], 'k-', lw=0.4)
            yticks.append((float(res) + float(end)) / 2)
            ylabels.append(prot)

        # plot the contact map
        if contactmap:
            tmp_array = np.zeros((nresx, nresy))
            for px in prot_listx:
                for py in prot_listy:
                    resx = resoffsety[px]
                    lengx = resendx[px] - 1
                    resy = resoffsety[py]
                    lengy = resendy[py] - 1
                    indexes_x = self.index_dict[px]
                    minx = min(indexes_x)
                    maxx = max(indexes_x)
                    indexes_y = self.index_dict[py]
                    miny = min(indexes_y)
                    maxy = max(indexes_y)
                    tmp_array[resx:lengx,resy:lengy] = self.contact_map[minx:maxx,miny:maxy]

            ax.imshow(tmp_array,
                      cmap=cm.binary,
                      origin='lower',
                      interpolation='nearest')

        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels, rotation=90)
        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels)

        # set the crosslinks
        already_added_xls = []
        for xl in self.cross_link_db:
            r1=int(xl[self.field_map["res1"]])
            c1=xl[self.field_map["prot1"]]
            r2=int(xl[self.field_map["res2"]])
            c2=xl[self.field_map["prot2"]]
            score=float(xl[self.field_map["score"]])

            try:
              mdist=self._get_distance(r1,c1,r2,c2)
              color = self._colormap(mdist)
            except KeyError:
              color="gray"

            try:
                pos1 = r1 + resoffsetx[c1]
            except:
                continue
            try:
                pos2 = r2 + resoffsety[c2]
            except:
                continue

            '''
            if not filter is None:
                xldb = self.external_csv_data[unique_identifier]
                xldb.update({"Distance": mdist})
                xldb.update({"Distance_stdv": stdv})

                if filter[1] == ">":
                    if float(xldb[filter[0]]) <= float(filter[2]):
                        continue

                if filter[1] == "<":
                    if float(xldb[filter[0]]) >= float(filter[2]):
                        continue

                if filter[1] == "==":
                    if float(xldb[filter[0]]) != float(filter[2]):
                        continue
            '''

            # everything below is used for plotting the diagonal
            # when you have a rectangolar plots
            pos_for_diagonal1 = r1 + resoffsetdiagonal[c1]
            pos_for_diagonal2 = r2 + resoffsetdiagonal[c2]
            if confidence_info:
                if confidence == '0.01':
                    markersize = 14 * scale_symbol_size
                elif confidence == '0.05':
                    markersize = 9 * scale_symbol_size
                elif confidence == '0.1':
                    markersize = 6 * scale_symbol_size
                else:
                    markersize = 15 * scale_symbol_size
            else:
                markersize = 5 * scale_symbol_size

            ax.plot([pos1],
                    [pos2],
                    'o',
                    c=color,
                    alpha=alphablend,
                    markersize=markersize)

            ax.plot([pos2],
                    [pos1],
                    'o',
                    c=color,
                    alpha=alphablend,
                    markersize=markersize)

        # plot requested residue pairs
        if display_residue_pairs:
           for rp in self.residue_pair_list:
               r1=rp[0]
               c1=rp[1]
               r2=rp[2]
               c2=rp[3]

               try:
                  dist=self._get_distance(r1,c1,r2,c2)
               except:
                  continue

               if dist<=40.0:
                 print rp
                 try:
                    pos1 = r1 + resoffsetx[c1]
                 except:
                    continue
                 try:
                    pos2 = r2 + resoffsety[c2]
                 except:
                    continue

                 ax.plot([pos1],
                         [pos2],
                         '+',
                         c="blue",
                         alpha=0.1,
                         markersize=markersize)

        # zadisplay and write to file
        fig.set_size_inches(0.002 * nresx, 0.002 * nresy)
        [i.set_linewidth(2.0) for i in ax.spines.itervalues()]
        if filename:
            plt.savefig(filename + ".pdf", dpi=300,transparent="False")
        plt.show()
