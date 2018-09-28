import RNA
from collections import Counter
#from SeedFeatures import *
from InteractionRichPresentation import *

class RichDuplex(object):

    def __init__(self, o_mir, o_mrna):
        self.mir = o_mir
        self.mrna = o_mrna
        self.duplex = RNA.duplexfold(self.mir, self.mrna)
        self.duplex_score = -self.duplex.energy

        (mir_pairing, mrna_pairing) = self.duplex.structure.split('&')
        self.mir_coor = (self.duplex.i - len(mir_pairing), self.duplex.i)
        self.mrna_coor = (self.duplex.j - 1, self.duplex.j + len(mrna_pairing) - 1)
        self.mir_idx = self.find_pairing(mir_pairing, '(')
        self.mrna_idx = self.find_pairing(mrna_pairing, ')')
        self.active_mir = self.mir[self.mir_coor[0]:self.mir_coor[1]]
        self.active_mrna = self.mrna[self.mrna_coor[0]:self.mrna_coor[1]]
        self.mrna_idx = self.mrna_idx[::-1]
        self.mir_len = len(self.active_mir)
        self.mrna_len = len(self.active_mrna)

        self.IRP = self.parse_interaction()


    def find_pairing(self, s, ch):
        return [i for i, ltr in enumerate(s) if ltr == ch]

    def parse_interaction (self):
        mir = self.active_mir
        mrna = self.active_mrna

        mrna_bulge = ""
        mrna_inter = ""
        mir_inter = ""
        mir_bulge = ""
        mir_i = 0
        mrna_i = self.mrna_len - 1
        if (self.mir_coor[0] > 0):
            mir_bulge+=self.mir[:self.mir_coor[0]]
            mir_inter+=" "*self.mir_coor[0]
            mrna_inter+=" "*self.mir_coor[0]
            mrna_bulge_additon=self.mrna[self.mrna_coor[1]:self.mrna_coor[1]+self.mir_coor[0]]
            mrna_bulge+=mrna_bulge_additon[::-1]


        for i in range(len(self.mir_idx)) :
            #deal with the bulge
            mir_bulge_idx = range (mir_i, self.mir_idx[i])
            mir_bulge+=mir[mir_i:self.mir_idx[i]]
            mrna_bulge_idx= range (mrna_i, self.mrna_idx[i], -1)
            mrna_bulge+=mrna[mrna_i:self.mrna_idx[i]: -1]
            c_pos = max (len(mrna_bulge_idx), len(mir_bulge_idx))
            mrna_inter += " " * c_pos
            mir_inter += " " * c_pos
            mrna_bulge+=" "*(c_pos - len(mrna_bulge_idx))
            mir_bulge+=" "*(c_pos - len(mir_bulge_idx))
            #deal with the interaction
            mir_bulge+=" "
            mir_inter+=mir[self.mir_idx[i]]
            mrna_bulge+=" "
            mrna_inter+=mrna[self.mrna_idx[i]]
            #update the idx
            mir_i=self.mir_idx[i] + 1
            mrna_i= self.mrna_idx[i] - 1
        #deal with the tail
        if (mir_i<=len(mir)) :
            mir_bulge+=mir[mir_i:]
        if (mrna_i >=0) :
            mrna_bulge+=mrna[mrna_i::-1]

        return InteractionRichPresentation (mrna_bulge, mrna_inter, mir_inter, mir_bulge)


    def tostring(self):
        classstr = ""
        classstr += " {} \n".format(self.duplex.structure)
        classstr += self.IRP.__str__()

        classstr = classstr + "active_mrna[-1] ({}): {} \n".format(len(self.active_mrna[::-1]), self.active_mrna[::-1])
        classstr = classstr + "active_mir      ({}): {} \n".format(len(self.active_mir), self.active_mir)
        classstr = classstr + "full_mir        ({}): {} \n".format(len(self.mir), self.mir)
        classstr = classstr + "offset = {} \n".format(self.mir_coor[0])





        #s = SeedFeatures(self)
        #classstr = classstr + "canonic = {} \n".format(s.canonic_seed)

   #     mir_seed, pairs_in_seed = self.extract_seed()
    #    classstr = classstr + "mir seed = {} \n".format(mir_seed)
     #   classstr = classstr + "pairs_in_seed = {} num_of_pairs = {} \n".format(pairs_in_seed, len(pairs_in_seed))

        return classstr


    def __str__(self):
        return self.tostring()
