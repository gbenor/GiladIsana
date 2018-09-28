class InteractionRichPresentation (object) :


    def __init__(self, mrna_bulge, mrna_inter, mir_inter, mir_bulge):
        self.mrna_bulge = mrna_bulge
        self.mrna_inter = mrna_inter
        self.mir_inter = mir_inter
        self.mir_bulge = mir_bulge
        self.mir_bulges_count, self.mrna_bulges_count = self.count_bulges()


    def mir_iterator (self):
        i=0
        while i<max (len(self.mir_inter), len(self.mir_bulge)) :
            if i<len(self.mir_inter) :
                if self.mir_inter[i]!=' ':
                    yield self.mir_inter[i]
            if i < len(self.mir_bulge):
                if self.mir_bulge[i] != ' ':
                    yield self.mir_bulge[i]
            i+=1

    def extract_seed (self) :
        def remove_bulge_from_mrna_start (c_seed) :
            #function to remove the option of leading bulge of mrna from the seed
            if c_seed.mrna_bulge[0]!=" " and c_seed.mrna_inter[0]==" " and c_seed.mir_inter[0]==" " and c_seed.mir_bulge[0]==" " :
                return InteractionRichPresentation (c_seed.mrna_bulge[1:], c_seed.mrna_inter[1:], c_seed.mir_inter[1:], c_seed.mir_bulge[1:])
            else:
                return c_seed

        i=7
        while True:
            i+=1
            num_of_mirna_nt = sum(c != ' ' for c in self.mir_inter[:i]) + sum(c != ' ' for c in self.mir_bulge[:i])
            if num_of_mirna_nt==8:
                break
        cur_seed = InteractionRichPresentation(self.mrna_bulge[:i], self.mrna_inter[:i], self.mir_inter[:i], self.mir_bulge[:i])
        return remove_bulge_from_mrna_start(cur_seed)

    def seed_without_bulges (self, seed) :
        mir = ""
        mrna = ""
        i = 0
        c = 0
        while (c<8):
            if (seed.mir_inter[i]!=" "):
                mir+=seed.mir_inter[i]
                mrna+=seed.mrna_inter[i]
                i+=1
                c+=1
            elif (seed.mir_bulge[i]!=" ") :
                mir+=seed.mir_bulge[i]
                mrna+=seed.mrna_inter[i]
                i += 1
                c += 1
            else:
                i+=1
        return mrna, mir



    def count_bulges (self) :
        mir_bulges_count = len(self.mir_bulge.split())
        mrna_bulges_count = len(self.mrna_bulge.split())
        return  mir_bulges_count, mrna_bulges_count

    def tostring(self):
        classstr = ""
        classstr = classstr + "target_bulge:       {}\n".format(self.mrna_bulge)
        classstr = classstr + "target_interaction: {}\n".format(self.mrna_inter)
        classstr = classstr + "mirna_interaction:  {}\n".format(self.mir_inter)
        classstr = classstr + "mirna_bulge:        {}\n".format(self.mir_bulge)
        classstr = classstr + "mrna_bulges_count: {} \nmir_bulges_count:  {}\n".format(self.mrna_bulges_count,self.mir_bulges_count)

        return classstr

    def replace_T_U (self) :
        self.mrna_bulge = self.mrna_bulge.upper().replace('T', 'U')
        self.mrna_inter = self.mrna_inter.upper().replace('T', 'U')
        self.mir_inter = self.mir_inter.upper().replace('T', 'U')
        self.mir_bulge = self.mir_bulge.upper().replace('T', 'U')

    def __str__(self):
        return self.tostring()


