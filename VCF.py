import pandas as pd
import os

class VCF:
    def __init__(self, HG38_FLAG):
        self.hg38 = HG38_FLAG
        
        ### further infos required for SQL output
        self.MEMBER_ID = None # interne Pat-ID, Teil 2 des Dateinamens
        self.BOGEN_NR = None # MGU Bogennr, Teil 3 des Dateinamens
        self.ERFMIT = None # interne Mitarbeite-ID, Teil 4 des Dateinamens
        self.ERFDAT = None # Zeitstempel des Uploads, Teil 5 des Dateinames
        
        self.processed = False
        self.header = [] # list of header lines 
        VARIANT_HEADER = ['chrom', 'pos_hg38', 'ref_hg38', 'alt_hg38', 'pos_hg19', 'ref_hg19', 'alt_hg19', 'gene',  'transcript', 'hgvsc', 'hgvsp', 'effect', 'annotation', 'class', 'gt', 'norm_fail', 'ref_fail', 'liftover_fail']
        self.variants = pd.DataFrame(columns=VARIANT_HEADER)
    
    def normalize(self, seq_dict):
        # https://genome.sph.umich.edu/wiki/Variant_Normalization
        if self.hg38:
            for i in range(len(self.variants)):
                try:
                    CHROM, POS, REF, ALT = self.variants['chrom'][i], int(self.variants['pos_hg38'][i]), self.variants['ref_hg38'][i], self.variants['alt_hg38'][i]
                    change_flag = True
                    while change_flag:
                        change_flag = False
                        if REF[-1] == ALT[-1]:
                            REF, ALT = REF[:-1], ALT[:-1]
                            change_flag = True
                        if not len(REF) or not len(ALT):
                            POS = POS -1
                            REF = seq_dict[CHROM][POS-1].upper() + REF
                            ALT = seq_dict[CHROM][POS-1].upper() + ALT
                            change_flag = True
                    while REF[0] == ALT[0] and len(REF) > 1 and len(ALT) > 1:
                        REF, ALT = REF[1:], ALT[1:]
                        POS += 1
                    self.variants.loc[i,'pos_hg38'] = str(POS)
                    self.variants.loc[i,'ref_hg38'] = REF
                    self.variants.loc[i,'alt_hg38'] = ALT
                    self.variants.loc[i, 'norm_fail'] = False
                except:
                    self.variants.loc[i, 'norm_fail'] = True
            
        else:
            for i in range(len(self.variants)):
                try:
                    CHROM, POS, REF, ALT = self.variants['chrom'][i], int(self.variants['pos_hg19'][i]), self.variants['ref_hg19'][i], self.variants['alt_hg19'][i]
                    change_flag = True
                    while change_flag:
                        change_flag = False
                        if REF[-1] == ALT[-1]:
                            REF, ALT = REF[:-1], ALT[:-1]
                            change_flag = True
                        if not len(REF) or not len(ALT):
                            POS = POS -1
                            REF = seq_dict[CHROM][POS-1].upper() + REF
                            ALT = seq_dict[CHROM][POS-1].upper() + ALT
                            change_flag = True
                    while REF[0] == ALT[0] and len(REF) > 1 and len(ALT) > 1:
                        REF, ALT = REF[1:], ALT[1:]
                        POS += 1
                    self.variants.loc[i,'pos_hg19'] = str(POS)
                    self.variants.loc[i,'ref_hg19'] = REF
                    self.variants.loc[i,'alt_hg19'] = ALT
                    self.variants.loc[i, 'norm_fail'] = False
                except:
                    self.variants.loc[i, 'norm_fail'] = True
    
    def liftover(self, seq_dict):
        # https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
        # http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz
        # https://genome.ucsc.edu/goldenPath/help/chain.html
        # DEBUG: https://liftover.broadinstitute.org/
        # TODO: reverse query sequence in chain specification
        
        SRC_DIR = os.path.dirname(os.path.realpath(__file__))
        CHAIN_FILE = SRC_DIR + '\\' + 'resources/hg38ToHg19.over.chain' if self.hg38 else SRC_DIR + '\\' + r'resources\hg19ToHg38.over.chain'
        
        for i in range(len(self.variants)):
            if self.hg38:
                CHROM, POS, REF, ALT = 'chr' + self.variants['chrom'][i], int(self.variants['pos_hg38'][i]), self.variants['ref_hg38'][i], self.variants['alt_hg38'][i]  
            else:
                CHROM, POS, REF, ALT = 'chr' + self.variants['chrom'][i], int(self.variants['pos_hg19'][i]), self.variants['ref_hg19'][i], self.variants['alt_hg19'][i] 
            print(CHROM, POS, REF, ALT)
            CHAIN_FLAG = False
            self.variants.loc[i,'liftover_fail'] = True
            #print(self.variants)
            with open(CHAIN_FILE) as infile:
                for line in infile:
                    if line.startswith('chain'):
                        ll = line.rstrip().split()
                        if ll[2] == CHROM and int(ll[5]) <= POS-1 < int(ll[6]):
                            CHAIN_FLAG = True
                            REF_POS = int(ll[5])
                            QUERY_POS = int(ll[10])
                            if ll[9] == '-':
                                # TODO: reverse query sequence in chain specification
                                self.variants.loc[i,'liftover_fail'] = True
                                CHAIN_FLAG = False
                    
                    elif CHAIN_FLAG:
                        ll = line.rstrip().split()
                        
                        if len(ll) == 1:
                            # treat last line in block
                            REF_POS += int(ll[0])
                            QUERY_POS += int(ll[0])
                            if REF_POS >= POS:
                                QUERY_POS = QUERY_POS - (REF_POS - POS)
                                #print('QUERY_POS final (last block)', QUERY_POS)
                                QUERY_SEQ = seq_dict[CHROM[3:]][QUERY_POS-1:(QUERY_POS + len(REF) -1)].upper()
                                #print('QUERY_SEQ', QUERY_SEQ)
                                if QUERY_SEQ == REF.upper():
                                    if self.hg38:
                                        self.variants.loc[i,'pos_hg19'] = QUERY_POS
                                        self.variants.loc[i,'ref_hg19'] = REF
                                        self.variants.loc[i, 'alt_hg19'] = ALT
                                    else:
                                        self.variants.loc[i,'pos_hg38'] = QUERY_POS
                                        self.variants.loc[i,'ref_hg38'] = REF
                                        self.variants.loc[i, 'alt_hg38'] = ALT
                                    self.variants.loc[i,'liftover_fail'] = False
                            CHAIN_FLAG = False
                        
                        else:
                            block = int(ll[0])
                            dref, dquery = int(ll[1]), int(ll[2])
                            REF_POS += block
                            QUERY_POS += block
                            if REF_POS >= POS:
                                QUERY_POS = QUERY_POS - (REF_POS - POS)
                                REF_POS = REF_POS - (REF_POS -POS)
                                #print('QUERY_POS final', QUERY_POS)
                                QUERY_SEQ = seq_dict[CHROM[3:]][QUERY_POS-1:(QUERY_POS + len(REF) -1)].upper()
                                #print('QUERY_SEQ', QUERY_SEQ)
                                if QUERY_SEQ == REF.upper():
                                    if self.hg38:
                                        self.variants.loc[i,'pos_hg19'] = QUERY_POS
                                        self.variants.loc[i,'ref_hg19'] = REF
                                        self.variants.loc[i, 'alt_hg19'] = ALT
                                    else:
                                        self.variants.loc[i,'pos_hg38'] = QUERY_POS
                                        self.variants.loc[i,'ref_hg38'] = REF
                                        self.variants.loc[i, 'alt_hg38'] = ALT
                                    self.variants.loc[i,'liftover_fail'] = False
                                CHAIN_FLAG = False
                            else:
                                REF_POS += dref
                                QUERY_POS += dquery
                                
                            
                            if CHAIN_FLAG and REF_POS >= POS:
                                CHAIN_FLAG =False

