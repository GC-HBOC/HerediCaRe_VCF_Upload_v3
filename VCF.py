import pandas as pd
import os
from datetime import datetime

class VCF:

    def __init__(self, HG38_FLAG, VERSION):
        self.hg38 = HG38_FLAG


        self.VERSION = VERSION
        ### infos for VCF_UPLOAD_META
        self.VCF_NAME = None
        self.MEMBER_ID = None # interne Pat-ID, Teil 2 des Dateinamens
        self.BOGEN_NR = None # MGU Bogennr, Teil 3 des Dateinamens
        self.ERFMIT = None # interne Mitarbeite-ID, Teil 4 des Dateinamens
        self.ERFDAT = None # Zeitstempel des Uploads, Teil 5 des Dateinames
        self.ERROR = []
        self.ERROR_LONG = []
        
        self.PARSE_N_IN_SOURCE = None

        self.header = [] # list of header lines 
        # lid = line ID, is required to calculate PARSE_ROWS_PROCESSED in final output
        VARIANT_HEADER = ['chrom', 'pos_hg38', 'ref_hg38', 'alt_hg38', 'pos_hg19', 'ref_hg19', 'alt_hg19', 'gene',  'transcript', 'hgvsc', 'hgvsp', 'effect', 'annotation', 'class', 'gt', 'norm_fail', 'ref_fail', 'liftover_fail', 'lid']
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
                    self.variants.loc[i,'pos_hg38'] = str(POS) #XXX
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
                    self.variants.loc[i,'pos_hg19'] = str(POS) # XXX
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
        
        #SRC_DIR = os.path.dirname(os.path.realpath(__file__))
        SRC_DIR = os.getcwd()
        CHAIN_FILE = SRC_DIR + '\\' + 'resources/hg38ToHg19.over.chain' if self.hg38 else SRC_DIR + '\\' + r'resources/hg19ToHg38.over.chain'
        
        for i in range(len(self.variants)):
            if self.hg38:
                CHROM, POS, REF, ALT = 'chr' + self.variants['chrom'][i], int(self.variants['pos_hg38'][i]), self.variants['ref_hg38'][i], self.variants['alt_hg38'][i]  
            else:
                CHROM, POS, REF, ALT = 'chr' + self.variants['chrom'][i], int(self.variants['pos_hg19'][i]), self.variants['ref_hg19'][i], self.variants['alt_hg19'][i] 
            #print(CHROM, POS, REF, ALT)
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
                                        self.variants.loc[i,'pos_hg19'] = str(QUERY_POS)
                                        self.variants.loc[i,'ref_hg19'] = REF
                                        self.variants.loc[i, 'alt_hg19'] = ALT
                                    else:
                                        self.variants.loc[i,'pos_hg38'] = str(QUERY_POS)
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
                                        self.variants.loc[i,'pos_hg19'] = str(QUERY_POS)
                                        self.variants.loc[i,'ref_hg19'] = REF
                                        self.variants.loc[i, 'alt_hg19'] = ALT
                                    else:
                                        self.variants.loc[i,'pos_hg38'] = str(QUERY_POS)
                                        self.variants.loc[i,'ref_hg38'] = REF
                                        self.variants.loc[i, 'alt_hg38'] = ALT
                                    self.variants.loc[i,'liftover_fail'] = False
                                CHAIN_FLAG = False
                            else:
                                REF_POS += dref
                                QUERY_POS += dquery
                                
                            
                            if CHAIN_FLAG and REF_POS >= POS:
                                CHAIN_FLAG =False


    def write_sql_output(self, outpath):


        ###  write meta output
        # check if (and how many) variants could be processed
        if self.variants['ref_fail'].any():
            PARSE_RESULT = '0'
            print('PARSE_RESULT ref_fail') #DEBUG
        elif len(self.variants.loc[(self.variants['norm_fail'] == False) & (self.variants['liftover_fail'] == False) & (self.variants['gene'].notna()) ]):
            PARSE_RESULT = '1'
        else:
            PARSE_RESULT = '0'

        db_entries_meta = '(UPLDATEI,MEMBER_ID,EVENT,ERFMIT,ERFDAT,REFGEN_UPLOAD,PARSE_DATE,PARSE_VERSION,PARSE_RESULT,PARSE_ROWS_IN_SOURCE,PARSE_ROWS_PROCESSED,PARSE_VARS_PROCESSED,ERROR_SHORT,ERROR_TEXT)'
        # In der bereits vorhandenen Tabelle VCF_UPLOAD haben wir folgende 2 Felder ergänzt:
            #PARSE_SUCCESFULL als einstelliges Zahlenfeld, zu befüllen mit 1 wenn die Verarbeitung erfolgreich war sonst 0.
            #ERROR_MSG als Textfeld (max 4000 Zeichen), zu Befüllen mit Fehlermeldung wenn die Verarbeitung nicht erfolgreich war
        db_entries_var = '(MEMBER_ID,BOGEN_NR,ERFMIT,ERFDAT, GEN2,REFSEQ,HGVS_DNA,HGVS_PROT,ART,PATH,CHROM,POS_HG19,REF_HG19,ALT_HG19,POS_HG38,REF_HG38,ALT_HG38,ZYGOT,PARSE_SUCCESSFUL,ERROR_MSG)'
        var_line_prefix =  "into VCF_UPLOAD " + db_entries_var + " values (" + ','.join([self.MEMBER_ID,self.BOGEN_NR,self.ERFMIT]) + ','
        
        with open(outpath, 'w') as outfile:
            outfile.write("Insert Into VCF_UPLOAD_META " + db_entries_meta + " Values ('")
            REFGEN = 'hg38' if self.hg38 else 'hg19'
            
            # to_date('13490910134119', 'dd.mm.yyyy hh24:mi:ss') => TO_DATE( '10.09.1349 13:41:19', 'dd.mm.yyyy hh24:mi:ss') # email from August 25, 26
            ERFDAT_str = '.'.join([self.ERFDAT[6:8], self.ERFDAT[4:6], self.ERFDAT[:4]]) + ' ' + ':'.join([ self.ERFDAT[8:10], self.ERFDAT[10:12], self.ERFDAT[12:14] ])

            #outfile.write(', '. join([self.VCF_NAME + "'", self.MEMBER_ID, self.BOGEN_NR, self.ERFMIT, "to_date('" + self.ERFDAT + "', 'dd.mm.yyyy hh24:mi:ss')", REFGEN, "to_date('" +  datetime.now().strftime("%Y%m%d%H%M%S") + "', 'dd.mm.yyyy hh24:mi:ss')" , "'"+ self.VERSION + "'"] ) + ', ')
            outfile.write(', '. join([self.VCF_NAME + "'", self.MEMBER_ID, self.BOGEN_NR, self.ERFMIT, "to_date('" + ERFDAT_str + "', 'dd.mm.yyyy hh24:mi:ss')", REFGEN, "to_date('" +  datetime.now().strftime("%d.%m.%Y %H%M%S") + "', 'dd.mm.yyyy hh24:mi:ss')" , "'"+ self.VERSION + "'"] ) + ', ')
            NVAR = len(self.variants.loc[(self.variants['norm_fail'] == False) & (self.variants['liftover_fail'] == False) & (self.variants['gene'].notna())]) 
            if not len(self.ERROR):
                outfile.write(', '. join([PARSE_RESULT, str(self.PARSE_N_IN_SOURCE), str(len(self.variants['lid'].unique())), str(NVAR), "'###'", "'###'" ]))
            else:
            
                outfile.write(', '. join([PARSE_RESULT, str(self.PARSE_N_IN_SOURCE), str(len(self.variants['lid'].unique())), str(NVAR), "'" + ';'.join(self.ERROR) + "'", "'" + ';'.join(self.ERROR_LONG) + "'" ]))
            outfile.write(")\n\n")



        ### write variant output
        #db_entries = '(UPLDATEI,MEMBER_ID,BOGEN_NR,ERFMIT,ERFDAT, GEN2,REFSEQ,HGVS_DNA,HGVS_PROT,ART,PATH,CHROM,POS_HG19,REF_HG19,ALT_HG19,POS_HG38,REF_HG38,ALT_HG38,ZYGOT,PARSE_N_IN_SOURCE,PARSE_N_PROCESSED,ERROR_SHORT,ERROR_TEXT)'
        # LINE_PREFIX = "into VCF_UPLOAD " + db_entries_var + " values"
        
            print('PARSE_RESULT', PARSE_RESULT)

            if PARSE_RESULT == '1':
                outfile.write("Insert all\n")


                for i in range(len(self.variants)):
                    tmp_var = self.variants.loc[i,:]
                    # '(MEMBER_ID,BOGEN_NR,ERFMIT,ERFDAT, GEN2,REFSEQ,HGVS_DNA,HGVS_PROT,ART,PATH,CHROM,POS_HG19,REF_HG19,ALT_HG19,POS_HG38,REF_HG38,ALT_HG38,ZYGOT,PARSE_SUCCESSFUL,ERROR_MSG)'
                    if tmp_var['norm_fail']:
                        pass
                    elif tmp_var['liftover_fail']:
                        pass
                    else:
                        VAR_ENTRIES = ["'" + _ + "'" for _ in [tmp_var['gene'],tmp_var['transcript'], tmp_var['hgvsc'] ]]
                        if tmp_var['hgvsp']: VAR_ENTRIES.append("'" + tmp_var['hgvsp'] +"'")
                        else: VAR_ENTRIES.append("'###'")
                        VAR_ENTRIES.append("'" + tmp_var['effect'] +"'")

                        annerr = ''
                        if tmp_var['annotation']:
                            # TODO set WARNING/ERROR (annerr) in MSG
                            annot, annerr = VCF.return_path_class(str(tmp_var['annotation']), str(tmp_var['class']))
                            print(str(tmp_var['annotation']), str(tmp_var['class']), ':', annot, annerr)
                            if annot and not annerr:
                                VAR_ENTRIES.append("'" + str(annot) + "'")
                            else:
                                 VAR_ENTRIES.append("'###'")
                        else:
                            VAR_ENTRIES.append("'###'")
                        VAR_ENTRIES.append(str(tmp_var['chrom']) )
                        ### XXX should not happen that POS, REF, ALT unavailable
                        if tmp_var['pos_hg19']:
                            VAR_ENTRIES.append(str(tmp_var['pos_hg19']) )
                        else: VAR_ENTRIES.append("'###'")
                        if tmp_var['ref_hg19']:
                            VAR_ENTRIES.append("'" + str(tmp_var['ref_hg19']) + "'" )
                        else: VAR_ENTRIES.append("'###'")
                        if tmp_var['alt_hg19']:
                            VAR_ENTRIES.append("'" + str(tmp_var['alt_hg19']) + "'")
                        else: VAR_ENTRIES.append("'###'")
                        if tmp_var['pos_hg38']:
                            VAR_ENTRIES.append(str(tmp_var['pos_hg38']) )
                        else: VAR_ENTRIES.append("'###'")
                        if tmp_var['ref_hg38']:
                            VAR_ENTRIES.append("'" + str(tmp_var['ref_hg38'])  + "'")
                        else: VAR_ENTRIES.append("'###'")
                        if tmp_var['alt_hg38']:
                            VAR_ENTRIES.append("'" + str(tmp_var['alt_hg38'])  + "'")
                        else: VAR_ENTRIES.append("'###'")
                        
                        gterr = ''
                        if tmp_var['gt']:
                            if tmp_var['gt'] == 1:
                                VAR_ENTRIES.append("'" + '0/1' + "'")
                            elif tmp_var['gt'] == 2:
                                VAR_ENTRIES.append("'" + '1/1' + "'")
                            else:
                                gterr = 'GT=' + str(tmp_var['gt']) + '?'
                        else:
                            VAR_ENTRIES.append("'###'")
                        


                        if annerr or gterr:
                            VAR_ENTRIES.append(str(0))
                            VAR_ENTRIES.append("'" + ';'.join([_ for _ in [annerr, gterr] if _ != '']) + "'")
                        else:
                            VAR_ENTRIES.append(str(1))
                            VAR_ENTRIES.append("'###'")

                    outfile.write(var_line_prefix + ','.join(VAR_ENTRIES))


                    outfile.write(')\n')


                outfile.write("select * from dual;\n")
                outfile.write("commit;\n")

    @staticmethod
    def return_path_class(annot, pathclass):
        if annot == "MT":
            if pathclass == "BEN":
                return (1, '')
            elif pathclass == "LBEN":
                return (2, '')
            elif pathclass == "UI":
                return (3,'')
            elif pathclass == "LPAT":
                return (4,'')
            elif pathclass == "PAT":
                return (5,'')
            elif pathclass == "UD":
                return ('','') # keine Angabe
            else:
                return ('', 'Unknown tag for MT annotation of pathogenicity: ' + pathclass)
        
        elif annot == "MutDB:Classification":
            if pathclass == "benign":
                return (1, '')
            elif pathclass == "likely benign":
                return (2, '')
            elif pathclass == "uncertain significance":
                return (3,'')
            elif pathclass == "likely pathogenic":
                return (4,'')
            elif pathclass == "pathogenic":
                return (5,'')
            elif pathclass == "undefined":
                return ('','') # keine Angabe
            else:
                return ('', 'Unknown tag for MT annotation of pathogenicity: ' + pathclass)

        elif annot == "CLASS":
            if pathclass == "1":
                return (1, '')
            elif pathclass == "1":
                return (2, '')
            elif pathclass == "3":
                return (3,'')
            elif pathclass == "4":
                return (4,'')
            elif pathclass == "5":
                return (5,'')
            elif pathclass == "":
                return ('','') # keine Angabe
            else:
                return ('', 'Unknown tag for CLASS annotation of pathogenicity: ' + pathclass)

        else:
            return('', 'Unknown tag for pathogenicity annotation: ' + pathclass )

            




"""             ## REF check failed
            if self.variants['ref_fail'].any():
                _ind = self.variants.loc[self.variants['ref_fail'] == True].index[0]
                if self.hg38:
                    _var = '-'.join([self.variants.loc[_ind,'chrom'], self.variants.loc[_ind,'pos_hg38'],self.variants.loc[_ind,'ref_hg38'], self.variants.loc[_ind,'alt_hg38'] ])
                else:
                    _var = '-'.join([self.variants.loc[_ind,'chrom'], self.variants.loc[_ind,'pos_hg19'],self.variants.loc[_ind,'ref_hg19'], self.variants.loc[_ind,'alt_hg19'] ])

                ENTRY_LIST = [self.VCF_NAME,self.MEMBER_ID, self.BOGEN_NR, self.ERFMIT, self.ERFDAT] # UPLDATEI,MEMBER_ID,BOGEN_NR,ERFMIT,ERFDAT,
                ENTRY_LIST = ENTRY_LIST + ['###','###','###','###','###','###','###','###','###','###','###','###','###','###'] # GEN2,REFSEQ,HGVS_DNA,HGVS_PROT,ART,PATH,CHROM,POS_HG19,REF_HG19,ALT_HG19,POS_HG38,REF_HG38,ALT_HG38,ZYGOT
                ENTRY_LIST = ENTRY_LIST + [self.PARSE_N_IN_SOURCE,0,"\'REF_CHECK_FAIL\'", "\'Reference check failed for variant " + _var + "\'" ] # PARSE_N_IN_SOURCE,PARSE_N_PROCESSED,ERROR_SHORT,ERROR_TEXT
                OUT = "into VCF_UPLOAD " + db_entries + " values ("

                OUT += ','.join([str(_) for _ in ENTRY_LIST]) + ')\n'
                outfile.write(OUT)
            
            ## there are variants to report!!
            elif len(self.variants.loc[(self.variants['norm_fail'] == False) & (self.variants['liftover_fail'] == False) & (self.variants['gene'].notna()) ]):
                pass

            else:
                pass
#                if self.hg38:
#                    _var = '-'.join([self.variants.loc[_ind,'chrom'], self.variants.loc[_ind,'pos_hg38'],self.variants.loc[_ind,'ref_hg38'], self.variants.loc[_ind,'alt_hg38'] ])
#                else:
#                    _var = '-'.join([self.variants.loc[_ind,'chrom'], self.variants.loc[_ind,'pos_hg19'],self.variants.loc[_ind,'ref_hg19'], self.variants.loc[_ind,'alt_hg19'] ])
#                OUT = "into VCF_UPLOAD (MEMBER_ID,BOGEN_NR,ERFMIT,ERFDAT,GEN2,HGVS_DNA) values ("
#                OUT += ','.join([self.MEMBER_ID, self.BOGEN_NR, self.ERFMIT, self.ERFDAT, "\'STATUS\'", "\'No valid variants to report\'"]) + ')\n'
#                outfile.write(OUT)




            outfile.write("select * from dual;\n")
            outfile.write("commit;\n") """

"""     def write_sql_meta_output(self, outpath):

        # check if (and how many) variants could be processed
        if self.variants['ref_fail'].any():
            PARSE_RESULT = '0'
        elif len(self.variants.loc[(self.variants['norm_fail'] == False) & (self.variants['liftover_fail'] == False) & (self.variants['gene'].notna()) ]):
            PARSE_RESULT = '1'
        else:
            PARSE_RESULT = '0'

        db_entries = '(UPLDATEI,MEMBER_ID,EVENT,ERFMIT,ERFDAT,REFGEN_UPLOAD,PARSE_DATE,PARSE_VERSION,PARSE_RESULT,PARSE_ROWS_IN_SOURCE,PARSE_ROWS_PROCESSED,PARSE_VARS_PROCESSED,ERROR_SHORT,ERROR_TEXT)'
        with open(outpath, 'w') as outfile:
            outfile.write("Insert Into VCF_UPLOAD_META " + db_entries + " Values ('")
            REFGEN = 'hg38' if self.hg38 else 'hg19'
            
            # to_date('13490910134119', 'dd.mm.yyyy hh24:mi:ss') => TO_DATE( '10.09.1349 13:41:19', 'dd.mm.yyyy hh24:mi:ss') # email from August 25, 26
            ERFDAT_str = '.'.join([self.ERFDAT[6:8], self.ERFDAT[4:6], self.ERFDAT[:4]]) + ' ' + ':'.join([ self.ERFDAT[8:10], self.ERFDAT[10:12], self.ERFDAT[12:14] ])

            #outfile.write(', '. join([self.VCF_NAME + "'", self.MEMBER_ID, self.BOGEN_NR, self.ERFMIT, "to_date('" + self.ERFDAT + "', 'dd.mm.yyyy hh24:mi:ss')", REFGEN, "to_date('" +  datetime.now().strftime("%Y%m%d%H%M%S") + "', 'dd.mm.yyyy hh24:mi:ss')" , "'"+ self.VERSION + "'"] ) + ', ')
            outfile.write(', '. join([self.VCF_NAME + "'", self.MEMBER_ID, self.BOGEN_NR, self.ERFMIT, "to_date('" + ERFDAT_str + "', 'dd.mm.yyyy hh24:mi:ss')", REFGEN, "to_date('" +  datetime.now().strftime("%d.%m.%Y %H%M%S") + "', 'dd.mm.yyyy hh24:mi:ss')" , "'"+ self.VERSION + "'"] ) + ', ')
            NVAR = len(self.variants.loc[(self.variants['norm_fail'] == False) & (self.variants['liftover_fail'] == False) & (self.variants['gene'].notna())]) 
            if not len(self.ERROR):
                outfile.write(', '. join([PARSE_RESULT, str(self.PARSE_N_IN_SOURCE), str(len(self.variants['lid'].unique())), str(NVAR), "'###'", "'###'" ]))
            else:
            
                outfile.write(', '. join([PARSE_RESULT, str(self.PARSE_N_IN_SOURCE), str(len(self.variants['lid'].unique())), str(NVAR), "'" + ';'.join(self.ERROR) + "'", "'" + ';'.join(self.ERROR_LONG) + "'" ]))
            outfile.write(")\n") """


