import argparse
import os # required to get *.py and resources folder
import shutil
import sys
import chardet
#import pysam
import subprocess
from VCF import VCF
import re
import fastapy


VERSION = "1.0.0"

SRC_DIR = os.path.dirname(os.path.realpath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument("input_folder", help="Path to folder with VCF files to parse")
parser.add_argument("-o", "--output_folder", default=SRC_DIR + '\\' +"Output", help="Output  folder  for  final  .txt  files  (default: Output)")
parser.add_argument("-sp", default=SRC_DIR + '\\' + r"resources\snpEff\snpEff.jar", help=r"Path to snpEff .jar (default: resources\snpEff\snpEff.jar)")
parser.add_argument("-jp", default='java', help="Path to java executable for running snpEff (default: java)") # "O:\microsoft-jdk-25.0.2-windows-x64\jdk-25.0.2+10\bin\java.exe" (MHH) or "V:\Bioinformatik\software\jdk-25.0.1.8-hotspot\bin\java.exe" (FBZ)
parser.add_argument("-d", "--debug_folder", default="Debug", help="Debug Folder. Contains processed, erroneous & normalized VCFs + Rejected Variants TSV file (default: Debug)")
parser.add_argument("-ram", default=4, type=int, help="Accessible  RAM  (GB)  for  java virtual machine (default: 4)")
parser.add_argument("-t", "--transcript_tsv", default= SRC_DIR + '\\' + r'resources\transcripts.tsv', help="Transcripts TSV file (default: resources\\transcripts.tsv)")
args = parser.parse_args()

##
sys.stderr.write('... reading transcripts file ' + args.transcript_tsv + '.\n')
TRANSCRIPTS = dict()
with open(args.transcript_tsv) as infile:
    c = 1
    for line in infile:
        #ABRAXAS1	NM_139076	1
        #ACD	NM_001082486	1
        ll = line.rstrip().split('\t')
        if len(ll) >= 3 and ll[1].startswith('NM_') and ll[2] in ['0', '1']:
            TRANSCRIPTS[ll[1].split('.')[0]] = (ll[0], int(ll[2])) 
        else:
            sys.stderr.write('Could not parse line ' + str(c) + ': ' + line.rstrip() + '\n' )
        c+=1

sys.stderr.write('... reading reference FASTA files.\n')
hg38records = fastapy.parse(os.path.dirname(os.path.realpath(__file__)) + '\\' + r'resources\ref\GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz')
#print(hg38records)
hg38_dict = dict()
for record in hg38records:
    if record.id in list(['chr' + str(_) for _ in range(1,23)]) + ['chrX', 'chrY']:
        hg38_dict[record.id[3:]] = record.seq

hg19records = fastapy.parse(os.path.dirname(os.path.realpath(__file__)) + '\\' + r'resources\ref\hs37d5.fa.gz')
hg19_dict = dict()
for record in hg19records:
    #print(record.id)
    if record.id in list([str(_) for _ in range(1,23)]) + ['X', 'Y']:
        hg19_dict[record.id] = record.seq
sys.stderr.write('Reference FASTA files parsed.\n')

try:
    VCFS = [_ for _ in os.listdir(args.input_folder)]
    print(VCFS)
except:
    sys.exit("Can not access input folder " + args.input_folder + " ...TERMINATING\n")


fname_pattern = r"^hg\d\d-\S+-\S+-\S+-\S+.vcf"
for VCF_FILE in VCFS:
    if not VCF_FILE.startswith('hg19-') and not VCF_FILE.startswith('hg38-'):
        sys.stderr.write('Can not process input file ' + VCF_FILE + '... prefix not hg19 or hg38\n')
    elif not VCF_FILE.endswith('.vcf'):
        sys.stderr.write('Can not process input file ' + VCF_FILE + '... file name does not end with .vcf\n')
        ### TODO Check name pattern <hg19|hg38>-<interne PatientenID>-<MGU-Bogennr>-<interne Mitarbeiter-ID>-<Zeitstempel>.vcf
        ### TODO check of timestamp, has to be YYYYMMDDHH24MISSS (??)
    elif not bool(re.match(fname_pattern, VCF_FILE)):
        sys.stderr.write('Can not process input file ' + VCF_FILE + '... file name does not match file name convention\n')
    else:

        rawdata = open(args.input_folder + '\\' + VCF_FILE, "rb").read()
        encoding = chardet.detect(rawdata)['encoding']
        print(VCF_FILE, encoding)
        #print(VCF)
        HG38_FLAG = False if VCF_FILE.startswith('hg19-') else True
        
        

        ## normalize VCF: ends up with a normalized hg38 VCF file in vcf_Normalisiert 
        if HG38_FLAG:
            vcf = VCF(True, VERSION=VERSION)
            #normalize_vcf(os.path.join(args.input_folder, VCF), hg38_dict, HG38_FLAG)
        if not HG38_FLAG:
            vcf = VCF(False, VERSION=VERSION)
            #normalize_vcf(os.path.join(args.input_folder, VCF), hg19_dict, HG38_FLAG)

        ### further infos required for SQL output
        #hg38-2155522-1-2026-20180910134119.vcf
        vcf.VCF_NAME = VCF_FILE
        vcf.MEMBER_ID = VCF_FILE.split('-')[1]
        vcf.BOGEN_NR = VCF_FILE.split('-')[2]
        vcf.ERFMIT =  VCF_FILE.split('-')[3]
        vcf.ERFDAT =  VCF_FILE.split('-')[4].split('.')[0]


        LINE_COUNTER, PROCESSED_COUNTER, LID = 0, 0, 0
        with open(os.path.join(args.input_folder, VCF_FILE), encoding=chardet.detect(rawdata)['encoding']) as infile:
            FAIL_FLAG = False
            for _l in infile:
                line = _l.rstrip().strip('"')
                if len(line):
                    if line.startswith('#'):
                        vcf.header.append(line)
                        #try:
                    else:
                        LINE_COUNTER +=1
                        ll = line.rstrip().split('\t')
                        if len(ll) not in [8, 10]:
                            sys.stderr.write('...invalid number of columns in VCF file ' + VCF_FILE + ': ' + str(len(ll)) + '\n')
                            FAIL_FLAG = True
                            vcf.ERROR.append('INVALID_NUMBER_OF_COLUMNS')
                            vcf.ERROR_LONG.append('INVALID_NUMBER_OF_COLUMNS: ' + l)
                            break
                        try:
                            CHROM, POS, REF, ALT, INFO = ll[0], ll[1], ll[3], ll[4], ll[7]
                        except:
                            sys.stderr.write('... can not parse: ' + line)
                        #Es besteht für Nutzerinnen und Nutzer die Möglichkeit, zusätzliche Informationen zur Klassifizierung der Pathogenität von Varianten in der INFO-Spalte 
                        #(Spalte 8) mithilfe der Schlagworte MutDB:Classification, CLASS oder MT zu hinterlegen Sind mehrere dieser Einträge für die selbe Variante vorhanden, 
                        # wird der MUtDB:Classification-Eintrag vor allen anderen und der CLASS-Eintrag vor dem MT-Eintrag priorisiert.
                        if "MutDB_Classification" in INFO:
                            ANNOT_TAG = "MutDB_Classification"
                        elif "CLASS" in INFO:
                            ANNOT_TAG = "CLASS"
                        elif "MT" in INFO:
                            ANNOT_TAG = "MT"
                        else:
                            ANNOT_TAG = None
                        # chek annotation for multi-ALT alleles 
                        nalt = len(ALT.split(','))
                        ANNOT = [_ for _ in INFO.split(';') if _.startswith(ANNOT_TAG + '=')] if ANNOT_TAG else []
                        if len(ANNOT):
                            ANNOT = ANNOT[0].split('=')[1]
                            if len(ANNOT.split(',')) != nalt:
                                sys.stderr.write('...unable to parse annotation for variant ' + '-'.join([CHROM,POS,REF,ALT]) + ' in VCF file ' + VCF_FILE + '\n')
                                FAIL_FLAG = True
                                vcf.ERROR.append('ANNOTATION_PARSE_ERROR')
                                vcf.ERROR_LONG.append('ANNOTAION_PARSE_ERROR for variant' + '-'.join([CHROM,POS,REF,ALT]) )
                        
                        
                        ## split variant in single-ALT variants
                        if not FAIL_FLAG:
                            # [1] Mitochondrial variants are ignored & prefix chr are removed (chr1 --> 1)
                            if CHROM.startswith('chr') or CHROM.startswith('Chr'): CHROM = CHROM[3:]
                            if CHROM == "23": CHROM = "X"
                            if CHROM == "24": CHROM = "Y"
                            if CHROM in list([str(_) for _ in range(1,23)]) + ['X', 'Y']:
                                for i in range(len(ALT.split(','))):
                                    _ALT = ALT.split(',')[i]
                                    if _ALT not in ['.', '*']:
                                        GT = ll[9].split(':')[0].count(str(i+1)) if len(ll)>8 else None
                                        varclass = ANNOT.split(',')[i] if len(ANNOT) else None
                                        # ['chrom', 'pos_hg38', 'ref_hg38', 'alt_hg38', 'pos_hg19', 'ref_hg19', 'alt_hg19', 'gene',  'transcript', 'hgvsc', 'hgvsp', 'effect', 'annotation', 'class', 'gt', 'norm_fail', 'ref_fail', 'liftover_fail']
                                        if HG38_FLAG:
                                            ## TODO REF check
                                            if REF.upper() == hg38_dict[CHROM][int(POS)-1:int(POS)+len(REF)-1].upper():
                                            
                                                VAR = [CHROM, POS, REF, _ALT, None, None, None, None, None, None, None, None, ANNOT_TAG, varclass, GT, None, False, None, LID]
                                            else:
                                                VAR = [CHROM, POS, REF, _ALT, None, None, None, None, None, None, None, None, ANNOT_TAG, varclass, GT, None, True, None, LID]
                                        else:
                                            if REF.upper() == hg19_dict[CHROM][int(POS)-1:int(POS)+len(REF)-1].upper():
                                                VAR = [CHROM, None, None, None, POS, REF, _ALT, None, None, None, None, None, ANNOT_TAG, varclass, GT, None, False, None, LID]
                                            else:
                                                VAR = [CHROM, None, None, None, POS, REF, _ALT, None, None, None, None, None, ANNOT_TAG, varclass, GT, None, True, None, LID]
                                        vcf.variants.loc[len(vcf.variants)] = VAR
                                    else:
                                        sys.stderr.write("Variant with ALT " + _ALT + " at " + CHROM + ':' + str(POS) + ' is ignored\n')
                        LID +=1

        print(vcf.variants)
        print('FAIL_FLAG:', FAIL_FLAG)
        vcf.PARSE_N_IN_SOURCE =  LINE_COUNTER
        
        # set FAIL_FLAG if any ref_fail == True
        if vcf.variants['ref_fail'].any(): 
            FAIL_FLAG = True
            vcf.ERROR.append('REF_FAIL_ERROR')
            tmp = vcf.variants.loc[vcf.variants['ref_fail']==True]
            if not vcf.hg38:
                vcf.ERROR_LONG.append('REF_FAIL_ERROR for variant ' + '-'.join([str(tmp['chrom'][0]), str(tmp['pos_hg19'][0]) , tmp['ref_hg19'][0], tmp['alt_hg19'][0]]) )
            else:
                vcf.ERROR_LONG.append('REF_FAIL_ERROR for variant ' + '-'.join([str(tmp['chrom'][0]), str(tmp['pos_hg38'][0]) , tmp['ref_hg38'][0], tmp['alt_hg38'][0]]) )
            
        ### NORMALIZATION
        if HG38_FLAG:
            vcf.normalize(hg38_dict)
        else:
            vcf.normalize(hg19_dict)
        
        ### LIFTOVER
        if HG38_FLAG:
            vcf.liftover(hg19_dict)
        else:
            vcf.liftover(hg38_dict)
        #print(vcf.variants)
            
        if FAIL_FLAG:
            # Create the directory if it doesn't exist
            os.makedirs('rejected_vcf_input', exist_ok=True)
            shutil.move(os.path.join(args.input_folder, VCF_FILE), os.path.join('rejected_vcf_input', VCF_FILE))
            
            sys.stderr.write('...unable to parse ' + VCF_FILE + '\n')
        
        else:
            print(vcf.variants)
            
            # write temporary snpEff input file
            os.makedirs('tmp', exist_ok=True)
            with open('tmp/' + VCF_FILE + '.tmp.vcf', 'w') as outfile:
                ##fileformat=VCFv4.2
                #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
                outfile.write('##fileformat=VCFv4.2\n')
                outfile.write('\t'.join(['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']) )
                for i in range(len(vcf.variants)):
                    if vcf.variants.loc[i,'norm_fail'] == vcf.variants.loc[i,'ref_fail'] == False and (HG38_FLAG or vcf.variants.loc[i,'liftover_fail'] == False):
                        CHROM, POS, REF, ALT = 'chr' + vcf.variants.loc[i,'chrom'], str(vcf.variants.loc[i,'pos_hg38']), vcf.variants.loc[i,'ref_hg38'], vcf.variants.loc[i,'alt_hg38']
                        outfile.write('\n'+ '\t'.join([CHROM, POS, '.', REF, ALT, '.', '.', '.']) )
            # run snpEff
            CMD = args.jp + ' -Xmx' + str(args.ram) + 'g -jar ' + args.sp + ' hg38 ' + 'tmp\\'  + VCF_FILE + '.tmp.vcf'
            #print(CMD)
            output = None
            try:
                output = subprocess.check_output(CMD, shell=True, text=True)
                sys.stderr.write('### running snpEff for', VCF_FILE)
                ## XXX DEBUG
                #subprocess.Popen(CMD + ' > snpeff_check\\' + VCF_FILE, shell=True ).wait()
                #print(output)
            except:
                sys.stderr.write('... Could not run ' + CMD + '\n')
            
            if output:
                #DEL_INDS = [] # store indices of variants not located within pre-defined transcripts
                IND_DICT = dict() # use dict to deal with variants located in different genes or transcripts 
                for line in output.split('\n'):
                    if not line.startswith('#') and len(line.split('\t')) >= 8:
                        ll = line.split('\t')
                        CHROM, POS, REF, ALT, INFO = ll[0][3:], ll[1], ll[3], ll[4], ll[7]
                        _inds = vcf.variants.index[(vcf.variants['chrom'] == CHROM) & (vcf.variants['pos_hg38'] == POS) & (vcf.variants['ref_hg38'] == REF) & (vcf.variants['alt_hg38'] == ALT)].tolist()


                        if len(_inds) == 1:
                            ind = _inds[0]
                            ### NOTE: due to self-generated VCF input, ANN is the only entry in INFO column 
                            print(line)
                            ANN = [_ for _ in INFO.split(';') if _.startswith('ANN=')][0]
                            ANN = [_ for _ in ANN[4:].split(',') if (_.split('|')[6].split('.')[0] in TRANSCRIPTS.keys())]
                            #print(ANN)
                            for ann in ANN:
                                _gene, _transcript = ann.split('|')[3], ann.split('|')[6]
                                _hgvsc, _hgvsp, _eff = ann.split('|')[9], ann.split('|')[10], ann.split('|')[1]
                                if ind not in IND_DICT.keys(): IND_DICT[ind] = dict()
                                IND_DICT[ind][_gene] = (_transcript, _hgvsc, _hgvsp, _eff)
                            #print(IND_DICT)

                        else:
                            ## deal with variant found less or more than once in vcf.variants
                            if not len(_inds):
                                sys.stderr.write("Couldn't identify variant " + '-'.join([CHROM, str(POS), REF, ALT]) + " from snpEff Output\n")
                            elif len(_inds) > 1:
                                sys.stderr.write("Several entries of variant " + '-'.join([CHROM, str(POS), REF, ALT]) + " in normalized VCF input ... Ignoring this variant.\n")
                                
                            #TODO variant not found or doubled            
                
                print(IND_DICT)

                ### treat variants located in different genes or transcripts
                _N = len(vcf.variants)
                for i in range(_N):
                    if i in IND_DICT:
                        ## ...otherwise variant is not in valid transcript from TRANSCRIPTS
                        K = list(IND_DICT[i].keys())
                        if len(K) == 1:
                            vcf.variants.loc[i,'gene'] = K[0]
                            vcf.variants.loc[i,'transcript'] = IND_DICT[i][K[0]][0]
                            vcf.variants.loc[i,'hgvsc'] = IND_DICT[i][K[0]][1]
                            if IND_DICT[i][K[0]][2]: vcf.variants.loc[i,'hgvsp'] = IND_DICT[i][K[0]][2]
                            vcf.variants.loc[i,'effect'] = IND_DICT[i][K[0]][3]

                        if len(K) > 1:
                            vcf.variants.loc[i,'gene'] = K[0]
                            vcf.variants.loc[i,'transcript'] = IND_DICT[i][K[0]][0]
                            vcf.variants.loc[i,'hgvsc'] = IND_DICT[i][K[0]][1]
                            if IND_DICT[i][K[0]][2]: vcf.variants.loc[i,'hgvsp'] = IND_DICT[i][K[0]][2]
                            vcf.variants.loc[i,'effect'] = IND_DICT[i][K[0]][3]
                            for g in K[1:]:
                                _i = len(vcf.variants)
                                row_to_copy = vcf.variants.loc[i].copy()
                                vcf.variants.loc[_i] = row_to_copy
                                vcf.variants.loc[_i,'gene'] = g
                                vcf.variants.loc[_i,'transcript'] = IND_DICT[i][g][0]
                                vcf.variants.loc[_i,'hgvsc'] = IND_DICT[i][g][1]
                                if IND_DICT[i][g][2]: vcf.variants.loc[_i,'hgvsp'] = IND_DICT[i][g][2]
                                vcf.variants.loc[_i,'effect'] = IND_DICT[i][g][3]



                        #vcf.variants.loc[i,'REFSEQ']
            
            #vcf.variants.to_csv('test.tsv', sep='\t', index=False)
            print(vcf.variants)

            os.makedirs(args.output_folder, exist_ok=True)
            vcf.write_sql_output(args.output_folder + '/' + VCF_FILE + '.txt')
            vcf.write_sql_meta_output(args.output_folder + '/' + VCF_FILE + '_meta.txt')

