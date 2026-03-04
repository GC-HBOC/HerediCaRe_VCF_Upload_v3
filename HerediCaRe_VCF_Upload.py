import argparse
import os # required to get *.py and resources folder
import shutil
import sys
import chardet
#import pysam
import subprocess
from VCF import VCF
import re

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
#print(TRANSCRIPTS)

#from easyfasta import load_fasta
import fastapy
#import easyfasta

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
            vcf = VCF(True)
            #normalize_vcf(os.path.join(args.input_folder, VCF), hg38_dict, HG38_FLAG)
        if not HG38_FLAG:
            vcf = VCF(False)
            #normalize_vcf(os.path.join(args.input_folder, VCF), hg19_dict, HG38_FLAG)

            ### further infos required for SQL output

        vcf.MEMBER_ID = VCF_FILE.split('-')[1]
        vcf.BOGEN_NR = VCF_FILE.split('-')[2]
        vcf.ERFMIT =  VCF_FILE.split('-')[3]
        vcf.ERFDAT =  VCF_FILE.split('-')[3]



        
        with open(os.path.join(args.input_folder, VCF_FILE), encoding=chardet.detect(rawdata)['encoding']) as infile:
            FAIL_FLAG = False
            for _l in infile:
                line = _l.rstrip().strip('"')
                if len(line):
                    if line.startswith('#'):
                        vcf.header.append(line)
                        #try:
                    else:
                        ll = line.rstrip().split('\t')
                        if len(ll) not in [8, 10]:
                            sys.stderr.write('...invalid number of columns in VCF file ' + VCF_FILE + ': ' + str(len(ll)) + '\n')
                            FAIL_FLAG = True
                            break
                        try:
                            CHROM, POS, REF, ALT, INFO = ll[0], ll[1], ll[3], ll[4], ll[7]
                        except:
                            print(ll)
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
                        
                        
                        ## split variant in single-ALT variants
                        if not FAIL_FLAG:
                            print(ll)
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
                                        print(VCF_FILE, CHROM, POS, REF, ALT, _ALT, GT, ANNOT_TAG, varclass )
                                        # ['chrom', 'pos_hg38', 'ref_hg38', 'alt_hg38', 'pos_hg19', 'ref_hg19', 'alt_hg19', 'gene',  'transcript', 'hgvsc', 'hgvsp', 'effect', 'annotation', 'class', 'gt', 'norm_fail', 'ref_fail', 'liftover_fail']
                                        if HG38_FLAG:
                                            ## TODO REF check
                                            if REF.upper() == hg38_dict[CHROM][int(POS)-1:int(POS)+len(REF)-1].upper():
                                            
                                                VAR = [CHROM, POS, REF, _ALT, None, None, None, None, None, None, None, None, ANNOT_TAG, varclass, GT, None, False, None]
                                            else:
                                                VAR = [CHROM, POS, REF, _ALT, None, None, None, None, None, None, None, None, ANNOT_TAG, varclass, GT, None, True, None]
                                        else:
                                            print('checking hg19:', CHROM, POS, REF, ALT)
                                            print(int(POS)-1, int(POS)+len(REF)-1)
                                            print(hg19_dict[CHROM][int(POS)-1:int(POS)+len(REF)-1].upper())
                                            if REF.upper() == hg19_dict[CHROM][int(POS)-1:int(POS)+len(REF)-1].upper():
                                                VAR = [CHROM, None, None, None, POS, REF, _ALT, None, None, None, None, None, ANNOT_TAG, varclass, GT, None, False, None]
                                            else:
                                                VAR = [CHROM, None, None, None, POS, REF, _ALT, None, None, None, None, None, ANNOT_TAG, varclass, GT, None, True, None]
                                        vcf.variants.loc[len(vcf.variants)] = VAR
                                    else:
                                        sys.stderr.write("Variant with ALT " + _ALT + " at " + CHROM + ':' + str(POS) + ' is ignored\n')
                        #if ',' in ALT:
                        #    print('XXX', line)
                        #    # XXX 1   109792734       .       A       AACG,AC .       .       .       GT      1/2
                             # XXX 12  121434630       .       C       CC,CTCATTCAT    .       .       .       GT      1/2
                                
                        #except:
                        #    FAIL_FLAG = True
        print(vcf.variants)
        print('FAIL_FLAG:', FAIL_FLAG)
        
        # set FAIL_FLAG if any ref_fail == True
        if vcf.variants['ref_fail'].any(): 
            print('FAILED!!!!') # TODO output 1st failed variant
            print(vcf.variants['ref_fail'])
            FAIL_FLAG = True
            
        ### NORMALIZATION
        if HG38_FLAG:
            vcf.normalize(hg38_dict)
        else:
            vcf.normalize(hg19_dict)
        print('NORMALIZED')
        
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
            CMD = args.jp + ' -Xmx' + str(args.ram) + 'g -jar ' + args.sp + ' hg38 ' + 'tmp\\'  + VCF_FILE + '.tmp.vcf'
            #print(CMD)
            output = None
            try:
                output = subprocess.check_output(CMD, shell=True, text=True)
                #print(output)
            except:
                sys.stderr.write('..Could not run ' + CMD + '\n')
            
            if output:
                DEL_INDS = [] # store indices of variants not located within pre-defined transcripts
                for line in output.split('\n'):
                    if not line.startswith('#') and len(line.split('\t')) >= 8:
                        ll = line.split('\t')
                        #print(ll)
                        CHROM, POS, REF, ALT, INFO = ll[0][3:], int(ll[1]), ll[3], ll[4], ll[7]
                        #print(CHROM, POS, REF, ALT)
                        _inds = vcf.variants.index[(vcf.variants['chrom'] == CHROM) & (vcf.variants['pos_hg38'] == POS) & (vcf.variants['ref_hg38'] == REF) & (vcf.variants['alt_hg38'] == ALT)].tolist()
                        print(_inds)
                        if len(_inds) == 1:
                            ### NOTE: due to self-generated VCF input, ANN is the only entry in INFO column
                            ANN = [_ for _ in INFO[4:].split(',') if (_.split('|')[6].split('.')[0] in TRANSCRIPTS.keys())]
                            print(ANN)
                            if len(ANN) == 1:
                                # 'gene',  'transcript', 'hgvsc', 'hgvsp', 'effect',
                                vcf.variants.loc[_inds[0], 'gene'] = ANN[0].split('|')[3]
                                vcf.variants.loc[_inds[0], 'transcript'] = ANN[0].split('|')[6]
                                vcf.variants.loc[_inds[0], 'hgvsc'] = ANN[0].split('|')[9]
                                vcf.variants.loc[_inds[0], 'hgvsp'] = ANN[0].split('|')[10]
                                vcf.variants.loc[_inds[0], 'effect'] = ANN[0].split('|')[1]



                            
                            elif len(ANN) > 1:
                                #TODO add new line in vcf.variants
                                pass
                            else:
                                DEL_INDS.append(_inds[0])
                        else:
                            pass
                            #TODO variant not found or doubled
            vcf.variants.to_csv('test.tsv', sep='\t', index=False)
            print(vcf.variants)

            os.makedirs(args.output_folder, exist_ok=True)
            vcf.write_sql_output(args.output_folder + '/' + VCF_FILE + '.txt')

