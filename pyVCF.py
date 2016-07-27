#!/usr/bin/env python3
# encoding: utf8

'''
VariantChecker

@author: Yannick Boursin

@license: GNU GPLv3

@organization: Gustave Roussy

@contact: yannick.boursin@gustaveroussy.fr

@created: 01/03/2016
'''

import argparse
from collections import defaultdict, Counter
import pysam

global vc

parser = argparse.ArgumentParser(description='')
parser.add_argument('--vcf', dest="vcf", nargs='+', help='Please give in input your VCF files following this guidelines: <Normal|Tumor>:<SampleID>:<Path to VCF file>')
parser.add_argument('--bam', dest="bam", nargs='+', help='Please give in input your BAM files following this guidelines: <Normal|Tumor>:<SampleID>:<Path to your BAM files>')
parser.add_argument('--variantCaller', dest="vc", required=True, help="Please help this poor program by providing insights about your variant caller (we could detect it but it wouldn't be fun). Available: varscan, mutect, gatk")
parser.add_argument('--output-type', dest="ot", required=False, default="quantumclone", help="Switches between sciclone and quantumclone. Type sciclone for sciclone output, quantumclone for quantumclone output (isn't that so easy ?)")
parser.add_argument('-o', '--out', dest="o", required=True, default="pyVCF.out", help="Please give some output filename")
parser.add_argument('--noindel', dest='noindel', required=False, default=True)
args = parser.parse_args()
vcf = args.vcf
bam = args.bam
vc = args.vc
ot = args.ot

noindel = args.noindel

sdico = defaultdict(lambda: defaultdict(list)) # Dico for samples

def parse_arguments(arg, sdico):
    dico = {"normal": defaultdict(str), "tumor": defaultdict(str)}
    for a in arg:
        a = a.split(':')
        sample_type = a[0].lower()
        sample_id = a[1]
        path = a[2]
        dico[sample_type][sample_id] = path
        sdico[sample_id][sample_type] += [path]
    return dico, sdico

def getPathDico(vcf, bam):
    dic = defaultdict(bool)
    for a in vcf:
        a = a.split(':')
        path=a[2]
        id=":".join([a[1], a[0].lower()])
        dic[path] = id
    for a in bam:
        a = a.split(':')
        path=a[2]
        id=":".join([a[1],a[0].lower()])
        dic[path] = id
    return dic

dico_paths = getPathDico(vcf, bam)
print(dico_paths)
dico_vcf, sdico = parse_arguments(vcf, sdico)
dico_bam, sdico = parse_arguments(bam, sdico)


if ot == "quantumclone":
    oH = open(args.o, "w")
else:
    oH = {}
    for k,v in sdico.items():
        for k2, v2 in v.items():
            id = '{}:{}'.format(k, k2)
            file = open(args.o + "_{}_{}".format(k, k2), "a")
            oH[id] = file


def inbam(bam, chr, pos):
    bamfile = pysam.AlignmentFile(bam, "rb")
    count = bamfile.count(reference=chr, start=pos-1, end=pos+1)
    if count == 0:
        return Counter(), chr, pos, 0
    bases = []
    for pilCol in bamfile.pileup(reference=chr, start=pos-1, end=pos+1):
        for pilRead in pilCol.pileups:
            if pilCol.pos == pos-1:
                if not pilRead.is_del and not pilRead.is_refskip:
                    bases.append(pilRead.alignment.query_sequence[pilRead.query_position])
                else:
                    if pilRead.is_del: bases.append('del')
                    if pilRead.is_refskip: bases.append('refskip')
                pilPos = pilCol.pos
    results = Counter(bases), chr, pos, count
    #print(results)
    return results

class Variant(object):
    def __init__(self, chr, pos, ref, alt, line=None, sample=None):
        self.chr = chr
        self.pos = pos 
        self.ref = ref
        self.alt = alt
        self.sample = sample
        self.depth = defaultdict(lambda: defaultdict(bool))
        self.bamDepth = defaultdict(lambda:defaultdict(bool))
        
    def __eq__(self, other):
        if self.chr == other.chr:
            if self.pos == other.pos:
                if self.alt == other.alt:
                    return True
        return False
    
    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}\n'.format(self.chr, self.pos, self.ref, self.alt, str(self.sample))
    
    def getInfos(self, sample, line):
        if vc == "varscan":
            if "tumor" in sample:
                gF = line[10]
            elif "normal" in sample:
                gF = line[9]
            gF = gF.split(':')
            self.depth[sample]["totdp"] = int(gF[2])
            self.depth[sample]["refdp"] = int(gF[3])
            self.depth[sample]["altdp"] = int(gF[4])
        elif vc == "mutect":
            if "tumor" in sample:
                gF = line[10]
            elif "normal" in sample:
                gF = line[9]
            gF = gF.split(':')
            self.depth[sample]["refdp"] = int(gF[0])
            self.depth[sample]["altdp"] = int(gF[1])
            self.depth[sample]["totdp"] = self.altdp + self.refdp
        elif vc == "gatk":
            gF = line[-1].rstrip('\n')
            #print(gF)
            gF = gF.split(':')
            ad = gF[1].split(',')
            gt = gF[0].split(',')
            #print(gt)
            altallele = int(gt[0].split('/')[-1])
            self.depth[sample]["refdp"] = int(ad[0])
            self.depth[sample]["altdp"] = int(ad[altallele])
            self.depth[sample]["totdp"] = int(gF[2])
        else:
            pass
        
    def pushBamInfos(self, bamId, **kwargs):
        #print(kwargs)
        for k, v in kwargs.items():
            self.bamDepth[bamId][k] = v
        #print(self.bamDepth)
        # Now we must id the refdepth and the altdepth in those. We will use our infos: v.ref, v.alt
        if not self.depth[bamId]['totdp']:
            self.depth[bamId]['totdp'] = self.bamDepth[bamId]['DP']
            if type(self.bamDepth[bamId][self.ref]) is int: 
                self.depth[bamId]['refdp'] = self.bamDepth[bamId][self.ref]
            else: self.depth[bamId]['refdp'] = 0
            if type(self.bamDepth[bamId][self.alt]) is int:
                self.depth[bamId]['altdp'] = self.bamDepth[bamId][self.alt]
            else: self.depth[bamId]['altdp'] = 0
        self.bamDepth[bamId]['totdp'] = self.bamDepth[bamId]['DP']
        if type(self.bamDepth[bamId][self.ref]) is int: self.bamDepth[bamId]['refdp'] = self.bamDepth[bamId][self.ref]
        else: self.bamDepth[bamId]['refdp'] = 0
        if type(self.bamDepth[bamId][self.alt]) is int: self.bamDepth[bamId]['altdp'] = self.bamDepth[bamId][self.alt]
        else: self.bamDepth[bamId]['altdp'] = 0

# Now we must parse VCF files and BAM files by sample and retain only variants that are covered in bams (dp à définir) and present in vcf
# Au premier passage, on importe tous les variants de tous les VCF, en traçant les samples.
# Ensuite, on verifiera que chacun de ces variants est couvert dans les bams
# Puis on sortira tous les variants selon le format voulu
variants = []
variants_dic = defaultdict(bool)
counter = 0
sampleDico = defaultdict(list)

for k,v in sdico.items():
    for k2, v2 in v.items():
        sampleDico[k].append(k2)
print(sampleDico)

for k,v in sdico.items():
    # k = sample_name
    # v = dict(tumor / Normal)
    # in v: list: (vcf then bam)
    for k2, v2 in v.items():
        # k2 = tumor or Normal
        # v2 = pathToVCF (index0) or pathToBam (index1)
        if len(v2) > 2: 
            exit("There cannot be more than two files associated with sample {0} for {1} condition".format(k, k2))
        cvcf = v2[0]
        cbam = v2[1]
        with open(cvcf, "r") as vcfh:
            for line in vcfh.readlines():
                if not line.startswith('#'):
                    counter +=1
                    if counter % 5000 == 0: print('Processed: {} variants'.format(counter))
                    line = line.split('\t')
                    chr = line[0]
                    pos = line[1]
                    ref = line[3]
                    alt = line[4]
                    if len(ref) > 1 or len(alt) > 1:
                        continue
                    var = Variant(chr, pos, ref, alt, sample=sampleDico)
                    var.getInfos("{}:{}".format(k,k2), line)
                    varhash = '{}-{}'.format(chr, pos)
                    present = variants_dic[varhash]
                    if not present:
                        variants.append(var)
                        variants_dic[varhash] = len(variants) - 1
                    else:
                        var = variants_dic[varhash]
                        var = variants[var]
                        #if k in var.sample.keys():
                        #    var.sample[k].append(k2)
                        #else:
                        #    var.sample[k] = [k2]
                        var.getInfos("{}:{}".format(k,k2), line)
        counter = 0
        print('Done {}:{}'.format(k, k2))
print('Total variants: {}'.format(len(variants)))


#print('For test purposes: keeping only first thousand variants')
#variants = variants[:100] # For test purpose.


bam_list = []
for k,v  in dico_bam.items():
    for k2,v2 in v.items():
        bam_list.append(v2)
#print(bam_list)

currated_variants = []
discarded_variants = []
counter = 0
total = len(variants) * len(bam_list)

for v in variants:
    if counter % 100 == 0: print('Requests in BAM file: {} / {}'.format(counter, total))
    chr = v.chr
    start = int(v.pos)
    filter = False
    for bam in bam_list:
        #print(bam)
        dp_res = inbam(bam, chr, start)
        #print(dp_res)
        dp = dp_res[3]
        #print('{}:{} {} {} {}'.format(v.chr, v.pos, v.ref, v.alt, v.depth))
        #print('{}:{} {} {}'.format(dp_res[1], dp_res[2], dp_res[0], dp_res[3]))
        try:
            toPush = dp_res[0]
            toPush['DP'] = dp_res[3]
        except:
            import pdb
            pdb.set_trace()
        #print(toPush)
        v.pushBamInfos(dico_paths[bam], **toPush) 
        #print('{}:{} {} {} {}'.format(v.chr, v.pos, v.ref, v.alt, v.depth))
        counter += 1
        if dp is None or dp == 0:
            filter = True
            print('Discarding current variant: {}-{} {}->{} as position is not covered enough in {}; {} reads'.format(v.chr, v.pos,v.ref, v.alt, dico_paths[bam], dp))
            discarded_variants.append(v)
            break
    if not filter:
        #print(v.chr, v.pos, v.ref, v.alt, v.depth, v.bamDepth)
        currated_variants.append(v)
        
        #print('Keeping this variant ! {}'.format(v))
print('Kept: {} variants'.format(len(currated_variants)))
print('Trash: {} variants'.format(len(discarded_variants)))

for v in currated_variants:
    # currated_variants is a list of variant objects
    for k, d in v.sample.items():
        #print(k,d, len(k), len(d))
        # v.sample is a defaultdict of lists with sampleIds in keys, and tumor/normal status in values
        for s in d:
            
            #d is a list of values "tumor/normal"
            # s == normal or s == tumor
            # k,s == sample, tumor 
            #print(v.depth)
            id = "{}:{}".format(k,s)
            dico_depth = v.depth[id]
            if ot == "quantumclone":
                print('{}:{}\t{}\t{}\t{}\t{}'.format(k,s, v.chr, v.pos, dico_depth["totdp"], dico_depth["altdp"]), file=oH)
            elif ot == "sciclone":
                file = oH[id]
                if dico_depth["totdp"] != 0 and dico_depth["totdp"] is not False:
                    freq = float(dico_depth["altdp"]) / (float(dico_depth["refdp"]) + float(dico_depth["altdp"]))
                else: freq = False
                print('{}\t{}\t{}\t{}\t{}'.format(v.chr, v.pos, dico_depth["refdp"], dico_depth["altdp"], freq), file=file)

if type(oH) is list:
    [x.close() for x in oH]
else:
    oH.close()
