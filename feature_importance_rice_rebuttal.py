import numpy as np
from scipy import stats
from scipy.stats import spearmanr
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns


from statsmodels.stats.multitest import fdrcorrection as fdr



from sklearn.feature_selection import mutual_info_regression


def getFeatures(ctrac1,loop1,ctrac2,loop2,ctrac3,loop3,ctrac4):
    ctracNum=len(ctrac1)
    imotif_seq=ctrac1+loop1+ctrac2+loop2+ctrac3+loop3+ctrac4
    imotif_len=len(ctrac1)+len(ctrac2)+len(ctrac3)+len(ctrac4)+len(loop1)+len(loop2)+len(loop3)

    all_A_density=imotif_seq.count("A")/imotif_len
    all_C_density=imotif_seq.count("C")/imotif_len
    all_G_density=imotif_seq.count("G")/imotif_len
    all_T_density=imotif_seq.count("T")/imotif_len

    allLoop_seq=loop1+loop2+loop3
    allLoop_len=len(allLoop_seq)
    allLoop_A_density=allLoop_seq.count("A")/allLoop_len
    allLoop_C_density=allLoop_seq.count("C")/allLoop_len
    allLoop_G_density=allLoop_seq.count("G")/allLoop_len
    allLoop_T_density=allLoop_seq.count("T")/allLoop_len

    loop2_A_density=loop2.count("A")/len(loop2)
    loop2_C_density=loop2.count("C")/len(loop2)
    loop2_G_density=loop2.count("G")/len(loop2)
    loop2_T_density=loop2.count("T")/len(loop2)
    
    side_loop_A_density=(loop1.count("A")+loop3.count("A"))/(len(loop1)+len(loop3))
    side_loop_C_density=(loop1.count("C")+loop3.count("C"))/(len(loop1)+len(loop3))
    side_loop_G_density=(loop1.count("G")+loop3.count("G"))/(len(loop1)+len(loop3))
    side_loop_T_density=(loop1.count("T")+loop3.count("T"))/(len(loop1)+len(loop3))

    loop1_len=len(loop1)
    loop2_len=len(loop2)
    loop3_len=len(loop3)
    maxLoop=max(len(loop1),len(loop2),len(loop3))
    minLoop=min(len(loop1),len(loop2),len(loop3))
    
    if loop1_len>loop3_len:
        max_side_loop=loop1
        min_side_loop=loop3
    else:
        max_side_loop=loop3
        min_side_loop=loop1
        
    max_side_loop_A_density=max_side_loop.count("A")/len(max_side_loop)
    max_side_loop_C_density=max_side_loop.count("C")/len(max_side_loop)
    max_side_loop_G_density=max_side_loop.count("G")/len(max_side_loop)
    max_side_loop_T_density=max_side_loop.count("T")/len(max_side_loop)
    min_side_loop_A_density=min_side_loop.count("A")/len(min_side_loop)
    min_side_loop_C_density=min_side_loop.count("C")/len(min_side_loop)
    min_side_loop_G_density=min_side_loop.count("G")/len(min_side_loop)
    min_side_loop_T_density=min_side_loop.count("T")/len(min_side_loop)
    return [
        ctracNum,
        imotif_len,
        all_A_density,
        all_C_density,
        all_G_density,
        all_T_density,
        allLoop_len,
        loop2_len,
        loop1_len+loop3_len,
        max(loop1_len,loop3_len),
        min(loop1_len,loop3_len),
        maxLoop,
        minLoop,
        allLoop_A_density,
        allLoop_C_density,
        allLoop_G_density,
        allLoop_T_density,
        loop2_A_density,
        loop2_C_density,
        loop2_G_density,
        loop2_T_density,
        side_loop_A_density,
        side_loop_C_density,
        side_loop_G_density,
        side_loop_T_density,
        max_side_loop_A_density,
        max_side_loop_C_density,
        max_side_loop_G_density,
        max_side_loop_T_density,
        min_side_loop_A_density,
        min_side_loop_C_density,
        min_side_loop_G_density,
        min_side_loop_T_density
    ]

def input_chr_file(file):
    chr_dic={}
    with open(file,"r") as f:
        sig=0
        for line in f:
            seq=line.strip()
            if seq[0]==">":
                if sig==1:
                    chr_dic[char_name]=Seq_total.upper()
                char_name=seq[1:]
                Seq_total=""
                sig=1
            else:
                Seq_total=Seq_total+seq
        chr_dic[char_name]=Seq_total.upper()
    return chr_dic

def read_TE(file):
    translation={}
    with open(file,"r") as f:
        for line in f:
            q=line.strip().split("\t")
            if q[0] != "Longest Isoforms":
                translation[q[0]]=float(q[-1])
    return translation

def read_im_result(file,sequence):
    rice_with_imotif={}
    with open(file,"r") as f:
        for line in f:
            q=line.strip().split("\t")
            if q[0]=="Putative i-motif id":
              continue
            name=q[0].split("|")
            if name[0][-1] =="+":
              if name[0][0:-1] not in sequence.keys():
                continue
              if name[0][0:-1] not in rice_with_imotif.keys():
                rice_with_imotif[name[0][0:-1]]=[]
              rice_with_imotif[name[0][0:-1]].append(name[-7:]+[float(q[-1]),len(sequence[name[0][0:-1]])-int(name[2])])
    return rice_with_imotif

def merge_sort(lst,loc):
    if len(lst)<= 1:
        return lst
    mid = len(lst)//2
    left = merge_sort(lst[:mid],loc) 
    right = merge_sort(lst[mid:],loc)
    merged = [] 
    while left and right:
        if left[0][loc] <= right[0][loc]:
            merged.append(left.pop(0))
        else:
            merged.append(right.pop(0))
    merged.extend(right if right else left) 
    return merged


    
name_dic={
    0:"C-tract length",
    1:"iM length",
    2:"A density in iM",
    3:"C density in iM",
    4:"G density in iM",
    5:"T density in iM",
    6:"Loop length",
    7:"Middle loop length",
    8:"Two side loop length",
    9:"Longest side loop length",
    10:"Shortest side loop length",
    11:"Longest loop length",
    12:"Shortest loop length",
    13:"A density in loop regions",
    14:"C density in loop regions",
    15:"G density in loop regions",
    16:"T density in loop regions",
    17:"A density in middle loop",
    18:"C density in middle loop",
    19:"G density in middle loop",
    20:"T density in middle loop",
    21:"A density in side loops",
    22:"C density in side loops",
    23:"G density in side loops",
    24:"T density in side loops",
    25:"A density in longest side loop",
    26:"C density in longest side loop",
    27:"G density in longest side loop",
    28:"T density in longest side loop",
    29:"A density in shortest side loop",
    30:"C density in shortest side loop",
    31:"G density in shortest side loop",
    32:"T density in shortest side loop",
    33:"iM strength (predicted)",
}

permutation=1000


sequence=input_chr_file("Rice_5UTR.fa")
rice_with_imotif=read_im_result("iM-seeker_final_prediction_Nip.txt",sequence)
translation=read_TE("polysome_and_RNAseq_one_isoform_RiceNip.txt")
result_file="all_the_statistics_feature_importance_rice_rebuttal.txt"


n=[]
x=[]
fil=[]
for i in rice_with_imotif.keys():
  #if len(rice_with_imotif[i])==1:
  if i in translation.keys():
    n.append(i)
    if len(rice_with_imotif[i])==1:
        tmp=rice_with_imotif[i][0]
        fil.append("-".join([tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6]]))
        x.append(i)
print(len(rice_with_imotif.keys()),len(n),len(x))
name=list(rice_with_imotif.keys())

X1=[]
trans=[]

for i in x:
    tmp=rice_with_imotif[i][0]
    if fil.count("-".join([tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6]]))>1: #filter one iM with multiple TE to aviod redundancy.
        continue
    features=getFeatures(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6])
    features.append(tmp[-2])
    X1.append(features)
    trans.append(translation[i])
print(len(trans),len(X1))

sort_trans=sorted(trans)
insert=round(len(sort_trans)/4)

sig=[]


significance_spearman={}
significance_pearson={}
for i in range(len(X1[0])):
    tmp=[]
    for j in X1:
        tmp.append(j[i])
    sp,spp=spearmanr(tmp,trans)
    pr,prp=pearsonr(tmp,trans)
    significance_spearman[i]=[sp,spp]

    significance_pearson[i]=[pr,prp]



utest_list=[]
significance_utest={}
significance_KStest={}
for i in range(len(X1[0])):
    tmp1,tmp2=[],[]
    for j in range(len(X1)):
        if trans[j]>=np.mean(trans):#sort_trans[-insert]:#
            tmp1.append(X1[j][i])
        else:

            tmp2.append(X1[j][i])

    tgt, tgp=stats.mannwhitneyu(tmp1,tmp2)

    significance_utest[i]=[i, tgt, tgp,np.mean(tmp1),np.mean(tmp2)]

    sig.append(i)
    utest_list.append(tgp)
    
print(len(tmp2),len(tmp1))

poo, v=fdr(utest_list)
print(poo, v)
for i in range(len(poo)):
    if poo[i]:
        print(name_dic[i],v[i])



permutation_utest=1000

permutation_record={}
for i in range(len(X1[0])):
    permutation_record[i]=[]
    
for rounds in range(permutation_utest):
    tmp_trans=np.random.permutation(trans)
    for i in range(len(X1[0])):
        tmp1,tmp2=[],[]
        for j in range(len(X1)):
            if tmp_trans[j]>=np.mean(tmp_trans):
                tmp1.append(X1[j][i])
            else:
                tmp2.append(X1[j][i])
        tgt, tgp=stats.mannwhitneyu(tmp1,tmp2)
        permutation_record[i].append(tgp)

utest_permutation_list=[]
for i in permutation_record.keys():
    utest_permutation_list.append(sum(1 for j in permutation_record[i] if j < utest_list[i])/permutation_utest)

for i in range(len(utest_permutation_list)):
    print(name_dic[i],utest_permutation_list[i])


X=X1


mi = mutual_info_regression(X, trans,random_state=12345)

print(mi)



p_mi={}
all_mi_permutation=[]
for rounds in range(permutation):
    tmp_trans=np.random.permutation(trans)
    mi_tmp = mutual_info_regression(X, tmp_trans,random_state=12345)
    all_mi_permutation.append(mi_tmp)

for i in range(len(all_mi_permutation[0])):
    more=0
    for j in all_mi_permutation:
        if j[i]>mi[i]:
            more=more+1
    if more/permutation<0.05 and utest_list[i]<0.05:
        print("\t".join([name_dic[i],str(mi[i]),str(more/permutation),str(significance_spearman[i][0]),
                         str(significance_spearman[i][1]),str(utest_list[i]),str(v[i]),str(utest_permutation_list[i])]))

