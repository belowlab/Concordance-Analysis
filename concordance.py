import sys
###quick note, script will only work in python3. Assumes ordered dictionary
file1_dic = {}
file2_dic = {}

sim1 = {} # dic[rsid] = [+/-,1/0] #1 indicates new block
sim2 = {}

file1 = open(sys.argv[1],"r")
file2 = open(sys.argv[2],"r")

pval_thresh = .05 #### CHANGE TO THE PVALUE THRESHOLDS AS NEEDED

file1.readline()
first = True
print("reading file 1")
for line in file1:
    spline = line.split()
    if(float(spline[8].rstrip())<pval_thresh): ##### CHANGE spline[index] to the correct pvalue column
        ###create sim dic
        if(first):
            pos = int(spline[2]) ####### CHANGE spline[index] to column with the basepair position
            if(float(spline[6])>0): ###### CHANGE spline[index] to the BETA column
                sign = "+"
                sim1[spline[0]] = [sign,1] ##### CHANGE spline[inxdex to the SNP column
            else:
                sign = "-"
                sim1[spline[0]] = [sign,1] ##### CHANGE spline[inxdex to the SNP column
            first = False
        else:
            if(float(spline[6])>0): ###### CHANGE spline[index] to the BETA column
                new_sign = "+"
            else:
                new_sign = "-"
            if((int(spline[2])-int(pos))<10000 and sign == new_sign):#if still in ld block (same direction of effect and within 10kb ####### CHANGE spline[index] to column with the basepair position
                sim1[spline[0]] = [sign,"0"] ##### CHANGE spline[inxdex] to the SNP column
                pos = spline[2] #update position ####### CHANGE spline[index] to column with the basepair position
            else:#if new ld block
                sim1[spline[0]] = [new_sign,"1"] #####AS: CHANGE spline[inxdex] to the SNP column
                sign = new_sign #update to the new sign
                pos = spline[2] #update position  ####### CHANGE spline[index] to column with the basepair position


        ###create true dictionary
        if(float(spline[6])>0): ###### CHANGE spline[index] to the BETA column
            file1_dic[spline[0]] = "+" ##### CHANGE spline[index] to the SNP column
        else:
            file1_dic[spline[0]] = "-" ##### CHANGE spline[index] to the SNP column



file2.readline()
first = True
print("reading file 2")
for line in file2:
    spline = line.split()
    if(float(spline[8].rstrip())<pval_thresh):  ##### CHANGE spline[index] to the correct pvalue column
        ###create sim dictionary
        if(first):###First read just establishes the starting position and sign
            pos = int(spline[2]) ####### CHANGE spline[index] to column with the basepair position
            if(float(spline[6])>0): ###### CHANGE spline[index] to the BETA column
                sign = "+"
                sim2[spline[0]] = [sign,1] ##### CHANGE spline[inxdex to the SNP column
            else:
                sign = "-"
                sim2[spline[0]] = [sign,1] ##### CHANGE spline[inxdex to the SNP column
            first = False
        else:
            if(float(spline[6])>0):  ###### CHANGE spline[index] to the BETA column
                new_sign = "+"
            else:
                new_sign = "-"
            if((int(spline[2])-int(pos))<10000 and sign == new_sign):#if still in ld block (same direction of effect and within 10kb.. #######AS: CHANGE spline[index] to column with the basepair position
                sim2[spline[0]] = [sign,"0"] ##### CHANGE spline[index] to the SNP column
                pos = spline[2] #update position  ####### CHANGE spline[index] to column with the basepair position
            else:#if new ld block
                sim2[spline[0]] = [new_sign,"1"] ##### CHANGE spline[inxdex to the SNP column
                sign = new_sign #update to the new sign
                pos = spline[2] #update position ####### CHANGE spline[index] to column with the basepair position




        ###create true dictionary
        if(float(spline[6])>0): ###### CHANGE spline[index] to the BETA column
            file2_dic[spline[0]] = "+" ##### CHANGE spline[index] to the SNP column
        else:
            file2_dic[spline[0]] = "-" ##### CHANGE spline[index]to the SNP column


###simulate concordance analysis to establish an expected percentage of overlap
sim = 25 #run the concordance analysis x times 
conc_rates = []
import random
for i in range(0,sim):
    print("sim:",i)
    ###randomize the dictionaries
    ##rand dic 1
    for rsid in sim1:
        if((sim1[rsid])[1]):###if its a new block... which the first one should be    
            if(random.randint(0,1)): ###flip a coin
                rand_sign = "+"
            else:
                rand_sign = "-"
            (sim1[rsid])[0] = rand_sign #update dictionary with random sign
        else: ###if its not a new block...
            (sim1[rsid])[0] = rand_sign #just update it to the previous randomly drawn sign

    ###rand dic 2
    for rsid in sim2:
        if((sim2[rsid])[1]): ###if its a new block... which the first one should be
            if(random.randint(0,1)): ###flip a coin
                rand_sign = "+"
            else:
                rand_sign = "-"
            (sim2[rsid])[0] = rand_sign #update dictionary with random sign
        else: ###if its not a new block...
            (sim2[rsid])[0] = rand_sign #just update it to the previous randomly drawn sign

    ###find concordance rate with randomized dictionaries
    total = 0
    overlap = 0
    for rsid in sim1:
        if(rsid in sim2):
            total+=1
            if((sim1[rsid])[0]==(sim2[rsid])[0]):
                overlap+=1
    conc_rates.append(overlap/total)

total = 0
overlap = 0
for rsid in file1_dic:
    if(rsid in file2_dic):
        total +=1
        if(file1_dic[rsid] == file2_dic[rsid]):
            overlap += 1
print(conc_rates)
expected_concordance = sum(conc_rates)/len(conc_rates)


print("total:",total)
print("overlap:",overlap)
print("concordance:", overlap/total)
print("expected concordance rate based on simluation:", expected_concordance)
print("pval threshold:",pval_thresh)

from scipy import stats
con_p = stats.binom_test(overlap, n=total, p=expected_concordance, alternative='greater')
print("pvalue: ",con_p)


file1.close()
file2.close()
