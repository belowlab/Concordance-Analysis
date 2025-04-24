import sys
###quick note, script will only work in python3. Assumes ordered dictionary

sim1 = {} # dic[rsid] = [+/-,1/0] #1 indicates new block
sim2 = {}

file1 = open(sys.argv[1],"r")
file2 = open(sys.argv[2],"r")

pval_thresh = float(sys.argv[3]) #change pvalue threshold as needed, used 0.005 for publication

###first summary statistics file
file1.readline()
first = True
print("reading file 1")
for line in file1:
    spline = line.split()
    if(True):
        ###create sim dic
        meets_pval = float(spline[8].rstrip())<pval_thresh ##### CHANGE spline[index] to the correct pvalue column
        if(first):
            pos = int(spline[2]) ####### CHANGE spline[index] to column with the basepair position
            if(float(spline[6])>0): ###### CHANGE spline[index] to the BETA column
                sign = "+"
                sim1[spline[0]] = [sign,1,meets_pval] ##### CHANGE spline[inxdex to the SNP column
            else:
                sign = "-"
                sim1[spline[0]] = [sign,1,meets_pval] ##### CHANGE spline[inxdex to the SNP column
            first = False
        else:
            if(float(spline[6])>0): ###### CHANGE spline[index] to the BETA column
                new_sign = "+"
            else:
                new_sign = "-"
            if((int(spline[2])-int(pos))<10000 and sign == new_sign): ##if still in ld block (same direction of effect and within 10kb ####### CHANGE spline[index] to column with the basepair position
                sim1[spline[0]] = [sign,0,meets_pval] ##### CHANGE spline[inxdex] to the SNP column
                pos = spline[2] #update position ####### CHANGE spline[index] to column with the basepair position
            else:#if new ld block
                sim1[spline[0]] = [new_sign,1,meets_pval] ###CHANGE spline[inxdex] to the SNP column
                sign = new_sign ##update to the new sign
                pos = spline[2] #update position  ####### CHANGE spline[index] to column with the basepair position

###Prune dictionary and update blocks
print("Pruning dictionary 1...")
update_next = False
for rsid in list(sim1.keys()):
    if(sim1[rsid][2]==False):
        if(sim1[rsid][1]==1):
            update_next = True
        del sim1[rsid]
    else:
        if(update_next):
            sim1[rsid][1]=1 ###if the start of a block gets removed, update the next position to become the new start of the block
        update_next = False

###second summary statistics file
file2.readline()
first = True
print("reading file 2")
for line in file2:
    spline = line.split()
    if(True):
        ###create sim dictionary
        meets_pval = float(spline[8].rstrip())<pval_thresh ##### CHANGE spline[index] to the correct pvalue column
        if(first):###First read just establishes the starting position and sign
            pos = int(spline[2]) ####### CHANGE spline[index] to column with the basepair position
            if(float(spline[6])>0):  ###### CHANGE spline[index] to the BETA column
                sign = "+"
                sim2[spline[0]] = [sign,1,meets_pval] ##### CHANGE spline[inxdex to the SNP column
            else:
                sign = "-"
                sim2[spline[0]] = [sign,1,meets_pval] ##### CHANGE spline[inxdex to the SNP column
            first = False
        else:
            if(float(spline[6])>0): ###### CHANGE spline[index] to the BETA column
                new_sign = "+"
            else:
                new_sign = "-"
            if((int(spline[2])-int(pos))<10000 and sign == new_sign): #if still in ld block (same direction of effect and within 10kb.. #######AS: CHANGE spline[index] to column with the basepair position
                sim2[spline[0]] = [sign,0,meets_pval] ##### CHANGE spline[index] to the SNP column
                pos = spline[2] #update position  ####### CHANGE spline[index] to column with the basepair position
            else: #if new ld block
                sim2[spline[0]] = [new_sign,1,meets_pval] ##### CHANGE spline[inxdex to the SNP column
                sign = new_sign #update to the new sign
                pos = spline[2] #update position ####### CHANGE spline[index] to column with the basepair position


print("Pruning dictionary 2...")
update_next = False
for rsid in list(sim2.keys()): 
    if(sim2[rsid][2]==False): #check if SNP does not meet the p-value threshold
        if(sim2[rsid][1]==1): #check if it's the start of an LD block
            update_next = True #mark to update the next SNP as the new start of the block
        del sim2[rsid] #delete the SNP from the dictionary
    else:
        if(update_next): 
            sim2[rsid][1]=1 ###if the start of a block gets removed, update the next position to become the new start of the block
        update_next = False #reset the marker

###calculate the concordance before we start shuffling the dictionaries
real_total = 0
real_overlap = 0

for rsid in sim1:
    if(sim1[rsid][1]==1): #check if it's the start of an LD block in sim1
        if(rsid in sim2): #check if the rsid also exists in sim2
            real_total+=1 #increase total
            if(sim1[rsid][0]==sim2[rsid][0]): #check if the sign of effect matches
                real_overlap+=1 #increase total
    else: #if not the start of an LD block
        if(rsid in sim2): #check if the rsid exists in sim2
            if(sim2[rsid][1]==1): #check if it's the start of an LD block in sim2
                real_total+=1 #increase total
                if(sim1[rsid][0]==sim2[rsid][0]):
                    real_overlap+=1 #increase total

###simulate concordance analysis to establish an expected percentage of overlap
sim = 25 #run the concordance analysis x times 
conc_rates = []
import random
for i in range(0,sim):
    print("sim:",i)
    ###randomize the dictionaries
    ##rand dic 1
    for rsid in sim1:
        if((sim1[rsid])[1]): ###if its a new block
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
            (sim2[rsid])[0] = rand_sign #just update it to the previous randomly drawn sigh

    ###find concordance rate with randomized dictionaries, block by block
    total = 0
    overlap = 0
    for rsid in sim1:
        if(sim1[rsid][1]==1):
            if(rsid in sim2):
                total+=1
                if(sim2[rsid][0]==sim1[rsid][0]):
                    overlap+=1
        else:
            if(rsid in sim2):
                if(sim2[rsid][1]==1):
                    total+=1
                    if(sim2[rsid][0]==sim1[rsid][0]):
                        overlap+=1

    conc_rates.append(overlap/total)



print(conc_rates)
expected_concordance = sum(conc_rates)/len(conc_rates)


print("total:",real_total)
print("overlap:",real_overlap)
print("concordance:", real_overlap/real_total)
print("expected concordance rate based on simluation:", expected_concordance)
print("pval threshold:",pval_thresh)

from scipy import stats
con_p = stats.binom_test(real_overlap, n=real_total, p=expected_concordance, alternative='greater')
print("pvalue: ",con_p)


f = open("23andMe_concordance.csv","a") #open results file
f.write("\n")
f.write(sys.argv[1])
f.write(",")
f.write(sys.argv[2])
f.write(",")
f.write(sys.argv[3])
f.write(",")
f.write(str(real_total))
f.write(",")
f.write(str(expected_concordance))
f.write(",")
f.write(str(real_overlap/real_total))
f.write(",")
f.write(str(con_p))
f.close()


file1.close()
file2.close()
print("pvalue: ",con_p)


file1.close()
file2.close()
