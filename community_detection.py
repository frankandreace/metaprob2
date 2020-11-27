#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 14:28:24 2020
#375302
@author: francescoandreace
"""

import csv
import numpy as np
import resource
import time
from scipy import sparse
import scipy as sp
from sknetwork.clustering import Louvain
import sys

print("Starting community detection.")
print("Reading input files...")
time__1 = time.time()
in_str =  str(sys.argv[1])# FILE .GFA 
reads_str = str(sys.argv[2])# STARTING FILE.FA
out_str = str(sys.argv[3])# OUTPUT FILE .FASTA
add_str = str(sys.argv[4])# CSV FILES



max_read_length = 25000 # MAX LENGTH OF OUPUT READS 
file = open(in_str,"r")
read = file.read().split("S\t")
file.close()
del file
del read[0]


read_12 = open(reads_str,"r")
read_fna_12 = read_12.read().split("\n")
read_12.close()
leng = int(len(read_fna_12)/4) #MAX INDEX 

    
generated_reads = [] # HERE I STORE ALL THE UNITIGS

cluster_csv = []
grp_numb = []

#creo due array inizializzati a 0 della lunghezza del numero delle read (per end). Assegnerò al posto di 0 il numero del gruppo a cui fan riferimento
#in questo modo so che se one[i] e two[i] sono da 2 gruppi diversi li devo concatenare. 
one = np.zeros(leng,dtype=np.uint32) # dove mi salvo in che gruppo è presente la read .1 -- '0' se in nessun gruppo
two = np.zeros(leng,dtype=np.uint32) # dove mi salvo in che gruppo è presente la read .2 -- '0' se in nessun gruppo


for i in range(len(read)):
    read[i] = read[i].split("a")
    grp =[]
    gen = read[i][0].split("\t")
    generated_reads.append(gen[1]) # MI SALVO L'UNITIG CORRISPONDENTE AL GRUPPO
    del gen
    del read[i][0]
    for k in range(len(read[i])):
        read[i][k] = read[i][k].split("\t")
    for s in range(len(read[i])):
        if  int(read[i][s][3][len(read[i][s][3])-1]) == 1: # se la read finisce con .1
            one[int(read[i][s][3][1:len(read[i][s][3])-2])-1] = i+1 #metto il valore i+1 (dei gruppi) al posto dello 0 . i+1 perchè parto dall'indice 0 e voglio distinguere il primo gruppo da quelle non raggruppate
            grp.append(int(read[i][s][3][1:len(read[i][s][3])-2]))
        else:            # se la read finisce con .2
            two[int(read[i][s][3][1:len(read[i][s][3])-2])-1] = i+1 
            grp.append(int(read[i][s][3][1:len(read[i][s][3])-2]))
   
    grp_numb.append(len(grp))
    cluster_csv.append(grp)
    
#del grp
kkk = len(read)
del read


#creo una lista di adiacenza controllando se due read paired end sono in grupppi diversi
adj_list = []
for i in range(len(one)):
    adj_list.append([])

mat_list_1 = []
mat_list_2 = []
dat = []
num_grpd = 0
num_diff = 0
num_not_grp = 0
left_out = []
#per tutta la lunghezza degli array, controllo se:
# A) entrambi sono diversi da 0 - se anche solo 1 dei 2 è ugula a 0 la read non è stata raggruppata. ora cerco solo se sono in 2 gruppi differenti
# B) se i gruppi sono differenti, aggiungo l'arco tra i due (bidirezionale / non orientato)
#per la lista di adiacenza controllo che non sia già stato inserito il nuemro del gruppo
myDic = {}
deg = np.zeros(kkk)

print("Grouping paired-end and creating graph...")
def isIn(lista, numero):
    for j in lista:
        if j == numero:
            return True
    return False

for i in range(len(one)):
    
    if one[i] == 0 and two[i] == 0:
        num_not_grp +=1 
        left_out.append(i)
    elif one[i] != 0 and two[i] != 0:
        num_grpd +=1
        if one[i] != two[i]:
            deg[one[i]-1] += 1
            deg[two[i]-1] += 1
            num_diff +=1
            
            if isIn(adj_list[one[i]-1],two[i]-1) == False:
                adj_list[one[i]-1].append(two[i]-1)
                #adj_list[two[i]-1].append(one[i]-1)
                mat_list_1.append(one[i]-1)
                
                mat_list_2.append(two[i]-1)
                
                dat.append(1)
                myDic[str(one[i]-1)+","+str(two[i]-1)] = len(mat_list_1)-1
            elif isIn(adj_list[two[i]-1],one[i]-1) == False:
                adj_list[two[i]-1].append(one[i]-1)
                mat_list_1.append(one[i]-1)
                
                mat_list_2.append(two[i]-1)
                
                dat.append(1)
                myDic[str(two[i]-1)+","+str(one[i]-1)] = len(mat_list_1)-1
            else:
                dat[int(myDic.get(str(one[i]-1)+","+str(two[i]-1)))] += 1
    else:
        x = max(one[i],two[i])
        cluster_csv[x-1].append(int(i+1))
        grp_numb[x-1] = grp_numb[x-1] + 1

in_place = 0       
for i in range(len(cluster_csv)):
    in_place += len(cluster_csv[i])


degree = []
for i in range(len(deg)):
    degree.append((int(deg[i]),i))
degree = sorted(degree, key=lambda x: x[0],reverse=True)

print("Clustering...")
mat = sparse.csr_matrix((dat,(mat_list_1,mat_list_2)),shape=(kkk,kkk))

#find the connected components and their label
n_components, labels = sp.sparse.csgraph.connected_components(csgraph=mat, directed=True, return_labels=True)


groups_id =[]
reads_out = []

for i in range(n_components):
    reads_out.append([])
    groups_id.append([])
for i in range(len(labels)):
    groups_id[labels[i]].append(i)
    reads_out[labels[i]].append(generated_reads[i])


def createDic(array):
    dic = {}
    for i in range(len(array)):
        dic[array[i]] = i
    return dic 

def createMatrix(adjlist,reads,dic,datt,dic1):
    a = list()
    b = list()
    c = list()
    for i in range(len(reads)):
        for k in range(len(adjlist[reads[i]])):
            a.append(dic.get(reads[i]))
            b.append(dic.get(adjlist[reads[i]][k]))
            c.append(dat[dic1.get(str(reads[i])+","+str(adjlist[reads[i]][k]))]) # HO MODIFICATO QUI : metto "1/"?
    #c = np.ones(len(a))            
    return a,b,c

#this method create the new groups based on the spectral clustering labels
def createClust(arrayOfLabels, dd):
    listOfLists = []
    for i in range(np.amax(arrayOfLabels)+1):
        listOfLists.append([])
    for i in range(len(arrayOfLabels)):
        listOfLists[arrayOfLabels[i]].append(dd[i])

    return listOfLists

#LA SCELTA DEL NUMERO MINIMO DI UNITIG AFFINCHE' AVVENGA IL CLUSTERING (minClusterGr) PER ORA E' TOTALMENTE ARBITRARIA... DA PENSARE SU CHE BASE SI PUO' RENDERE PIU' SENSATA.
minClusterGr = 100
myList =[] 

#CLUSTERING GRUPPI GROSSI
for i in range(len(groups_id)):
    if len(groups_id[i]) > minClusterGr:
        dic = createDic(groups_id[i])
        a,b,c = createMatrix(adj_list,groups_id[i],dic,dat,myDic)

        matn = sparse.csr_matrix((c,(a,b)),shape=(len(groups_id[i]),len(groups_id[i])))
        louvain = Louvain(random_state=0)
        out = louvain.fit_transform(matn)
            

        myList += createClust(out,groups_id[i])
   
   

len_finals = []

out_file = open(out_str,"w")
for i in range(len(myList)):
    x_l = 0
    agglomerated = ""
    j = 0
    while j <len(degree):  ## LE READ CREATE DAI SUPERGRUPPI LE CHIAMO G
        if len(myList[i]) == 0:
            break
        if degree[j][1] in myList[i]:
            k = myList[i].index(degree[j][1])
            len_new = len(generated_reads[myList[i][k]])+1
            if x_l+len_new < max_read_length:
                x_l += len_new
                agglomerated += str(generated_reads[int(myList[i][k])])+"N"  
                del degree[j]
            elif x_l == 0:
                agglomerated += str(generated_reads[int(myList[i][k])])[:max_read_length]
                del degree[j]
                break
            else:
                break
 

        j += 1
    if len(agglomerated) != 0:
        len_finals.append(len(agglomerated)) 

        out_file.write(">G"+str(i)+"\n")
        out_file.write(agglomerated[:len(agglomerated)-1]+"\n")

print("Generating output files...")
act = 0
for i in range(len(groups_id)): ## LE READ CREATE DAI GRUPPI NORMALI LE CHIAMO U
    if len(groups_id[i]) < minClusterGr:
        agglomerated = ""
        for j in range(len(groups_id[i])):
            if len(agglomerated)+len(generated_reads[int(groups_id[i][j])]) < max_read_length:
                agglomerated += str(generated_reads[int(groups_id[i][j])])+"N"
            else:
                if len(agglomerated)== 0:
                    agglomerated = generated_reads[int(groups_id[i][j])][:max_read_length]
                else:
                    break
        
        out_file.write(">U"+str(act)+"\n")
        out_file.write(agglomerated[:len(agglomerated)-1]+"\n")
        act+=1

for i in range(len(left_out)):    
    out_file.write(read_fna_12[left_out[i]*2]+"\n")
    out_file.write(read_fna_12[left_out[i]*2+1]+"N"+read_fna_12[(left_out[i]+leng)*2+1]+"\n")

out_file.close()

#CONTA FINALE DELLE READS IN OGNI GRUPPO
number_final_G = []
number_final_U = []

for i in range(len(myList)):
    summ = 0
    for j in range(len(myList[i])):
        summ += grp_numb[myList[i][j]]
    number_final_G.append(summ)
for i in range(len(groups_id)):
    if len(groups_id[i]) < minClusterGr:
        summ = 0
        for j in range(len(groups_id[i])):
            summ += grp_numb[groups_id[i][j]]
        number_final_U.append(summ)

with open(add_str+"_G.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(number_final_G)
with open(add_str+"_U.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(number_final_U)

#CREO CSV CON, IN GR_UN, TUTTE LE READS IN OGNI GRUPPO DI UNITIGS E, IN SN_UN, TUTTE LE READS IN OGNI UNITIG RIMASTO SINGOLO O IN UN GRUPPO PICCOLO.
gr_un = []
for i in range(len(myList)):
    gr_un.append([])
    for j in range(len(myList[i])):
        gr_un[i].extend(cluster_csv[myList[i][j]])
del myList
sn_un = []
for i in range(len(groups_id)):
    if len(groups_id[i]) < minClusterGr:
        sn_un.append([])
        for j in range(len(groups_id[i])):
            sn_un[len(sn_un)-1].extend(cluster_csv[groups_id[i][j]])
del groups_id

with open(add_str+"_groups.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(gr_un)

with open(add_str+"_unitigs.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(sn_un)

time__2 = time.time()
print("Elapsed time: "+str(time__2-time__1))
print("Max memory: "+str(round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000)))