import numpy as np
import scipy as sp
import sys
import pandas as pd
import multiprocessing as mp
import csv

# utils #
def isIn(group,elem):
    for el in group:
        if el == elem:
            return True
    return False

def max_value(inputlist):
    return max([max(sublist) for sublist in inputlist])

def num_sp(precisions):
    sp = []
    for el in precisions:
        sp.append(el[0])
    return len(set(sp))

#### GFA FILE HANDLING ####

def elab_gfa_single(gfa):
    """
    args: 1)gfa file;
    returns: 1) utgs; 2) reads in each UTG; 
    """
    #saving groups
    unitigs = [] #SAVES UNITIGS
    unitigs_gfa = [] # SAVES READS IN EACH UTG

    ### READING GFA AND GETTING UTGS AND READS INSIDE
    
    for i in range(len(gfa)):
        gfa[i] = gfa[i].split("a") # a is the first char of every read line in each utg
        grp =[]
        gen = gfa[i][0].split("\t") #separating tab values
        try:
            unitigs.append(gen[1]) #getting UTG
            #print(gen[1])
            #print('-')
        except IndexError:
            print(i)
            print(gen)
            break
        del gen
        del gfa[i][0]

        #for every line in part of gfa that represents one utgs, split and get reads
        for k in range(len(gfa[i])): 
            gfa[i][k] = gfa[i][k].split("\t") 
        for s in range(len(gfa[i])):
            read_num = int(gfa[i][s][3][1:len(gfa[i][s][3])-2])
            if isIn(grp,read_num) == False:
                grp.append(read_num)

        unitigs_gfa.append(grp)

    return unitigs,unitigs_gfa

def read_gfa_single(gfafile):
    """
    args: path to gfa file
    returns:  1) utgs; 2) reads in each UTG; 
    """
    f = open(gfafile,"r")
    gfas = f.read().split("S\tutg")
    f.close()
    del f
    del gfas[0]

    unitigs, reads_in = elab_gfa_single(gfas)
    return unitigs, reads_in

def elab_gfa(gfa):
    """
    args: 1)gfa file;
    returns: 1) utgs; 2) reads in each UTG; 3) adjacency list between utgs; 
    """
    #saving groups
    unitigs = [] #SAVES UNITIGS
    unitigs_gfa = [] # SAVES READS IN EACH UTG
    adj_utgs = [[],[],[]] #adjacency list between utgs
    weigth_list = []

    ### READING GFA AND GETTING UTGS AND READS INSIDE
    
    for i in range(len(gfa)):
        gfa[i] = gfa[i].split("a") # a is the first char of every read line in each utg
        grp =[]
        weig = []
        gen = gfa[i][0].split("\t") #separating tab values
        try:
            unitigs.append(gen[1]) #getting UTG
            #print(gen[1])
            #print('-')
        except IndexError:
            print(i)
            print(gen)
            break
        del gen
        del gfa[i][0]

        #for every line in part of gfa that represents one utgs, split and get reads
        for k in range(len(gfa[i])): 
            gfa[i][k] = gfa[i][k].split("\t") 
        for s in range(len(gfa[i])):
            read_num = int(gfa[i][s][3][1:len(gfa[i][s][3])-2])
            if isIn(grp,read_num) == False:
                grp.append(read_num)
                w = int(gfa[i][s][5].split('\n')[0])
                weig.append(w)

        unitigs_gfa.append(grp)
        weigth_list.append(weig)
    #CREATING UTGS ADJACENCY LIST

    max_read = max_value(unitigs_gfa)
    check = np.zeros(max_read,dtype=np.uint32)

    for i in range(len(unitigs_gfa)):
        for j in range(len(unitigs_gfa[i])): 
            read = unitigs_gfa[i][j]
            unitigs_gfa[i][j] -= 1
            if check[read-1] == 0:
                check[read-1] = i+1
            else:
                point = len(adj_utgs[0])-1 #pointer at the last position
                wasIn = False #check if edge was already in

                #while the first el in adj list is the one we are in, look for the same edge:
                # if found: add weight 
                # else: look down
                if point >= 0:
                    while adj_utgs[0][point] == i and point>=0:
                        if adj_utgs[1][point] == check[read-1]-1:
                            adj_utgs[2][point] += 1
                            wasIn = True
                            break
                        point -= 1
                
                if wasIn==False:
                    adj_utgs[0].append(i)
                    adj_utgs[1].append(check[read-1]-1)
                    adj_utgs[2].append(1)

    return unitigs,unitigs_gfa,adj_utgs, weigth_list

def read_gfa(gfafile):
    """
    args: path to gfa file
    returns:  1) utgs; 2) reads in each UTG; 3) adjacency list between utgs; 
    """
    f = open(gfafile,"r")
    gfas = f.read().split("S\tutg")
    f.close()
    del f
    del gfas[0]

    unitigs, reads_in, adj, weigth_list = elab_gfa(gfas)
    return unitigs, reads_in, adj, weigth_list

# UTGS HANDLING
def adjTaxas(taxas, newTaxa):#taxas = list of lists containing a taxa and the relative counter, newTaxa is the taxa evaluated
    pos = -1
    for i in range(len(taxas)):
        if taxas[i][0] == newTaxa:
            taxas[i][1] += 1
            pos = i
            break
    if pos == -1: #if the taxa is new
        taxas.extend([[newTaxa,1]])
    #I will try to sort the list from the taxa with higher hits to the one with lower in order to minimize access.
    else: # check if I need to reoder the list.
        #taxas.sort( key=lambda x : x[1], reverse=True)
        for i in range(pos):
            if taxas[pos-i][1] > taxas[pos-i-1][1]:
                temp = taxas[pos-i-1]
                taxas[pos-i-1] = taxas[pos-i]
                taxas[pos-i] = temp

def get_precision(listoflists,taxa,minel,lev):
    result = []
    for i in range(len(listoflists)):
        if len(listoflists[i])>minel:
            temp = []

            for el in listoflists[i]:
                adjTaxas(temp,int(taxa[int(el)][lev]))
            result.append([temp[0][0],temp[0][1]/len(listoflists[i]),len(listoflists[i])])
    return result

def get_metrics(listoflists,taxa,minel,lev):
    result = []
    for i in range(len(listoflists)):
        result_temp = []
        if len(listoflists[i])>minel:
            temp = []

            for el in listoflists[i]:
                adjTaxas(temp,int(taxa[int(el)][lev]))
            for j in range(len(temp)):
                result_temp.append([temp[j][0],temp[j][1],len(listoflists[i])])
        result.append(result_temp)
    return result

def get_prec_recall_f1(metrics,n_species,n_unclustered,offset):   
    # precision
    prec_reads_up = 0
    prec_reads_down = 0
    for clus in metrics:
        #clus[0] è il gruppo di reads che rappresenta la specie dominante nel cluster; 
        # clus[0][1] mi da il # di reads della specie dominante; 
        # clus[0][1] mi da il # di reads nel cluster
        prec_reads_up += int(clus[0][1])
        prec_reads_down += int(clus[0][2])
    
    #recall 
    rec_reads_up = 0
    rec_reads_down = 0
    actual_reads = np.zeros(n_species)
    for clus in metrics:
        for specie in clus:
            if int(specie[1]) > actual_reads[int(specie[0])-offset]:
                actual_reads[int(specie[0])-offset] = int(specie[1])
    rec_reads_up = actual_reads.sum()
    rec_reads_down = prec_reads_down + n_unclustered

    precision = prec_reads_up/prec_reads_down

    recall = rec_reads_up/rec_reads_down

    f1 = (2*precision*recall)/(precision+recall)

    return round(precision,3),round(recall,3),round(f1,3)


def get_prec(metrics,n_species,n_unclustered,offset):   
    # precision
    prec_reads_up = 0
    prec_reads_down = 0
    for clus in metrics:
        #clus[0] è il gruppo di reads che rappresenta la specie dominante nel cluster; 
        # clus[0][1] mi da il # di reads della specie dominante; 
        # clus[0][1] mi da il # di reads nel cluster
        prec_reads_up += int(clus[0][1])
        prec_reads_down += int(clus[0][2])

    precision = prec_reads_up/prec_reads_down

    return round(precision,3)




def mean_prec(precisions):
    prec_up,prec_down = 0,0
    for i in range(len(precisions)):
        prec_up += precisions[i][1]*precisions[i][2]
        prec_down += precisions[i][2]
    return prec_up/prec_down

def get_unique_elements(grps):
    bucket = []
    for grp in grps:
        bucket += grp
    bucket = list(set(bucket))
    return bucket


def get_label_list(result):
    label = []
    for el in result:
        label.append(el[0])
    return label

def get_edge_prec(edg1,edg2,label):
    edg = []
    for i in range(len(edg1)):
        if label[edg1[i]] == label[edg2[i]]:
            edg.append(1)
        else: 
            edg.append(0)
    return edg

def get_div(str1):
    ext = get_ext(str1)
    if ext == 'fa':
        return 2
    else:
        return 4

# filename extension
def get_ext(str1):
    if str1[-2:] == 'fa':
        return 'fa'
    elif str1[-2:] == 'fq':
        return 'fq'
    else:
        sys.exit("Errore nel formato di ingresso dei file. Deve essere .fa o .fq")


def get_utg_clus_prec(clusters,labels):
    result = []
    for clus in clusters:
        temp = []
        for el in clus:
            adjTaxas(temp,labels[el])
        result.append([temp[0][0],temp[0][1]/len(clus),len(clus)])
    return result

def get_grpd_utgs(utgs,deg,clusters,maxLen):
    """
    args: unitigs sequences list, degree array, list of list of clusters, maximum length of grouped unitigs
    returns: list of grouped unitigs based on untig degree in the network 
    """
    aggr_utgs = []
    for clus in clusters:
        #creo la lista con i degrees di ogni utg nel gruppo.

        curr_deg = np.zeros(len(clus),dtype=np.int32)
        for i in range(len(clus)):
            curr_deg[i] =  deg[clus[i]]
        
        #ora ciclo while su lunghezz utg finale
        totLen = 0
        big_utg = ""
        while totLen < maxLen:

            #selezione argmax
            x = np.argmax(curr_deg) #trovo argmax nell'array
            if curr_deg[x] == -1: #se il max è -1 vuol dire che ho messo tutti gli elementi del cluster
                break
            current = clus[x] # mi salvo che utg è
            curr_deg[x] = -1 # pongo degree = -1 per non considerarlo più

            #creazione big utg
            if len(utgs[current]) < maxLen-totLen:
                big_utg += utgs[current] +'N'
                totLen += len(utgs[current])+1
            else:
                big_utg += utgs[current][:maxLen-totLen]+'N'
                totLen += len(utgs[current][:maxLen-totLen])+1
            
        aggr_utgs.append(big_utg[:-1])
    
    return aggr_utgs

def get_grpd_utgs_qual(grpd_prec,deg,clusters,maxLen,labe):
    """
    args: unitigs sequences list, degree array, list of list of clusters, maximum length of grouped unitigs
    returns: list of grouped unitigs based on untig degree in the network 
    """
    chosen_utgs = []
    for clus in clusters:
        #creo la lista con i degrees di ogni utg nel gruppo.

        curr_deg = np.zeros(len(clus),dtype=np.int32)
        for i in range(len(clus)):
            curr_deg[i] =  deg[clus[i]]
        
        #ora ciclo while su lunghezz utg finale
        totLen = 0
        chos = []
        while totLen < maxLen:

            #selezione argmax
            x = np.argmax(curr_deg) #trovo argmax nell'array
            if curr_deg[x] == -1: #se il max è -1 vuol dire che ho messo tutti gli elementi del cluster
                break
            current = clus[x] # mi salvo che utg è
            curr_deg[x] = -1 # pongo degree = -1 per non considerarlo più

            #creazione big utg
            chos.append(current)
        
        chosen_utgs.append(chos)
    i = 0
    chos_prec = []
    for grp in chosen_utgs:
        ch = 0
        for el in grp:
            if int(labe[el]) == int(grpd_prec[i][0]):
                ch +=1
        chos_prec.append(ch/len(grp))
        i += 1

    return chos_prec

def get_grpd_utgs_precise(grpd_prec,utgs,deg,clusters,maxLen,labe):
    """
    args: unitigs sequences list, degree array, list of list of clusters, maximum length of grouped unitigs
    returns: list of grouped unitigs based on untig degree in the network but checked
    """
    aggr_utgs = []
    k = 0
    for clus in clusters:
        #creo la lista con i degrees di ogni utg nel gruppo.

        curr_deg = np.zeros(len(clus),dtype=np.int32)
        for i in range(len(clus)):
            curr_deg[i] =  deg[clus[i]]
        
        #ora ciclo while su lunghezz utg finale
        totLen = 0
        big_utg = ""
        while totLen < maxLen:

            #selezione argmax
            x = np.argmax(curr_deg) #trovo argmax nell'array
            if curr_deg[x] == -1: #se il max è -1 vuol dire che ho messo tutti gli elementi del cluster
                break
            if int(labe[clus[x]]) == int(grpd_prec[k][0]):

                current = clus[x] # mi salvo che utg è
                curr_deg[x] = -1 # pongo degree = -1 per non considerarlo più

                #creazione big utg
                if len(utgs[current]) < maxLen-totLen:
                    big_utg += utgs[current] +'N'
                    totLen += len(utgs[current])+1
                else:
                    big_utg += utgs[current][:maxLen-totLen]+'N'
                    totLen += len(utgs[current][:maxLen-totLen])+1
            else:
                curr_deg[x] = -1 # pongo degree = -1 per non considerarlo più
        k += 1
        aggr_utgs.append(big_utg[:-1])
    return aggr_utgs



def get_outputcsv(clusters,grps):
    """
    args: clusters of utgs, reads in utgs
    returns: list of lists, in each list there are all the reads in one cluster
    """
    csv_out = []
    for clus in clusters:
        one_cl = []
        for el in clus:
            one_cl.extend(grps[el])
        csv_out.append(list(set(one_cl)))

    return csv_out

def get_range(utgrange,clust,utgs):
    popl = []
    for i in range(len(utgs)):
        utg = utgs[i]
        if len(utg) > utgrange[0] and len(utg) < utgrange[1]:
            popl.append(i)
    new_clust = []
    new_utgs = []
    for el in popl:
        new_clust.append(clust[el])
        new_utgs.append(utgs[el])
    return new_utgs,new_clust
    

def get_leftOut(grps,n_seq):
    all_seq = np.zeros(n_seq,dtype=np.uint8)
    for grp in grps:
        for el in grp:
            all_seq[el] = 1
    
    lef_out = []
    for i in range(len(all_seq)):
        if all_seq[i] == 0:
            lef_out.append(i)

    return lef_out
    
def get_clusters(elems):
    tot_c = max(elems)
    final_c = []
    left_out = 0
    for i in range(tot_c+1):
        final_c.append([])

    for i in range(len(elems)):
        if elems[i] != -1:
            final_c[elems[i]].append(i)
        else:
            left_out += 1
    return final_c,left_out

### FASTA FILE CREATION ###

def add_grps(grpd):
    grps = ''
    i = 1
    for el in grpd:
        first_line ='>G'+str(i)+'\n'
        grps +=first_line+el+'\n'
        i+=1
    return grps

def add_utgs(grpd,maxlen):
    grps = ''
    i = 1
    for el in grpd:
        first_line ='>G'+str(i)+'\n'
        if len(el) > maxlen:
            grps +=first_line+el[:maxlen]+'\n'
        else:
            grps +=first_line+el+'\n'
        i+=1
    return grps


def add_reads(left_out,reads_fq,div,leng):
    lOut = ''
    for i in range(len(left_out)):    
        first_line = '>'+reads_fq[left_out[i]*div][1:]+'\n'
        lOut += first_line+reads_fq[left_out[i]*div+1]+'N'+reads_fq[(left_out[i]+leng)*div+1]+'\n'

    return lOut

def add_reads_single(left_out,reads_fq,div):
    lOut = ''
    for i in range(len(left_out)):    
        first_line = '>'+reads_fq[left_out[i]*div][1:]+'\n'
        lOut += first_line+reads_fq[left_out[i]*div+1]+'\n'

    return lOut


##### PAF FILE HANDLING #######

def process_frame_list(daf):
    daf = daf.loc[:,[0,5,9]]
    quer = daf.loc[:,0].map(lambda x: int(x[1:-2])-1).tolist()
    del daf[0]
    tar = daf.loc[:,5].map(lambda x: int(x[1:-2])-1).tolist()
    del daf[5]
    wei = daf.loc[:,9].astype(np.uint8).tolist()
    del daf
    return [quer,tar,wei]

def openpaflist(filename,size,n_jobs):
    ctx = mp.get_context('spawn')
    pool = ctx.Pool(n_jobs)
    data = pool.imap(process_frame_list,pd.read_csv(filename,chunksize=size,sep='\t',header=None))
    pool.close()
    pool.join()
    return data

def get_paf(paf,dim,n_cores):
    data = openpaflist(paf,dim,n_cores)
    #print('CSV MULTIREADER',time.time()-t1)
    df = [[],[],[]]
    
    for d in data:
        #print(type(d))
        df[0].extend(d[0])
        df[1].extend(d[1])
        df[2].extend(d[2])
    
    del data
    return df

def create_dic(utgs):
    d = {}
    for i in range(len(utgs)):
        for j in utgs[i]:
            d[j] = i
    return d


def consensus_utgs(dic_utgs,clusters):
    consensus = []
    for clus in clusters:
        temp = []
        for el in clus:
            try:
                x = int(dic_utgs[el])
                adjTaxas(temp, x)
            except:
                continue
        try:
            consensus.append(int(temp[0][0]))
        except:
            consensus.append(-1)
    return consensus





if __name__ == '__main__':

    gfa = '/Users/francescoandreace/Downloads/SRR1804065_m70.gfa'
    utgs, grps, adj = read_gfa(gfa)
    print(len(utgs))
    print(utgs[0])
    print(len(utgs[0]))
