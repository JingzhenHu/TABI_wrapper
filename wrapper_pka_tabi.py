# Works for Python 2.7.9
# Command line: python wrapper_pKa PDBID
# I keep some print comment in case intermediate check
import sys
import os
import re
import requests
import math
import json
from copy import deepcopy
# import matplotlib.pyplot as plt

# this is the struture for store a line in PDB
class PDB:
    def __init__(self,fix,index,atom,aminoAcid,chain,aaNumber,x,y,z,c,r,atom_nospecific):
        self.fix = fix
        self.index = index
        self.atom = atom
        self.aminoAcid = aminoAcid
        self.chain = chain
        self.aaNumber = aaNumber
        self.x = x
        self.y = y
        self.z = z
        self.c = c
        self.r = r
        self.atom_nospecific = atom_nospecific

# this is the struture for store a line in PQR
class PQR:
    def __init__(self,fix,index,atom,aminoAcid,aaNumber,x,y,z,charge,radius):
        self.fix = fix
        self.index = index
        self.atom = atom
        self.aminoAcid = aminoAcid
        self.aaNumber = aaNumber
        self.x = x
        self.y = y
        self.z = z
        self.charge = charge
        self.radius = radius

# Given specific PDBID, read the online and split into three parts:
# 1. header - the content before the first line beginning with 'ATOM...'
# 2. pdb_list - all the lines beginning with 'ATOM...', spliting each line into the PDB data structure
# 3. trailor - the content after the last line beginning with 'ATOM...'
def downloadPDB( PDBID ):
    pdb_list = []
    urlName = ''.join(['https://files.rcsb.org/view/',PDBID,'.pdb'])
    page = requests.get(urlName).text.split('\n')
    line_index = 0
    atom_line = []
    for line in page:
        if re.match(r'^(ATOM).+$', line):
            line = line.split()
            pdb_list.append(
                PDB(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11]))
            atom_line.append(line_index)
        line_index += 1
    header = []
    trailor = []
    for index in range(0,len(page)):
        if index < atom_line[0]:
            header.append(page[index])
        if index > atom_line[len(atom_line)-1]:
            trailor.append(page[index])
    return pdb_list,header,trailor

# Located in the path, write out the PDB file based on the passing in list of PDB data structure
def writePDB(path, datalist, header, trailor):
    if os.path.isfile(path):
        return
    else:
        newpdb = open(path,'w')
        for line in header:
            newpdb.write(line)
            newpdb.write('\n')
        for obj in datalist:
            newpdb.write('{:<6}{:>5}  {:<4}'.format(obj.fix,obj.index,obj.atom))
            newpdb.write('{:<3} {:<1}{:>4}    '.format(obj.aminoAcid,obj.chain,obj.aaNumber))
            newpdb.write(' {:>7,.3f} {:>7,.3f} {:>7,.3f} '.format(float(obj.x),float(obj.y),float(obj.z)))
            newpdb.write(' {:<5,.2f}{:>5,.2f}'.format(float(obj.c),float(obj.r)))
            newpdb.write('           {:<1}'.format(obj.atom_nospecific))
            newpdb.write('\n')
        for line in trailor:
            newpdb.write(line)
            newpdb.write('\n')
        newpdb.close()
        return

# Located in the path, write out the PQR file based on the passing in list of PQR data structure
def writePQR(path, datalist):
    if os.path.isfile(path):
        return
    else:
        newpqr = open(path,'w')
        for obj in datalist:
            newpqr.write('{:<8}{:<8}{:<9}'.format(obj.fix,obj.index,obj.atom))
            newpqr.write('{:<8}{:<5}'.format(obj.aminoAcid,obj.aaNumber))
            newpqr.write('{:<10,.5f}{:<10,.5f}{:<10,.5f}'.format(float(obj.x),float(obj.y),float(obj.z)))
            newpqr.write('{:<9,.5f}'.format(float(obj.charge)))
            newpqr.write('{:<8,.5f}\n'.format(float(obj.radius)))
        newpqr.close()
        return

# change the activated site to its on-pronated or off-pronated name according to the parameter
# para: res-RSDID(can be on or off); status-('on' or 'off')
def titrResDict(res,status):
    titr_res_dict_on = {'AR0':'ARG','ASP':'ASH','GLU':'GLH','HIE':'HIP','HIS':'HIP','LYN':'LYS','TYM':'TYR'} # 'CYX':'CYM',
    titr_res_dict_off = {'ARG':'AR0','ASH':'ASP','GLH':'GLU','HIP':'HIE','HIS':'HIE','LYS':'LYN','TYR':'TYM'}
    for site in titr_res_dict_on:
        if (res==site):
            if (status == 'on'):
                return titr_res_dict_on[site]
            if (status == 'off'):
                return res
    for site in titr_res_dict_off:
        if (res==site):
            if (status == 'off'):
                return titr_res_dict_off[site]
            if (status == 'on'):
                return res

# give the contents map with necessary information, write or revise the usrdata.in in the same dir
def writeUseData(contents):
    # for mibpb3,
    # par = ['fname','dcel','den', 'epsp','epsw', 'bulk_strength', 'icg', 'isf','ibd','irg','inl']#,'twob_pos']
    # for tabipb
    par = ['fname', 'den', 'epsp', 'epsw', 'bulk_strength', 'order', 'maxparnode', 'mac']
    usrda = open("usrdata.in",'w')
    # for mibpb3 it should be range(11)
    # for tabi it should be range(8)
    for index in range(8):
        line = "".join([par[index], " ", str(contents[par[index]]), "\n"])
        usrda.write(line)
    # usrda.write(''.join([par[11],' ',contents[par[11]][0],' ',contents[par[11]][1]]))
    usrda.close()

# give the output filename aftering running the solver
# search whether the filename is completed or not and then grasp and return the free energy
def fetchEnergy(name):
    result = open(name,'r')
    endJudge = 'Total'
    search = '-----Free energy=:'
    running = True
    while running:
        if endJudge in result.read():
            running = False
    result.seek(0)
    energyStr = ''
    for line in result:
       if search in line:
           # print(line)
           energyStr = line[20:39]
           break
    result.close()
    return float(energyStr)

# combination all the steps from generating pqr from pdb, running solver and grasp the free energy
# via calling the functions above
def run_script(out_path,name,usrdata_map):
    path = ''.join(['test_proteins/', name,'.pqr'])
    pdb2pqr = ''.join(['python ','pdb2pqr/pdb2pqr.py ','--ff=PARSE ',''.join([out_path,name,'.pdb ']),path])
    #call pdb2pqr file to generate pqr file
    os.system(pdb2pqr)
    out_path = ''.join([out_path, name, '_out'])
    usrdata_map['fname'] = name
    writeUseData(usrdata_map)
    # tabipb.exe
    os.system(''.join(['./tabipb.exe', ' > ', out_path]))
    # mibpb3.exe
    # os.system(''.join(['./mibpb3.exe', ' > ', out_path]))
    return fetchEnergy(out_path)

# combination all the steps from writing pqr based on datalist, running solver and grasp the free energy
# via calling the functions above
def run_with_pqr(protonated,number,pqr_list,outpath,usrdata_map):
    name = ''.join(['aalone_',PBDID,'_',protonated,number])
    writePQR(''.join(['test_proteins/',name,'.pqr']),pqr_list)
    out_path = ''.join([outpath, name, '_out'])
    # if (os.path.isfile(path) == False):
    usrdata_map['fname'] = name
    writeUseData(usrdata_map)
    # tabipb.exe
    os.system(''.join(['./tabipb.exe', ' > ', out_path]))
    # mibpb3.exe
    # os.system(''.join(['./mibpb3.exe', ' > ', out_path]))
    return fetchEnergy(out_path)

# Given certain titratable site, compute the intrinsic_pKa
def intristicPKA(PBDID,outpath,E_all_mute,mute_list,AA,number,header,trailor,usrdata_map):
    deal_list = deepcopy(mute_list)
    protonated = titrResDict(AA,'on')
    depronated = titrResDict(AA,'off')
    model_pKa = {'ARG':12.0,'ASP':4.0,'CYS':9.5,'GLU':4.4,'HIS':6.3,'LYS':10.4,'TYR':9.6}
    for atom in deal_list:
        if (atom.aminoAcid == depronated and atom.aaNumber == number):
            atom.aminoAcid = protonated
    writePDB(''.join([outpath,PBDID,'_',protonated,number,'.pdb']),deal_list,header,trailor)
    # Free energy of protonated one titritable site on protein
    E_p_H =run_script(outpath,''.join([PBDID,'_',protonated,number]),usrdata_map)
    aa_alone = open(''.join(['test_proteins/',''.join([PBDID,'_',protonated,number]),'.pqr']),'r')
    pqr_list = []
    for line in aa_alone:
        if re.match(r'^(ATOM).+$', line):
            line = line.split()
            if line[3] == protonated and line[4] == number:
                pqr_list.append(
                PQR(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9]))
    # Free energy of protonated one titritable site off protein
    E_s_H = run_with_pqr(protonated,number,pqr_list,outpath,usrdata_map)
    pqr_all_mute = open(''.join(['test_proteins/',''.join([PBDID,'_all_mute']),'.pqr']),'r')
    pqr_list_un = []
    for line in pqr_all_mute:
        if re.match(r'^(ATOM).+$', line):
            line = line.split()
            if line[3] == depronated and line[4] == number:
                pqr_list_un.append(
                PQR(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9]))
    # Free energy of deprotonated one titritable site off protein
    E_s = run_with_pqr(depronated,number,pqr_list_un,outpath,usrdata_map)
    if AA in model_pKa:
        pKa_0 = model_pKa[AA]
        return pKa_0 - ((E_p_H - E_all_mute) - (E_s_H - E_s))*4.182/(2.5*2.303) # whether need to mult by 4.182
    if protonated in model_pKa:
        pKa_0 = model_pKa[protonated]
        return pKa_0 - ((E_p_H - E_all_mute) - (E_s_H - E_s))*4.182/(2.5*2.303) # whether need to mult by 4.182

# call the intristicPKA function the number of titratable site's times to get the list of intrinsic_pKa
# if we have computed the list already and stored them, then just read the data rather than recomputing them
def run_intristic(PBDID,outpath,titr_site,E_all_mute,pdb_all_mute,header,trailor,usrdata_map):
    intristic_file = ''.join([outpath,'intristic_pKa'])
    intristic_pKa = {}
    if (os.path.isfile(intristic_file) == False):
        intristic_pKa = {}
        for site in titr_site:
            for i in titr_site[site]:
                intristic_pKa[''.join([site, i])] = intristicPKA(PBDID,outpath, E_all_mute, pdb_all_mute, site, i, header,
                                                             trailor, usrdata_map)
        intris = open(''.join([outpath,'intristic_pKa']),'w')
        for pKa in intristic_pKa:
            intris.write('{:<5} {:<8,.3f}'.format(pKa,intristic_pKa[pKa]))
            intris.write('\n')
        intris.close()
    else:
        intri = open(intristic_file,'r')
        for line in intri:
            value = line.split()
            print(value[0],value[1])
            intristic_pKa[value[0]] = value[1]
    return intristic_pKa

# Given certain titratable site pairs s1 and s2, prepare the PQR file and compute the corresponding site-site Energy
def site_site_interaction(s1,s2,PBDID,usrdata_map,out_path):
    pqr_mute_name = ''.join(['test_proteins/', ''.join([PBDID, '_all_mute']), '.pqr'])
    pqr_all_mute = open(pqr_mute_name, 'r')
    s1_aa = s1[0:3]
    s1_num = s1[3:]
    s2_aa = s2[0:3]
    s2_num = s2[3:]
    s1_protonated = titrResDict(s1_aa,'on')
    s1_unprotonated = titrResDict(s1_aa, 'off')
    s2_protonated = titrResDict(s2_aa, 'on')
    s2_unprotonated = titrResDict(s2_aa, 'off')
    pqr_s1 = []
    s1_portion = open(''.join(['test_proteins/','aalone_',PBDID,'_',s1_protonated,s1_num,'.pqr']),'r')
    for line in s1_portion:
        line = line.split()
        pqr_s1.append(PQR(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9]))
    pqr_s2 = []
    s2_portion = open(''.join(['test_proteins/', 'aalone_', PBDID, '_', s2_protonated, s2_num,'.pqr']), 'r')
    for line in s2_portion:
        line = line.split()
        pqr_s2.append(PQR(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9]))
    pqr_site_site = []
    index_s1 = 0
    index_s2 = 0
    for line in pqr_all_mute:
        if re.match(r'^(ATOM).+$', line):
            line = line.split()
            if line[3] == s1_unprotonated and line[4] == s1_num and index_s1 == 0:
                for atom in pqr_s1:
                    pqr_site_site.append(atom)
                index_s1 += 1
            elif line[3] == s2_unprotonated and line[4] == s2_num and index_s2 == 0:
                for atom in pqr_s2:
                    pqr_site_site.append(atom)
                index_s2 += 1
            if line[3] != s1_unprotonated or line[4] != s1_num:
                if line[3] != s2_unprotonated or line[4] != s2_num:
                    pqr_site_site.append(PQR(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], 0.0, line[9]))
    pqr_all_mute.close()

    # all_mute with s1 and s2 on
    name = ''.join([PBDID,'_',s1,'_',s2])
    writePQR(''.join(['test_proteins/',name,'.pqr']),pqr_site_site)
    usrdata_map['fname'] = name
    writeUseData(usrdata_map)
    out_path_ = ''.join([out_path,name,'_out'])
    #if (os.path.isfile(out_path_) == False):
    os.system(''.join(['./tabipb.exe', ' > ', out_path_]))
    # os.system(''.join(['./mibpb3.exe', ' > ', out_path_]))

    # all_mute with s1 on only
    name_1 = ''.join([PBDID,'_',s1,'_s_s'])
    if (os.path.isfile(name_1) == False):
        pqr_all_mute = open(pqr_mute_name, 'r')
        pqr_site_1 = []
        index_s1 = 0
        for line in pqr_all_mute:
            if re.match(r'^(ATOM).+$', line):
                line = line.split()
                if line[3] == s1_unprotonated and line[4] == s1_num and index_s1 == 0:
                    for atom in pqr_s1:
                        pqr_site_1.append(atom)
                    index_s1 += 1
                if line[3] != s1_unprotonated or line[4] != s1_num:
                    pqr_site_1.append(PQR(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], 0.0, line[9]))
        pqr_all_mute.close()

        writePQR(''.join(['test_proteins/', name_1, '.pqr']), pqr_site_1)
        usrdata_map['fname'] = name_1
        writeUseData(usrdata_map)
        out_path_1 = ''.join([out_path,name_1,'_out'])
        #if (os.path.isfile(out_path_1) == False):
        os.system(''.join(['./tabipb.exe', ' > ', out_path_1]))
        # os.system(''.join(['./mibpb3.exe', ' > ', out_path_1]))

    # all_mute with s2 on only
    name_2 = ''.join([PBDID,'_',s2,'_s_s'])
    if (os.path.isfile(name_2) == False):
        pqr_all_mute = open(pqr_mute_name, 'r')
        pqr_site_2 = []
        index_s2 = 0
        for line in pqr_all_mute:
            if re.match(r'^(ATOM).+$', line):
                line = line.split()
                if line[3] == s2_unprotonated and line[4] == s2_num and index_s2 == 0:
                    for atom in pqr_s2:
                        pqr_site_2.append(atom)
                    index_s2 += 1
                if line[3] != s2_unprotonated or line[4] != s2_num:
                    pqr_site_2.append(PQR(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], 0.0, line[9]))
        pqr_all_mute.close()

        writePQR(''.join(['test_proteins/', name_2, '.pqr']), pqr_site_2)
        usrdata_map['fname'] = name_2
        writeUseData(usrdata_map)
        out_path_2 = ''.join([out_path,name_2,'_out'])
        #if (os.path.isfile(out_path_2) == False):
        os.system(''.join(['./tabipb.exe', ' > ', out_path_2]))
        # os.system(''.join(['./mibpb3.exe', ' > ', out_path_2]))

    return fetchEnergy(out_path_) - fetchEnergy(out_path_1) - fetchEnergy(out_path_2)

# loop and call the site_site_interaction, store them in a matrix
# if the data has existed, just read them in
def interaction_Energy(n,site_site,PBDID,usrdata_map,out_path):
    Matrix = [[0 for x in range(n)] for y in range(n)]
    name = ''.join([out_path,'site_site'])
    if (os.path.isfile(name) == False):
        s_s_interact = open(name,'w')
        for i in range(0,n):
            for j in range(0,n):
                if (i == j):
                    Matrix[i][j] = None
                # only calculate one-half site-site without dividing one-half later in computing statue-energy
                if (i < j):
                    Matrix[i][j] = site_site_interaction(site_site[i],site_site[j],PBDID,usrdata_map,out_path)
                    print(site_site[i],site_site[j],Matrix[i][j])
                    s_s_interact.write('{:<5} {:<5} {:<8,.3f}'.format(site_site[i],site_site[j],Matrix[i][j]))
                    s_s_interact.write('\n')
        s_s_interact.close()
    else:
        s_s_interact = open(name, 'r')
        content = s_s_interact.read().split('\n')
        amount = len(content)
        for i in range(0,amount-1):
            com = content[i].split()
            for i in range(0,n):
                for j in range(0,n):
                    if com[0] == site_site[i] and com[1] == site_site[j]:
                        Matrix[i][j] = float(com[2].replace(',','')) # convert the readin without comma into float type
        s_s_interact.close()
    return Matrix

# the function computes all the energy transfers (total number: 2^(the # of titratable sites))
# and store them in a map with keys are in a binary number in which 1 in ith position represents the
# ith titratable site is porotonated and 0 means unprotonated.
# the order of titratable sites is fixed as the site_site list generated before computing site_site.
def compute_statue_energy(PH, outpath, site_site, active_titr_len, intristic_pKa, E_interact):
    status_Energy = {}
    name = ''.join([outpath, 'Status_Energy_', str(PH)])
    # if (os.path.isfile(name) == False):
    for i in range(2 ** active_titr_len):
        energy = 0.0
        com = '{0:b}'.format(i)
        diff = active_titr_len - len(com)
        if diff != 0:
            add = ''
            for i in range(0, diff):
                add += '0'
            com = add + com
        # print(com)
        index = 0
        coeff = -2.5 * 2.303 / 4.182
        for each in com:
            if each == '1':
                site = site_site[index]
                pKa_0 = intristic_pKa[site]
                energy += (float(pKa_0) - PH)
            index += 1
        energy *= coeff
        E_site = 0
        site_1_index = 0
        for each in com:
            if each == '1':
                site_2_index = 0
                for second in com:
                    if second == '1' and site_2_index != site_1_index:
                        E_site += E_interact[site_1_index][site_2_index]
                    site_2_index += 1
            site_1_index += 1
        # E_site = E_site / 2
        energy += E_site
        status_Energy[com] = energy
    '''
    code for writing out Statue Energy if possible to check 
        Status_E = open(name, 'w')
        for key in status_Energy:
            Status_E.write(key)
            Status_E.write(' ')
            Status_E.write(str(status_Energy[key]))
            Status_E.write('\n')
    else:
        Status_E = open(name, 'r')
        content = Status_E.read().split('\n')
        # print('length: ', len(content))
        for line in content:
            if line != '':
                two = line.split()
                status_Energy[two[0]] = float(two[1])
    '''
    return status_Energy

# compute the numerator of the last equation in paper
# substitution with the invariable denominator and judge the answer relative to 1/2 with err_tol(error tolerance)
# if the result is satisfied with the standard, then keep it.
def is_pKa_okay(status_Energy, denominator, active_titr_len, site_site, pKa, PH, err_tol,table):
    for i in range(0,active_titr_len):
        # pKa[site_site[i]] = 0.0
        energy = 0
        for key in status_Energy:
            if key[i] == '1':
                energy += math.exp(-status_Energy[key]*4.182/2.5)
        table.append(energy/denominator)
    return

def main(PBDID):
    pdb_list, header, trailor = downloadPDB(PBDID)
    titr_res_name = ['ARG','ASP','GLU','HIS','LYS','TYR'] # 'CYS'
    # make a folder containing downloaded pdb files, transformed pqr file, energy outputs
    mkdir = ''.join(['mkdir -p ', PBDID])
    # os.system('ulimit -s 65532')
    os.system(mkdir)
    # read usrdata file into map
    usrdata = open('usrdata.in','r').read().split()
    usrdata[1] = PBDID
    usrdata_map = {}
    # for mibpb3 should be range(11)
    # for tabpib should be range(8)
    for i in range(8):
        usrdata_map[usrdata[2*i]] = usrdata[2*i+1]
    # twob_pos = [usrdata[23],usrdata[24]]
    # usrdata_map[usrdata[22]] = twob_pos

    # collect all the titratable sites in each protein
    titr_site = {}
    for titr in titr_res_name:
        # the reason to use set is that the RSDID is corresponding to atom which shows repeatedly
        titr_site[titr] = set()
    for atom in pdb_list:
        for titr in titr_res_name:
            if titr == atom.aminoAcid:
                titr_site[titr].add(atom.aaNumber)
    # change the map's value(set) into a list, then we can iterate item in list with order
    for res in titr_site:
        to_list = []
        for ele in titr_site[res]:
            to_list.append(ele)
        titr_site[res] = to_list
    print(titr_site)

    # get the all mute pdb (turn off all titratable sites), as we search for those site before,
    # I only rename the RSDID as their unprotonated status's name here and then get the free energy
    pdb_all_mute = pdb_list
    for site in titr_site:
        for i in titr_site[site]:
            for line in pdb_all_mute:
                if line.aaNumber == i:
                    line.aminoAcid = titrResDict(site,'off')
    # for obj in pdb_all_mute:
    #     print(obj.aaNumber, obj.aminoAcid)
    outpath = ''.join([PBDID, '/'])
    writePDB(''.join([outpath,PBDID,'_all_mute','.pdb']),pdb_all_mute,header,trailor)
    E_all_mute = run_script(outpath,''.join([PBDID,'_all_mute']),usrdata_map)
    # print('Energy_all_mute: ', E_all_mute)

    intristic_pKa = run_intristic(PBDID,outpath,titr_site,E_all_mute,pdb_all_mute,header,trailor,usrdata_map)
    # site-site interaction energies prepare
    site_site = []
    for site in titr_site:
        for i in titr_site[site]:
            site_site.append(''.join([site,i]))
    active_titr_len = len(site_site)
    E_interact = interaction_Energy(active_titr_len,site_site,PBDID,usrdata_map, outpath)
    for ele in site_site:
        print(ele)
    # dimension = len(site_site)
    # for i in range(0,dimension):
    #     for j in range(0,dimension):
    #         print(E_interact[i][j])
    res = []
    pKa = {}
    PH = 0.8 # hard code now, later changing to a dense dataset and loop the following steps
    # PH_arr = []
    obs_status_eng = open(''.join([outpath,'status_eng']),'w')
    while PH < 20:
        PH += 0.2
        # PH_arr.append(PH)
        status_Energy = compute_statue_energy(PH, outpath, site_site, active_titr_len, intristic_pKa, E_interact)
        obs_status_eng.write(str(PH))
        obs_status_eng.write('\n')
        table = []
        denominator = 0
        count1 = 0
        count2 = 0
        for key in status_Energy:
            jud = math.exp(-status_Energy[key]*4.182/2.5)
            if (jud < 1e-15 ):
                count1 += 1
            elif (jud >= 1e-15 and jud < 1e-5):
                count2 += 1
            denominator += math.exp(-status_Energy[key]*4.182/2.5)
        # print('denominator: ', denominator) # the number is out of range of double storage
        err_tol = 0.05
        is_pKa_okay(status_Energy, denominator, active_titr_len, site_site, pKa, PH, err_tol,table)
        res.append(table)
    # print('pKa_satisfied: \n')
    # print(json.dumps(pKa))

    # print('table for each pH: \n')
    # print(json.dumps(res))
    # print('\n')

    rerange = {}
    for item in site_site:
        rerange[item] = []
    for each in res:
        index = 0
        for item in each:
            rerange[site_site[index]].append(item)
            index += 1

    print('table for each item: \n')
    print(json.dumps(rerange))

    result = open(''.join([outpath,'RES.out']),'w')
    index = 0
    for item in rerange:
        result.write(str(item))
        index += 1
        result.write('\n')
        for each in rerange[item]:
            result.write(str(each))
            result.write('\n')
    result.close()

if __name__ == "__main__":
    PBDID = sys.argv[1]
    main(PBDID)
