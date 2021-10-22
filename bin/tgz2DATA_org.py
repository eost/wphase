#!/usr/bin/env python

import sacpy
import tarfile
from os       import walk, unlink, makedirs, chdir
from os.path  import join, exists, isdir, abspath, basename
from sys      import argv, stderr, exit
from shutil   import copy, rmtree
from tempfile import mkdtemp

def RM(name):
    if exists(name):
        if isdir(name): rmtree(name)
        else:           unlink(name)
    return

def populates_tmpdir(tgz):
    tgz     = abspath(tgz)
    tmpdir  = mkdtemp(dir='/tmp', prefix='_', suffix='.tmp')
    chdir(tmpdir)
    tar = tarfile.open(tgz)
    tar.extractall()
    tar.close()
    return abspath(tmpdir)

def adds_one_zero(pzin,pzout):
    #adds a zero and writes in Z-P-C order
    lipzc,liout = [],[]
    for line in open(pzin).readlines():
        if    len(line) == 0: continue
        elif  line[0] == '*': liout.append(line)
        else: lipzc.append(line)
    # Zeros
    for j in range(len(lipzc)):
        if lipzc[j][0:5] == 'ZEROS':
            n = int(lipzc[j].split()[1])
            liout.append('ZEROS %d\n0.000000e+00 0.000000e+00\n'%(n+1))
            for k in range(n): liout.append(lipzc[j+k+1])
            break
    # Poles
    for j in range(len(lipzc)):
        if lipzc[j][0:5] == 'POLES':
            n = int(lipzc[j].split()[1])
            for k in range(n+1): liout.append(lipzc[j+k])
            break
    # constant
    for j in range(len(lipzc)):
        if lipzc[j][0:8] == 'CONSTANT':
            liout.append(lipzc[j])
            break
    # Writes output file
    fd = open(pzout,'w')
    for line in liout: fd.write(line)
    fd.close()
    return

def avail_sacs_pzs_wilber(didi):
    chdir(didi)
    sacs,pzs = [],[]
    for root,dirs,files in walk(didi):
        for name in files:
            if name.endswith('.SAC') or name.endswith('.sac'):
                sacs.append(join(root,name))
            if name.startswith('SACPZ') or name.startswith('SAC_PZs'):
                pzs.append(join(root,name))
    return sacs, pzs

def avail_sacs_pzs_nied(didi):
    chdir(didi)
    sacs,pzs = [],[]
    for root,dirs,files in walk(didi):
        for name in files:
            if   name.endswith('.zp'): pzs.append(join(root,name))
            else:
                try:    s = sacpy.Sac(join(root,name))
                except: continue
                sacs.append(name)
    return sacs, pzs

def make_dicts(sacs, pzs, mode):
    sac_dict, start_dict, pz_dict, pz_pref = {}, {}, {}, {}
    if mode == 'wilber':
        for sac in sacs:
            s = sacpy.Sac(sac,datflag=False)
            if s.khole.strip() == '--' or s.khole.strip() == '-12345': s.khole = ''
            mykey = '%s.%s.%s.%s'%(s.knetwk,s.kstnm,s.khole,s.kcmpnm)
            sac_dict[mykey] = sac
            start_dict[mykey] = '%04d.%03d.%02d.%02d.%02d.%03d0'%(s.nzyear,s.nzjday,s.nzhour,s.nzmin,s.nzsec,s.nzmsec)
        for pz in pzs:
            n = 0
            lines = open(pz).readlines()
            for line in lines:
                if len(line) < 1: continue
                if line[0] != '*': continue
                segms = line[1:].strip().split(':')
                if len(segms) < 2: continue
                key, val = line[1:].strip().split(':')[:2]
                key = key.strip().split()[0]
                val = val.strip()
                if   key == 'NETWORK':
                    knetwk = val; n += 1
                elif key == 'STATION':
                    kstnm  = val; n += 1
                elif key == 'CHANNEL':
                    kcmpnm = val; n += 1
                elif key == 'LOCATION':
                    khole  = val; n += 1
                else: continue
                if n > 3:
                    pz_dict['%s.%s.%s.%s'%(knetwk,kstnm,khole,kcmpnm)] = pz
                    break
    elif mode == 'nied':
        knetwk = 'BO'
        khole  = ''
        for sac in sacs:
            s = sacpy.Sac(sac,datflag=False)
            mykey = '%s.%s.%s.%s'%(knetwk,s.kstnm,khole,s.kcmpnm)
            sac_dict[mykey] = sac
            start_dict[mykey] = '%04d.%03d.%02d.%02d.%02d.%03d0'%(s.nzyear,s.nzjday,s.nzhour,s.nzmin,s.nzsec,s.nzmsec)
        for pz in pzs:
            kstnm,kcmpnm = basename(pz).split('.')[0].split('_')
            pz_dict['%s.%s.%s.%s'%(knetwk,kstnm,khole,kcmpnm)] = pz
    return sac_dict, start_dict, pz_dict
    
def populates_destdir(sac_dict, start_dict, pz_dict, mode, destdir):
    # populating destdir with sacs and pzs
    pz_ids = list(pz_dict.keys())
    if exists(destdir): RM(destdir)
    makedirs(destdir)
    for mykey in sac_dict.keys():
        if mykey not in pz_ids:
            #stderr.write('No available PZ file for sac file %s\n'%sac_dict[mykey])
            continue
        sacname  = join(destdir,'%s.%s.X.SAC'%(start_dict[mykey],mykey))
        if mode == 'wilber':
            copy(sac_dict[mykey],sacname)
            knetwk,kstnm,khole,kcmpnm = mykey.split('.')
            if khole == '': khole = '--'
            pzname  = 'SAC_PZs_%s_%s_%s_%s_%s_%s'%(knetwk,kstnm,kcmpnm,khole,start_dict[mykey],start_dict[mykey])
            pzname  = join(destdir,pzname)
            copy(pz_dict[mykey], pzname)
        elif mode == 'nied':
            t = sacpy.Sac(sac_dict[mykey])
            t.knetwk = 'BO'
            t.khole  = '-12345'
            t.write(sacname)
            pzname  = 'SAC_PZs_%s_%s_%s_%s_%s_%s'%(t.knetwk,t.kstnm,t.kcmpnm,'--',start_dict[mykey],start_dict[mykey])
            pzname  = join(destdir,pzname)
            adds_one_zero(pz_dict[mykey], pzname)
        else:
            stderr.write('Only "wilber" and "nied" modes are coded\n')
            exit(1)
    return

def verif_syntax(argv):
    syntx_msg  = 'Syntax: %s tgz_file\n'%basename(argv[0])
    syntx_msg += '        %s tgz_file mode\n'%basename(argv[0])
    syntx_msg += '        %s tgz_file mode destdir\n'%basename(argv[0])
    syntx_msg += 'mode:     wilber/nied;  default wilber\n'
    syntx_msg += 'destdir:                default "DATA_org"\n'       
    destdir    = abspath('DATA_org')
    if   len(argv) == 2: mode    = 'wilber'
    elif len(argv) == 3: mode    = argv[2]
    elif len(argv) == 4:
        destdir = abspath(argv[3])
        mode    = argv[2]
    else:
        stderr.write(syntx_msg)
        exit(1)
    if mode != 'wilber' and mode != 'nied':
        stderr.write(syntx_msg)
        exit(1)
    return argv[1], mode, destdir

## Main
# This script takes a compressed tar ball either from "wilber"
#  (sac, little endian) or "nied" (sac and pzs)
# and creates a directory with sacs and pzs named with the rdseed conventions.

tgz, mode, destdir = verif_syntax(argv)
tmpdir   = populates_tmpdir(tgz)
if   mode == 'wilber': time_series_files,pzs = avail_sacs_pzs_wilber(tmpdir)
elif mode ==   'nied': time_series_files,pzs = avail_sacs_pzs_nied(  tmpdir)
sac_dict, start_dict, pz_dict   = make_dicts(time_series_files, pzs, mode)
populates_destdir(sac_dict, start_dict, pz_dict, mode, destdir)
RM(tmpdir)
