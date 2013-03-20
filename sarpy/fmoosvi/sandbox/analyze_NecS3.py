# -*- coding: utf-8 -*-

## NecS3 Analysis File

# Imports

import sarpy
import numpy
from sarpy.fmoosvi import wrappers
import time

#NecS3_BSB1_FLASH = sarpy.Experiment('NecS3').find_scan_by_protocol('07_bSB1mapFLASH')
#NecS3_BSB1_MSME = sarpy.Experiment('NecS3').find_scan_by_protocol('BSB1map\-MSME_bas')
NecS3_exp = sarpy.Experiment('NecS3')
NecS3_LLscans = NecS3_exp.find_scan_by_protocol('04_ubcLL+')

start_time = time.time()

for scan in NecS3_LLscans:
    
#    BSB1map_MSME = wrappers.calculate_BSB1map(NecS3_exp, protocol_name = 'BSB1map\-MSME_bas')
#    BSB1map_FLASH = wrappers.calculate_BSB1map(NecS3_exp, protocol_name = '07_bSB1mapFLASH')
    T1map_LL = wrappers.calculate_T1map(NecS3_LLscans, protocol_name = '04_ubcLL+') 
    
#    T1map_LL_MSME = wrappers.calculateT1map(NecS3_exp, protocol_name = '04_ubcLL+', flip_angle_map = BSB1map_MSME)
#    T1map_LL_FLASH = wrappers.calculateT1map(NecS3_exp, protocol_name = '04_ubcLL+', flip_angle_map = BSB1map_FLASH)

#    scan.store_adata(key='BSB1map_MSME', data = BSB1map_MSME)
#    scan.store_adata(key='BSB1map_FLASH', data = BSB1map_FLASH)
#    scan.store_adata(key='T1map_LL', data = T1map_LL)
#    scan.store_adata(key='T1map_LL_MSME', data = T1map_LL_MSME)
#    scan.store_adata(key='T1map_LL_FLASH', data= T1map_LL_FLASH)        
    
end_time = time.time()
    
print 'This run took {0} seconds.'.format(round(end_time - start_time))
