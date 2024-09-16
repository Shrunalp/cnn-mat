import numpy as np
from Reformat_DREAM3D_Data import b2w, np_to_block

sz = "128"

cell = ["BI_cell_data"+sz+".npy" , "EI_BG_cell_data"+sz+".npy", "EI_cell_data"+sz+".npy",
		"EI_Cu_cell_data"+sz+".npy", "RI_cell_data"+sz+".npy"]
pq = ["BI_quat_data"+sz+".npy" , "EI_BG_quat_data"+sz+".npy", "EI_quat_data"+sz+".npy",
		"EI_Cu_quat_data"+sz+".npy", "RI_quat_data"+sz+".npy"]

for i in range(len(cell)):
	np_to_block(cell[i], pq[i], Block_Size=(4,512,512,int(sz)))

