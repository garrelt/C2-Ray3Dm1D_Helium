import numpy as np
import pylab as pl
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import ImageGrid
import sys

path = '/home/klee/myscratch/HeliumHKTest1.5/results/'


#Test 1 (constant vs adaptive)

#LTE_choice = 'do_LTE_'
#LTE_choice = 'no_LTE_'
 
#mesh = "mesh_3_"
#mesh = "mesh_10_"
#mesh = "mesh_20_"
mesh = "mesh_30_"
#mesh = "mesh_50_"

timestep = "t1_"
#timestep = "t2_"
#timestep = "t4_"
#timestep = "t5_"
#timestep = "t10_"
#timestep = "t20_"
#timestep = "t50_"
#timestep = "t100_"
#timestep = "t200_"
#timestep = "t500_"
#timestep = "t1000_"
#timestep = "t2000_"
#timestep = "t5000_"
#timestep = "t10000_"
#timestep = "t20000_"

#f_function = "f1_"
#f_function = "f2_"
#f_function = "f3_"
#f_function = "f4_"
#f_function = "f5_"
#f_function = "f6_"
#f_function = "f7_"
#f_function = "f8_"
#f_function = "f9_"
f_function = "f10_"
#f_function = "f11_"
#f_function = "f12_"
#f_function = "f13_"
#f_function = "f14_"
#f_function = "f15_"
#f_function = "f16_"
#f_function = "f17_"
#f_function = "f18_"
#f_function = "f19_"
#f_function = "f20_"

#source_choice = 'bb_'
source_choice = 'pl_'

'''
ct_do_LTE_input_file1 = path+'constant_do_LTE_'+timestep+source_choice+mesh+'Temper3d_8.884.bin'
ct_do_LTE_input_file2 = path+'constant_do_LTE_'+timestep+source_choice+mesh+'Temper3d_8.771.bin'
ct_do_LTE_input_file3 = path+'constant_do_LTE_'+timestep+source_choice+mesh+'Temper3d_8.662.bin'
ct_do_LTE_input_file4 = path+'constant_do_LTE_'+timestep+source_choice+mesh+'Temper3d_8.556.bin'
ct_do_LTE_input_file5 = path+'constant_do_LTE_'+timestep+source_choice+mesh+'Temper3d_8.452.bin'
ct_do_LTE_input_file6 = path+'constant_do_LTE_'+timestep+source_choice+mesh+'Temper3d_8.351.bin'
ct_do_LTE_input_file7 = path+'constant_do_LTE_'+timestep+source_choice+mesh+'Temper3d_8.253.bin'
ct_do_LTE_input_file8 = path+'constant_do_LTE_'+timestep+source_choice+mesh+'Temper3d_8.157.bin'

ct_no_LTE_input_file1 = path+'constant_no_LTE_'+timestep+source_choice+mesh+'Temper3d_8.884.bin'
ct_no_LTE_input_file2 = path+'constant_no_LTE_'+timestep+source_choice+mesh+'Temper3d_8.771.bin'
ct_no_LTE_input_file3 = path+'constant_no_LTE_'+timestep+source_choice+mesh+'Temper3d_8.662.bin'
ct_no_LTE_input_file4 = path+'constant_no_LTE_'+timestep+source_choice+mesh+'Temper3d_8.556.bin'
ct_no_LTE_input_file5 = path+'constant_no_LTE_'+timestep+source_choice+mesh+'Temper3d_8.452.bin'
ct_no_LTE_input_file6 = path+'constant_no_LTE_'+timestep+source_choice+mesh+'Temper3d_8.351.bin'
ct_no_LTE_input_file7 = path+'constant_no_LTE_'+timestep+source_choice+mesh+'Temper3d_8.253.bin'
ct_no_LTE_input_file8 = path+'constant_no_LTE_'+timestep+source_choice+mesh+'Temper3d_8.157.bin'
'''
at_do_LTE_input_file1 = path+'adaptive_do_LTE_'+f_function+source_choice+mesh+'Temper3d_8.884.bin'
at_do_LTE_input_file2 = path+'adaptive_do_LTE_'+f_function+source_choice+mesh+'Temper3d_8.771.bin'
at_do_LTE_input_file3 = path+'adaptive_do_LTE_'+f_function+source_choice+mesh+'Temper3d_8.662.bin'
at_do_LTE_input_file4 = path+'adaptive_do_LTE_'+f_function+source_choice+mesh+'Temper3d_8.556.bin'
at_do_LTE_input_file5 = path+'adaptive_do_LTE_'+f_function+source_choice+mesh+'Temper3d_8.452.bin'
at_do_LTE_input_file6 = path+'adaptive_do_LTE_'+f_function+source_choice+mesh+'Temper3d_8.351.bin'
at_do_LTE_input_file7 = path+'adaptive_do_LTE_'+f_function+source_choice+mesh+'Temper3d_8.253.bin'
at_do_LTE_input_file8 = path+'adaptive_do_LTE_'+f_function+source_choice+mesh+'Temper3d_8.157.bin'

at_no_LTE_input_file1 = path+'adaptive_no_LTE_'+f_function+source_choice+mesh+'Temper3d_8.884.bin'
at_no_LTE_input_file2 = path+'adaptive_no_LTE_'+f_function+source_choice+mesh+'Temper3d_8.771.bin'
at_no_LTE_input_file3 = path+'adaptive_no_LTE_'+f_function+source_choice+mesh+'Temper3d_8.662.bin'
at_no_LTE_input_file4 = path+'adaptive_no_LTE_'+f_function+source_choice+mesh+'Temper3d_8.556.bin'
at_no_LTE_input_file5 = path+'adaptive_no_LTE_'+f_function+source_choice+mesh+'Temper3d_8.452.bin'
at_no_LTE_input_file6 = path+'adaptive_no_LTE_'+f_function+source_choice+mesh+'Temper3d_8.351.bin'
at_no_LTE_input_file7 = path+'adaptive_no_LTE_'+f_function+source_choice+mesh+'Temper3d_8.253.bin'
at_no_LTE_input_file8 = path+'adaptive_no_LTE_'+f_function+source_choice+mesh+'Temper3d_8.157.bin'

mesh_file1 = np.loadtxt(path+'mesh_8.884.bin')
mesh_file2 = np.loadtxt(path+'mesh_8.771.bin')
mesh_file3 = np.loadtxt(path+'mesh_8.662.bin')
mesh_file4 = np.loadtxt(path+'mesh_8.556.bin')
mesh_file5 = np.loadtxt(path+'mesh_8.452.bin')
mesh_file6 = np.loadtxt(path+'mesh_8.351.bin')
mesh_file7 = np.loadtxt(path+'mesh_8.253.bin')
mesh_file8 = np.loadtxt(path+'mesh_8.157.bin')
'''
cd1_data = open(ct_do_LTE_input_file1, 'rb')
cd2_data = open(ct_do_LTE_input_file2, 'rb')
cd3_data = open(ct_do_LTE_input_file3, 'rb')
cd4_data = open(ct_do_LTE_input_file4, 'rb')
cd5_data = open(ct_do_LTE_input_file5, 'rb')
cd6_data = open(ct_do_LTE_input_file6, 'rb')
cd7_data = open(ct_do_LTE_input_file7, 'rb')
cd8_data = open(ct_do_LTE_input_file8, 'rb')

cn1_data = open(ct_no_LTE_input_file1, 'rb')
cn2_data = open(ct_no_LTE_input_file2, 'rb')
cn3_data = open(ct_no_LTE_input_file3, 'rb')
cn4_data = open(ct_no_LTE_input_file4, 'rb')
cn5_data = open(ct_no_LTE_input_file5, 'rb')
cn6_data = open(ct_no_LTE_input_file6, 'rb')
cn7_data = open(ct_no_LTE_input_file7, 'rb')
cn8_data = open(ct_no_LTE_input_file8, 'rb')
'''
ad1_data = open(at_do_LTE_input_file1, 'rb')
ad2_data = open(at_do_LTE_input_file2, 'rb')
ad3_data = open(at_do_LTE_input_file3, 'rb')
ad4_data = open(at_do_LTE_input_file4, 'rb')
ad5_data = open(at_do_LTE_input_file5, 'rb')
ad6_data = open(at_do_LTE_input_file6, 'rb')
ad7_data = open(at_do_LTE_input_file7, 'rb')
ad8_data = open(at_do_LTE_input_file8, 'rb')

an1_data = open(at_no_LTE_input_file1, 'rb')
an2_data = open(at_no_LTE_input_file2, 'rb')
an3_data = open(at_no_LTE_input_file3, 'rb')
an4_data = open(at_no_LTE_input_file4, 'rb')
an5_data = open(at_no_LTE_input_file5, 'rb')
an6_data = open(at_no_LTE_input_file6, 'rb')
an7_data = open(at_no_LTE_input_file7, 'rb')
an8_data = open(at_no_LTE_input_file8, 'rb')
'''
mesh = np.fromfile(cd1_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cd1_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cd1 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cd1_data.close()

mesh = np.fromfile(cd2_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cd2_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cd2 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cd2_data.close()

mesh = np.fromfile(cd3_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cd3_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cd3 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cd3_data.close()

mesh = np.fromfile(cd4_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cd4_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cd4 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cd4_data.close()

mesh = np.fromfile(cd5_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cd5_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cd5 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cd5_data.close()

mesh = np.fromfile(cd6_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cd6_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cd6 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cd6_data.close()

mesh = np.fromfile(cd7_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cd7_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cd7 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cd7_data.close()

mesh = np.fromfile(cd8_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cd8_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cd8 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cd8_data.close()

mesh = np.fromfile(cn1_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cn1_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cn1 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cn1_data.close()

mesh = np.fromfile(cn2_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cn2_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cn2 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cn2_data.close()

mesh = np.fromfile(cn3_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cn3_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cn3 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cn3_data.close()

mesh = np.fromfile(cn4_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cn4_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cn4 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cn4_data.close()

mesh = np.fromfile(cn5_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cn5_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cn5 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cn5_data.close()

mesh = np.fromfile(cn6_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cn6_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cn6 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cn6_data.close()

mesh = np.fromfile(cn7_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cn7_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cn7 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cn7_data.close()

mesh = np.fromfile(cn8_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cn8_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cn8 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cn8_data.close()
'''
mesh = np.fromfile(ad1_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(ad1_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
ad1 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
ad1_data.close()

mesh = np.fromfile(ad2_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(ad2_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
ad2 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
ad2_data.close()

mesh = np.fromfile(ad3_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(ad3_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
ad3 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
ad3_data.close()

mesh = np.fromfile(ad4_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(ad4_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
ad4 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
ad4_data.close()

mesh = np.fromfile(ad5_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(ad5_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
ad5 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
ad5_data.close()

mesh = np.fromfile(ad6_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(ad6_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
ad6 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
ad6_data.close()

mesh = np.fromfile(ad7_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(ad7_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
ad7 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
ad7_data.close()

mesh = np.fromfile(ad8_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(ad8_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
ad8 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
ad8_data.close()

mesh = np.fromfile(an1_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(an1_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
an1 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
an1_data.close()

mesh = np.fromfile(an2_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(an2_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
an2 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
an2_data.close()

mesh = np.fromfile(an3_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(an3_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
an3 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
an3_data.close()

mesh = np.fromfile(an4_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(an4_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
an4 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
an4_data.close()

mesh = np.fromfile(an5_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(an5_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
an5 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
an5_data.close()

mesh = np.fromfile(an6_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(an6_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
an6 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
an6_data.close()

mesh = np.fromfile(an7_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(an7_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
an7 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
an7_data.close()

mesh = np.fromfile(an8_data, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(an8_data, count=mesh_x*mesh_y*mesh_z, dtype='float64')
an8 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
an8_data.close()

pl.figure()
pl.xlim([mesh_file8[0], mesh_file1[len(mesh_file1)-1]])
pl.plot(mesh_file8[:],ad1[14,14,:],color='b',linestyle='dashed')
pl.plot(mesh_file1[:],an1[14,14,:],color='b')
pl.plot(mesh_file8[:],ad4[14,14,:],color='r',linestyle='dashed')
pl.plot(mesh_file1[:],an4[14,14,:],color='r')
pl.plot(mesh_file8[:],ad8[14,14,:],color='g',linestyle='dashed')
pl.plot(mesh_file1[:],an8[14,14,:],color='g')
pl.xlabel('Distance(cm)',fontsize=15)
pl.ylabel('Temperature(K)',fontsize=15)
pl.title('')
#pl.savefig('results/Test1_bb.pdf')
#pl.savefig('results/Test1_pl.pdf')
pl.show()

#Test 1 (f value and percentage error)



'''
ct_input_file8 = path+'constant_no_LTE_t20000_bb_Temper3d_8.157.bin'
cf8 = open(ct_input_file8, 'rb')
mesh = np.fromfile(cf8, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cf8, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cbb = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cf8.close()

at_input_bb_f1 = path+'adaptive_no_LTE_f1_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f1, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf1 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f2 = path+'adaptive_no_LTE_f2_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f2, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf2 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f3 = path+'adaptive_no_LTE_f3_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f3, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf3 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f4 = path+'adaptive_no_LTE_f4_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f4, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf4 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f5 = path+'adaptive_no_LTE_f5_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f5, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf5 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f6 = path+'adaptive_no_LTE_f6_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f6, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf6 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f7 = path+'adaptive_no_LTE_f7_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f7, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf7 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f8 = path+'adaptive_no_LTE_f8_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f8, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf8 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f9 = path+'adaptive_no_LTE_f9_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f9, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf9 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f10 = path+'adaptive_no_LTE_f10_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f10, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf10 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f11 = path+'adaptive_no_LTE_f11_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f11, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf11 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f12 = path+'adaptive_no_LTE_f12_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f12, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf12 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f13 = path+'adaptive_no_LTE_f13_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f13, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf13 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f14 = path+'adaptive_no_LTE_f14_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f14, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf14 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f15 = path+'adaptive_no_LTE_f15_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f15, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf15 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f16 = path+'adaptive_no_LTE_f16_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f16, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf16 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f17 = path+'adaptive_no_LTE_f17_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f17, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf17 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f18 = path+'adaptive_no_LTE_f18_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f18, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf18 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f19 = path+'adaptive_no_LTE_f19_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f19, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf19 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_bb_f20 = path+'adaptive_no_LTE_f20_bb_Temper3d_8.157.bin'
temp = open(at_input_bb_f20, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
abbf20 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

ct_input_file8 = path+'constant_no_LTE_t20000_pl_Temper3d_8.157.bin'
cf8 = open(ct_input_file8, 'rb')
mesh = np.fromfile(cf8, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cf8, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cpl = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cf8.close()

at_input_pl_f1 = path+'adaptive_no_LTE_f1_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f1, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf1 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f2 = path+'adaptive_no_LTE_f2_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f2, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf2 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f3 = path+'adaptive_no_LTE_f3_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f3, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf3 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f4 = path+'adaptive_no_LTE_f4_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f4, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf4 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f5 = path+'adaptive_no_LTE_f5_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f5, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf5 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f6 = path+'adaptive_no_LTE_f6_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f6, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf6 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f7 = path+'adaptive_no_LTE_f7_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f7, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf7 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f8 = path+'adaptive_no_LTE_f8_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f8, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf8 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f9 = path+'adaptive_no_LTE_f9_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f9, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf9 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f10 = path+'adaptive_no_LTE_f10_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f10, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf10 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f11 = path+'adaptive_no_LTE_f11_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f11, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf11 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f12 = path+'adaptive_no_LTE_f12_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f12, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf12 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f13 = path+'adaptive_no_LTE_f13_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f13, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf13 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f14 = path+'adaptive_no_LTE_f14_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f14, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf14 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f15 = path+'adaptive_no_LTE_f15_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f15, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf15 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f16 = path+'adaptive_no_LTE_f16_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f16, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf16 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f17 = path+'adaptive_no_LTE_f17_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f17, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf17 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f18 = path+'adaptive_no_LTE_f18_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f18, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf18 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f19 = path+'adaptive_no_LTE_f19_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f19, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf19 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

at_input_pl_f20 = path+'adaptive_no_LTE_f20_pl_Temper3d_8.157.bin'
temp = open(at_input_pl_f20, 'rb')
mesh = np.fromfile(temp, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(temp, count=mesh_x*mesh_y*mesh_z, dtype='float64')
aplf20 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
temp.close()

integer = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
bb_percentage = np.zeros(20)
pl_percentage = np.zeros(20)
cbb_sum = np.sum(cbb)
cpl_sum = np.sum(cpl)

hihi = abs(cbb[:,:,:]-abbf1[:,:,:])
print 'bb f1', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[0] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf2[:,:,:])
print 'bb f2', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[1] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf3[:,:,:])
print 'bb f3', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[2] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf4[:,:,:])
print 'bb f4', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[3] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf5[:,:,:])
print 'bb f5', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[4] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf6[:,:,:])
print 'bb f6', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[5] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf7[:,:,:])
print 'bb f7', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[6] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf8[:,:,:])
print 'bb f8', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[7] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf9[:,:,:])
print 'bb f9', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[8] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf10[:,:,:])
print 'bb f10', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[9] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf11[:,:,:])
print 'bb f11', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[10] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf12[:,:,:])
print 'bb f12', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[11] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf13[:,:,:])
print 'bb f13', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[12] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf14[:,:,:])
print 'bb f14', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[13] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf15[:,:,:])
print 'bb f15', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[14] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf16[:,:,:])
print 'bb f16', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[15] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf17[:,:,:])
print 'bb f17', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[16] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf18[:,:,:])
print 'bb f18', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[17] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf19[:,:,:])
print 'bb f19', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[18] = np.sum(hihi)/cbb_sum*100
hihi = abs(cbb[:,:,:]-abbf20[:,:,:])
print 'bb f20', np.sum(hihi)/cbb_sum*100, '%'
bb_percentage[19] = np.sum(hihi)/cbb_sum*100

hihi = abs(cpl[:,:,:]-aplf1[:,:,:])
print 'pl f1', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[0] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf2[:,:,:])
print 'pl f2', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[1] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf3[:,:,:])
print 'pl f3', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[2] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf4[:,:,:])
print 'pl f4', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[3] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf5[:,:,:])
print 'pl f5', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[4] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf6[:,:,:])
print 'pl f6', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[5] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf7[:,:,:])
print 'pl f7', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[6] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf8[:,:,:])
print 'pl f8', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[7] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf9[:,:,:])
print 'pl f9', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[8] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf10[:,:,:])
print 'pl f10', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[9] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf11[:,:,:])
print 'pl f11', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[10] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf12[:,:,:])
print 'pl f12', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[11] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf13[:,:,:])
print 'pl f13', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[12] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf14[:,:,:])
print 'pl f14', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[13] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf15[:,:,:])
print 'pl f15', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[14] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf16[:,:,:])
print 'pl f16', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[15] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf17[:,:,:])
print 'pl f17', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[16] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf18[:,:,:])
print 'pl f18', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[17] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf19[:,:,:])
print 'pl f19', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[18] = np.sum(hihi)/cpl_sum*100
hihi = abs(cpl[:,:,:]-aplf20[:,:,:])
print 'pl f20', np.sum(hihi)/cpl_sum*100, '%'
pl_percentage[19] = np.sum(hihi)/cpl_sum*100


pl.figure()
pl.plot(integer,bb_percentage,color='b',label='black body')
pl.plot(integer,pl_percentage,color='r',label='power-law')
pl.xlabel('f parameter')
pl.ylabel('absolute percentage difference')
pl.xlim(0,20)
pl.legend(loc=1)
pl.show()
'''

#Test 1 (do_LTE vs no_LTE)
'''
#timestep = "t1_"
#timestep = "t2_"
#timestep = "t5_"
#timestep = "t10_"
#timestep = "t20_"
#timestep = "t50_"
timestep = "t100_"
#timestep = "t200_"
#timestep = "t500_"
#timestep = "t1000_"
#timestep = "t2000_"
#timestep = "t5000_"
#timestep = "t10000_"
#timestep = "t20000_"

#f_function = "f1_"
#f_function = "f2_"
#f_function = "f3_"
#f_function = "f4_"
f_function = "f5_"
#f_function = "f6_"
#f_function = "f7_"
#f_function = "f8_"
#f_function = "f9_"
#f_function = "f10_"
#f_function = "f11_"
#f_function = "f12_"
#f_function = "f13_"
#f_function = "f14_"
#f_function = "f15_"
#f_function = "f16_"
#f_function = "f17_"
#f_function = "f18_"
#f_function = "f19_"
#f_function = "f20_"

source_choice = 'bb_'
#source_choice = 'pl_'

ct_do_LTE_input_file8 = path+'constant_do_LTE_'+timestep+source_choice+'Temper3d_8.157.bin'
ct_no_LTE_input_file8 = path+'constant_no_LTE_'+timestep+source_choice+'Temper3d_8.157.bin'
at_temper_do_LTE_input_file8 = path+'adaptive_do_LTE_'+f_function+source_choice+'Temper3d_8.157.bin'
#at_xHeII_do_LTE_input_file8 = path+'adaptive_do_LTE_'+f_function+source_choice+'xfrac3dHe2_8.157.bin'
at_temper_no_LTE_input_file8 = path+'adaptive_no_LTE_'+f_function+source_choice+'Temper3d_8.157.bin'
#at_xHeII_no_LTE_input_file8 = path+'adaptive_no_LTE_'+f_function+source_choice+'xfrac3dHe2_8.157.bin'

mesh_file8 = np.loadtxt(path+'mesh_8.157.bin')

cf_do_LTE_8 = open(ct_do_LTE_input_file8, 'rb')
cf_no_LTE_8 = open(ct_no_LTE_input_file8, 'rb')
af_temper_do_LTE_8 = open(at_temper_do_LTE_input_file8, 'rb')
#af_xHeII_do_LTE_8 = open(at_xHeII_do_LTE_input_file8, 'rb')
af_temper_no_LTE_8 = open(at_temper_no_LTE_input_file8, 'rb')
#af_xHeII_no_LTE_8 = open(at_xHeII_no_LTE_input_file8, 'rb')

mesh = np.fromfile(cf_do_LTE_8, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cf_do_LTE_8, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cd_do_LTE_8 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cf_do_LTE_8.close()

mesh = np.fromfile(cf_no_LTE_8, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(cf_no_LTE_8, count=mesh_x*mesh_y*mesh_z, dtype='float64')
cd_no_LTE_8 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
cf_no_LTE_8.close()

mesh = np.fromfile(af_temper_do_LTE_8, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(af_temper_do_LTE_8, count=mesh_x*mesh_y*mesh_z, dtype='float64')
ad_temper_do_LTE_8 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
af_temper_do_LTE_8.close()

#mesh = np.fromfile(af_xHeII_do_LTE_8, count=6, dtype='int32')
#mesh_x, mesh_y, mesh_z = mesh[1:4]
#data = np.fromfile(af_xHeII_do_LTE_8, count=mesh_x*mesh_y*mesh_z, dtype='float64')
#ad_xHeII_do_LTE_8 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
#af_xHeII_do_LTE_8.close()

mesh = np.fromfile(af_temper_no_LTE_8, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(af_temper_no_LTE_8, count=mesh_x*mesh_y*mesh_z, dtype='float64')
ad_temper_no_LTE_8 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
af_temper_no_LTE_8.close()

#mesh = np.fromfile(af_xHeII_no_LTE_8, count=6, dtype='int32')
#mesh_x, mesh_y, mesh_z = mesh[1:4]
#data = np.fromfile(af_xHeII_no_LTE_8, count=mesh_x*mesh_y*mesh_z, dtype='float64')
#ad_xHeII_no_LTE_8 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
#af_xHeII_no_LTE_8.close()

pl.figure()
pl.xlim([mesh_file8[0], mesh_file8[len(mesh_file8)-1]])
#pl.plot(mesh_file8[:],cd_do_LTE_8[0,0,:],color='b')
#pl.plot(mesh_file8[:],cd_no_LTE_8[0,0,:],color='g')
pl.plot(mesh_file8[:],ad_temper_do_LTE_8[14,14,:],color='b')
pl.plot(mesh_file8[:],ad_temper_no_LTE_8[14,14,:],color='g')
#pl.plot(mesh_file8[:],ad_xHeII_no_LTE_8[14,14,:],color='g')
#for i in range(len(mesh_file8)):
#  if ad_xHeII_no_LTE_8[0,0,i] >0.99:
#    pl.axvline(x=mesh_file8[i],color='r',ls='dashed')
pl.xlabel('Distance(cm)',fontsize=15)
pl.ylabel('Temperature(K)',fontsize=15)
pl.title('')
'''
pl.show()

