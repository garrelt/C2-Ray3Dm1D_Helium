import numpy as np
import pylab as pl
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import ImageGrid
import sys

input_file = sys.argv[1]
if (input_file == 'ct'):
    input_file = 'constant_Temper3d'
elif (input_file == 'cx'):
    input_file = 'constant_xfrac3d'  
elif (input_file == 'c1'):
    input_file = 'constant_xfrac3dHe1'  
elif (input_file == 'c2'):
    input_file = 'constant_xfrac3dHe2'  
elif (input_file == 'at'):
    input_file = 'adaptive_Temper3d'
elif (input_file == 'ax'):
    input_file = 'adaptive_xfrac3d'  
elif (input_file == 'a1'):
    input_file = 'adaptive_xfrac3dHe1'  
elif (input_file == 'a2'):
    input_file = 'adaptive_xfrac3dHe2'  
elif (input_file == 't'):
    input_file_c = 'Temper3d'
elif (input_file == 'x'):
    input_file = 'xfrac3d'  
elif (input_file == '1'):
    input_file = 'xfrac3dHe1'  
elif (input_file == '2'):
    input_file = 'xfrac3dHe2'    
elif (input_file == 'i'):
    input_file = 'IonRates3'
elif (input_file == 's'):
    input_file = 'sca'
elif (input_file == 'l'):
    input_file = 'LTE'

path = '/home/klee/myscratch/HeliumHK/results/'

input_file1 = path+input_file+'_8.942.bin'
input_file2 = path+input_file+'_8.884.bin'
input_file3 = path+input_file+'_8.827.bin'
input_file4 = path+input_file+'_8.771.bin'
input_file5 = path+input_file+'_8.716.bin'
input_file6 = path+input_file+'_8.662.bin'
input_file7 = path+input_file+'_8.608.bin'
input_file8 = path+input_file+'_8.556.bin'

mesh_file1 = np.loadtxt(path+'mesh_8.942.bin')
mesh_file2 = np.loadtxt(path+'mesh_8.884.bin')
mesh_file3 = np.loadtxt(path+'mesh_8.827.bin')
mesh_file4 = np.loadtxt(path+'mesh_8.771.bin')
mesh_file5 = np.loadtxt(path+'mesh_8.716.bin')
mesh_file6 = np.loadtxt(path+'mesh_8.662.bin')
mesh_file7 = np.loadtxt(path+'mesh_8.608.bin')
mesh_file8 = np.loadtxt(path+'mesh_8.556.bin')

#Test 4
'''
input_file1 = path+input_file+'_0.125.bin'
input_file2 = path+input_file+'_0.250.bin'
input_file3 = path+input_file+'_0.375.bin'
input_file4 = path+input_file+'_0.500.bin'
input_file5 = path+input_file+'_0.625.bin'
input_file6 = path+input_file+'_0.750.bin'
input_file7 = path+input_file+'_0.875.bin'
input_file8 = path+input_file+'_1.000.bin'
'''

'''
mesh_file1 = np.loadtxt(path+'mesh_0.125.bin')
mesh_file2 = np.loadtxt(path+'mesh_0.250.bin')
mesh_file3 = np.loadtxt(path+'mesh_0.375.bin')
mesh_file4 = np.loadtxt(path+'mesh_0.500.bin')
mesh_file5 = np.loadtxt(path+'mesh_0.625.bin')
mesh_file6 = np.loadtxt(path+'mesh_0.750.bin')
mesh_file7 = np.loadtxt(path+'mesh_0.875.bin')
mesh_file8 = np.loadtxt(path+'mesh_1.000.bin')
'''
#Test 4
'''
#input_file9 = path+'ifront.dat'
#input_file8 = '/home/klee/myscratch/HydrogenHK/3D/results/short_temper.bin'
#input_file8 = '/home/klee/myscratch/HydrogenHK/3D/results/long_temper.bin'
#input_file8 = '/home/klee/myscratch/HydrogenHK/3D/results/pyramid_temper.bin'
#input_file8 = '/home/klee/myscratch/HydrogenHK/3D/results/pyramid_xfrac_adap200.bin'
#input_file8 = '/home/klee/myscratch/HydrogenHK/3D/results/pyramid_temper_adap200.bin'
#input_file8 = '/home/klee/myscratch/HydrogenHK/3D/results/short_xfrac_adap200.bin'
#input_file8 = '/home/klee/myscratch/HydrogenHK/3D/results/short_temper_adap200.bin'
'''

f1 = open(input_file1, 'rb')
f2 = open(input_file2, 'rb')
f3 = open(input_file3, 'rb')
f4 = open(input_file4, 'rb')
f5 = open(input_file5, 'rb')
f6 = open(input_file6, 'rb')
f7 = open(input_file7, 'rb')
f8 = open(input_file8, 'rb')

mesh = np.fromfile(f1, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(f1, count=mesh_x*mesh_y*mesh_z, dtype='float64')
d1 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
f1.close()

mesh = np.fromfile(f2, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(f2, count=mesh_x*mesh_y*mesh_z, dtype='float64')
d2 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
f2.close()

mesh = np.fromfile(f3, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(f3, count=mesh_x*mesh_y*mesh_z, dtype='float64')
d3 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
f3.close()

mesh = np.fromfile(f4, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(f4, count=mesh_x*mesh_y*mesh_z, dtype='float64')
d4 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
f4.close()

mesh = np.fromfile(f5, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(f5, count=mesh_x*mesh_y*mesh_z, dtype='float64')
d5 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
f5.close()

mesh = np.fromfile(f6, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(f6, count=mesh_x*mesh_y*mesh_z, dtype='float64')
d6 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
f6.close()

mesh = np.fromfile(f7, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(f7, count=mesh_x*mesh_y*mesh_z, dtype='float64')
d7 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
f7.close()

mesh = np.fromfile(f8, count=6, dtype='int32')
mesh_x, mesh_y, mesh_z = mesh[1:4]
data = np.fromfile(f8, count=mesh_x*mesh_y*mesh_z, dtype='float64')
d8 = data.reshape((mesh_x, mesh_y, mesh_z), order='Fortran')
f8.close()

'''
pl.figure()
pl.xlim([mesh_file1[0], mesh_file1[29]])
pl.plot(mesh_file1[:],d1[29,29,:])
pl.plot(mesh_file1[:],d2[29,29,:])
pl.plot(mesh_file1[:],d3[29,29,:])
pl.plot(mesh_file1[:],d4[29,29,:])
pl.plot(mesh_file1[:],d5[29,29,:])
pl.plot(mesh_file1[:],d6[29,29,:])
pl.plot(mesh_file1[:],d7[29,29,:])
pl.plot(mesh_file1[:],d8[29,29,:])
pl.show()
'''

fig = pl.figure(figsize=(16.0,5.0))
grid = ImageGrid(fig, 111, # similar to subplot(111)
                nrows_ncols = (3, 8),
                axes_pad=0.1, # pad between axes in inch.
                cbar_mode='single',
                )

im=grid[0].imshow(d1[0,:,:],norm=LogNorm(vmin=100,vmax=100000))
grid.cbar_axes[0].colorbar(im)
grid[1].imshow(d2[0,:,:],norm=LogNorm(vmin=100,vmax=100000))
grid[2].imshow(d3[0,:,:],norm=LogNorm(vmin=100,vmax=100000))
grid[3].imshow(d4[0,:,:],norm=LogNorm(vmin=100,vmax=100000))
grid[4].imshow(d5[0,:,:],norm=LogNorm(vmin=100,vmax=100000))
grid[5].imshow(d6[0,:,:],norm=LogNorm(vmin=100,vmax=100000))
grid[6].imshow(d7[0,:,:],norm=LogNorm(vmin=100,vmax=100000))
grid[7].imshow(d8[0,:,:],norm=LogNorm(vmin=100,vmax=100000))
grid[8].imshow(d1[:,0,:],norm=LogNorm(vmin=100,vmax=100000))
grid[9].imshow(d2[:,0,:],norm=LogNorm(vmin=100,vmax=100000))
grid[10].imshow(d3[:,0,:],norm=LogNorm(vmin=100,vmax=100000))
grid[11].imshow(d4[:,0,:],norm=LogNorm(vmin=100,vmax=100000))
grid[12].imshow(d5[:,0,:],norm=LogNorm(vmin=100,vmax=100000))
grid[13].imshow(d6[:,0,:],norm=LogNorm(vmin=100,vmax=100000))
grid[14].imshow(d7[:,0,:],norm=LogNorm(vmin=100,vmax=100000))
grid[15].imshow(d8[:,0,:],norm=LogNorm(vmin=100,vmax=100000))
grid[16].imshow(d1[:,:,0],norm=LogNorm(vmin=100,vmax=100000))
grid[17].imshow(d2[:,:,0],norm=LogNorm(vmin=100,vmax=100000))
grid[18].imshow(d3[:,:,0],norm=LogNorm(vmin=100,vmax=100000))
grid[19].imshow(d4[:,:,0],norm=LogNorm(vmin=100,vmax=100000))
grid[20].imshow(d5[:,:,0],norm=LogNorm(vmin=100,vmax=100000))
grid[21].imshow(d6[:,:,0],norm=LogNorm(vmin=100,vmax=100000))
grid[22].imshow(d7[:,:,0],norm=LogNorm(vmin=100,vmax=100000))
grid[23].imshow(d8[:,:,0],norm=LogNorm(vmin=100,vmax=100000))
pl.show()


'''
pl.figure(figsize=(16.0,5.0))
pl.subplot(3,8,1)
pl.imshow(d1[29,:,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,2)
pl.imshow(d2[29,:,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,3)
pl.imshow(d3[29,:,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,4)
pl.imshow(d4[14,:,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,5)
pl.imshow(d5[14,:,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,6)
pl.imshow(d6[14,:,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,7)
pl.imshow(d7[14,:,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,8)
pl.imshow(d8[14,:,:],norm=LogNorm(vmin=100,vmax=100000))

pl.colorbar()

pl.subplot(3,8,9)
pl.imshow(d1[:,14,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,10)
pl.imshow(d2[:,14,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,11)
pl.imshow(d3[:,14,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,12)
pl.imshow(d4[:,14,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,13)
pl.imshow(d5[:,14,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,14)
pl.imshow(d6[:,14,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,15)
pl.imshow(d7[:,14,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,16)
pl.imshow(d8[:,14,:],norm=LogNorm(vmin=100,vmax=100000))

pl.colorbar()

pl.subplot(3,8,17)
pl.imshow(d1[:,:,14],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,18)
pl.imshow(d2[:,:,14],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,19)
pl.imshow(d3[:,:,14],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,20)
pl.imshow(d4[:,:,14],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,21)
pl.imshow(d5[:,:,14],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,22)
pl.imshow(d6[:,:,14],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,23)
pl.imshow(d7[:,:,14],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,24)
pl.imshow(d8[:,:,14],norm=LogNorm(vmin=100,vmax=100000))

pl.colorbar()

pl.show()
'''

#pl.imshow(d8[:,:,14])
#pl.imshow(d8[:,14,:])
#pl.imshow(d8[14,:,:])



#pl.imshow(dh[:,:,14])
#pl.show()

#pl.semilogy(data1D[:,0],d8[:,14,14])

'''
pl.plot(data1D[:,0],d8[:,14,14])

pl.show()
'''


pl.figure()

pl.plot(mesh_file1[:],d1[0,0,:])
pl.plot(mesh_file1[:],d2[0,0,:])
pl.plot(mesh_file1[:],d3[0,0,:])
pl.plot(mesh_file1[:],d4[0,0,:])
pl.plot(mesh_file1[:],d5[0,0,:])
pl.plot(mesh_file1[:],d6[0,0,:])
pl.plot(mesh_file1[:],d7[0,0,:])
pl.plot(mesh_file1[:],d8[0,0,:])
pl.show()
'''
pl.plot(data1D[:,0],d1[14,14,:])
pl.plot(data1D[:,0],d2[14,14,:])
pl.plot(data1D[:,0],d3[14,14,:])
pl.plot(data1D[:,0],d4[14,14,:])
pl.plot(data1D[:,0],d5[14,14,:])
pl.plot(data1D[:,0],d6[14,14,:])
pl.plot(data1D[:,0],d7[14,14,:])
pl.plot(data1D[:,0],d8[14,14,:])


pl.plot(data1D[:,0],d1[46,46,:])
pl.plot(data1D[:,0],d2[46,46,:])
pl.plot(data1D[:,0],d3[46,46,:])
pl.plot(data1D[:,0],d4[46,46,:])
pl.plot(data1D[:,0],d5[46,46,:])
pl.plot(data1D[:,0],d6[46,46,:])
pl.plot(data1D[:,0],d7[46,46,:])
pl.plot(data1D[:,0],d8[46,46,:])
pl.title('hihi')
'''
#pl.plot(data1D[:,0],d7[14,14,:])
#pl.plot(data1D[:,0],d8[14,14,:])
#pl.plot(data1D[:,0],d1[14,:,14])
#pl.plot(data1D[:,0],d1[:,14,14])

#pl.show()
