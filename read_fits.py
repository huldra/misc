from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

def plot_distribution(arr):
    """
    This function plots the value contained in the given
    array against the associated (1D)index
    """
    n = np.linspace(0,arr.size, arr.size)
    plt.plot(n,arr.flatten())
    plt.show()

# ------------------BIAS---------------------- #
nb = 11 # Number of bias files
# Reading weird numbered file separately
bias = fits.getdata('MOzi260054.fits',1)
for i in range(nb-1):
    bias += fits.getdata('MOzi26002%d.fits'%i,1)
bias = bias/nb
# -------------------------------------------- #

# -----------------FLATFIELDS----------------- #

nf = 3 # number of flatfields for each filter
flats_halpha = np.zeros(bias.shape)
flats_ufilter = np.zeros(bias.shape)
flats_rfilter = np.zeros(bias.shape)
flats_gfilter = np.zeros(bias.shape)

for i in range(nf):
    flats_halpha += fits.getdata('MOzi26003%d.fits'%(i+2),1)
    flats_ufilter += fits.getdata('MOzi26003%d.fits'%(i+5),1)
    flats_rfilter += fits.getdata('MOzi26004%d.fits'%(i+1),1)
    flats_gfilter += fits.getdata('MOzi25003%d.fits'%(i+1),1)

flats_halpha  = (flats_halpha/nf) - bias
flats_ufilter = (flats_ufilter/nf) - bias
flats_rfilter = (flats_rfilter/nf) - bias
flats_gfilter = (flats_gfilter/nf) - bias
# --------------------------------------------- #

# Avoiding zero division by setting the zero elements to 1
flats_halpha += abs(flats_halpha)<1e-15
flats_ufilter += abs(flats_ufilter)<1e-15
flats_rfilter += abs(flats_rfilter)<1e-15
flats_gfilter += abs(flats_gfilter)<1e-15


# --------------Science Pictures--------------- #
h_alpha = (fits.getdata('MOzi260048.fits',1) - bias)/flats_halpha
ufilter = (fits.getdata('MOzi260050.fits',1) - bias)/flats_ufilter
rfilter = (fits.getdata('MOzi260051.fits',1) - bias)/flats_rfilter
gfilter = (fits.getdata('MOzi260049.fits',1) - bias)/flats_gfilter
# --------------------------------------------- #

# Creating RGB-array
m = rfilter.shape[0]
n = rfilter.shape[1]
image = np.zeros((m,n,3))
for i in range(m):
    for j in range(n):
        image[i][j][:] = (rfilter[i][j],gfilter[i][j],ufilter[i][j])

# Normalizing
mean = np.zeros(3)
for i in range(3):
    mean[i] = np.sum(image[:][:][i])/(m*n)
for i in range(3):
    image[:][:][i] /= mean[i]

# Plotting each color separately
plt.imshow(image[:,:,0], 
        cmap='gray', 
        vmin=0.0, 
        vmax=0.08, 
        origin='lower')
plt.show()

plt.imshow(image[:,:,1], 
        cmap='gray', 
        vmin=0.0, 
        vmax=0.07, 
        origin='lower')
plt.show()

plt.imshow(image[:,:,2], 
        cmap='gray', 
        vmin=0.0, 
        vmax=0.02, 
        origin='lower')
plt.show()

# Showing RGB-image (Not sure how to adjust the colors (argh!))
plt.imshow(image[400:-50,400:-100,:],origin='lower') 
plt.show()
