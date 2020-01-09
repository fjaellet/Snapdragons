# snapdragons (Barcelona flavour)

Snapdragons is a population synthesis &amp; mock survey data generation code - it creates synthetic stars from an N-body simulation of a galaxy. 
Please cite the Snapdragons paper ([Hunt et al. 2015, MNRAS, 450, 2132](http://adsabs.harvard.edu/abs/2015MNRAS.450.2132H) if you use it for publications.

This version is used in some of the hands-on session for the [Gaia-COST school in Barcelona ](https://github.com/Santiastro1/Gaia_School_BCN). It differs only in some superficial respects from the [original code](https://github.com/JASHunt/Snapdragons). Main contributors apart from [Jason Hunt](https://github.com/JASHunt) are [Santi Roca](https://github.com/Santiastro1) and [Mercé Romero](https://github.com/mromerog). 

## How does this work?

It's that simple: From the input: 
* N-body model (with positions, masses, ages, and metallicities) 
* isochrone tables
* (2D/3D) extinction model
* initial mass function
* limiting magnitude
* Gaia error model
snapdragons creates the output: 
* Synthetic stellar data

The code is in written in fortran. It has a Makefile in the main directory, which compiles the source code in src. To compile just type 'make' in the main directory.
The main module is called src/PopSynth.F, several options can be specified in Mockinioptions.h

You must dowload the extinction maps separately from https://1drv.ms/u/s!AnMgd4KLTt9j6HtlSWaRd9kP5Rem?e=QcCqax and place them in the ini folder (not inside UBV)

--------

## Changing input parameters

### Input simulation

The input file needs to be an ascii table with the following structure:
    * x, y, z(kpc), vx, vy, vz(km/s), metallicity(z), age(Gyr), mass(solar) [, other columns]
The name of the input file is hidden in PopSynth.F - so it might be easiest to replace the default input file (ini/input.dat) manually/using a python wrapper program.

### Input parameters

Mockinioptions.h lets you choose some basic things: 
* magnitude limit (in the V band): ChosenLim=14.0d0 
* whether to convolve the results with the Gaia error model: error=1  
* whether to output an ascii or a binary file: Binary=0

Don't forget to recompile before rerunning.

## Output

The default snapdragons output is stored in output/ as ascii or binary file without header.. There is, however, a python script in py/ that will convert this output to standard FITS file in the same directory:
    
    ``cd py``
    ``python snapdragons2fits.py``

Output columns are: 

* 1 true ra
* 2 true dec
* 3 true parallax
* 4 true pmra
* 5 true pmdec
* 6 true rv
* 7 observed ra
* 8 observed dec
* 9 observed parallax
* 10 observed pmra
* 11 observed pmdec
* 12 observed  rv
* 13 ra standard deviation 
* 14 dec standard deviation 
* 15 parallax standard deviation 
* 16 pmra standard deviation 
* 17 pmdec standard deviation 
* 18 rv standard deviation 
* 19 G Magnitude 
* 20 G_bp - G_rp
* 21 observed G Magnitude 
* 22 observed G_bp - G_rp 
* 23 G Magnitude Gaia standard deviation 
* 24 G_bp - G_rp Gaia standard deviation 
* 25 apparent V mag
* 26 distance from sun      
* 27 x 
* 28 y 
* 29 z 
* 30 vx 
* 31 vy 
* 32 vz 
* 33 mass (solar mass)
* 34 G_RVS
* 35 V-Igen
* 36 extinction
* 37 absolute V mag
* 38 absolute V-i 
* 39 metallicity z*10^6
* 40 Age (log10years)
* 41 observed distance
* 42 observed x		
* 43 observed y 
* 44 observed z 
* 45 observed vx 
* 46 observed vy 
* 47 observed vz 
* 48 observed galactocentric vr
* 49 observed vrot
* 50 galactocentric radius
* 51 galactocentric vr
* 52 vrot
* 53 observed galactocentric radius
* 54 pmra (km/s)
* 55 pmdec (km/s)
* 56 observed pmra (km/s)
* 57 observed pmdec (km/s)
* 58 l
* 59 b
* 60 observed l
* 61 observed b
* 62 vl 
* 63 vb
* 64 observed vl
* 65 observed vb
* 66 u mag
* 67 b mag
* 68 r mag
* 69 j mag
* 70 h mag
* 71 k mag
* 72 Error on U
* 73 Error on V
* 74 error on W
* 75 Teff
* 76 logg
* 77 observed Teff
* 78 observed logg 
* 79 Teff uncertainty
* 80 logg uncertainty









