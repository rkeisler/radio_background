import numpy as np
import ipdb
import matplotlib.pylab as pl
pl.ion()
datadir = 'data/'
from astropy.io import fits
import healpy as hp
import cPickle as pickle




def play(sample='CMASS'):
    lss = get_sdss_map(sample=sample)
    #lss = get_sdss_random_map(sample=sample)
    mask_lss = get_sdss_mask(sample=sample, quick=True)
    radio, mask_radio = get_haslam()
    #mask_gal = get_mask_gal(percentage_keep=40)
    mask = mask_lss*mask_radio

    # mask out the high-RA region
    ind=np.arange(len(mask))
    theta, phi = hp.pix2ang(512, ind)
    #phi_deg = phi*180./np.pi
    #wh = np.where((phi_deg>180.)&(phi_deg<270.))[0]
    #mask[wh] = 0.
    
    # Remove mean or median across mask area.
    whmask = np.where(mask>0.1)[0]
    lss -= np.median(lss[whmask])
    radio -= np.median(radio[whmask])
    # Take cross spectrum.
    mask_factor = np.mean(mask**2.)
    lmax=1000    
    cl = hp.anafast(lss*mask, map2=radio*mask, lmax=lmax)/mask_factor
    # Plot results.
    nl = len(cl)
    l=np.arange(nl)
    pl.figure(1); pl.clf()
    expon = 3
    pl.plot(l, cl*l**expon)
    pl.plot([0,lmax],[0,0],'k--')
    # Print some summary stats.
    y = cl*l**expon
    deltal=100.
    snsq=[]
    fsky = np.mean(mask)
    for i in range(int(lmax/deltal)):
        imin=i*deltal+50
        imax=(i+1)*deltal+50
        ytmp=y[imin:imax]
        ltmp=int(np.mean(l[imin:imax]))
        this_sn = ytmp.mean()/ytmp.std()*np.sqrt(deltal*fsky)
        print ltmp, this_sn
        snsq.append(this_sn**2.)
    print np.sqrt(np.sum(snsq))
    
    ipdb.set_trace()
    


def ra_dec_to_hpix(ra, dec, nside=512):
    from astropy.coordinates import FK5
    from astropy import units as u
    phi = ra*np.pi/180.
    theta = np.pi/2.-dec*np.pi/180.
    ind_hpix = hp.ang2pix(nside, theta, phi)
    nmap = np.bincount(ind_hpix, minlength=hp.nside2npix(nside))
    return nmap

    
def download_sdss():
    from os import system
    from urllib import urlretrieve
    files=['galaxy_DR10v8_CMASS_North.fits','galaxy_DR10v8_CMASS_South.fits',
           'galaxy_DR10v8_LOWZ_North.fits','galaxy_DR10v8_LOWZ_South.fits',
           'random1_DR10v8_CMASS_North.fits', 'random1_DR10v8_CMASS_South.fits', 
           'random1_DR10v8_LOWZ_North.fits', 'random1_DR10v8_LOWZ_South.fits', 
           'random2_DR10v8_CMASS_North.fits', 'random2_DR10v8_CMASS_South.fits', 
           'random2_DR10v8_LOWZ_North.fits', 'random2_DR10v8_LOWZ_South.fits']
    basepath = 'http://data.sdss3.org/sas/dr10/boss/lss/'
    for file in files:
        url = basepath + file + '.gz'
        savename = datadir + file + '.gz'
        urlretrieve(url, savename)
        system('gunzip '+savename)


def download_haslam():
    from urllib import urlretrieve
    basepath = 'http://lambda.gsfc.nasa.gov/data/foregrounds/haslam/'
    files = ['lambda_haslam408_dsds.fits','lambda_haslam408_nofilt.fits']
    for file in files:
        url = basepath + file
        savename = datadir + file
        urlretrieve(url, savename)

        

def get_haslam(nside_out=512, coordinates='eq', quick=True):
    savename = datadir+'haslam_'+coordinates+'_%i.pkl'%nside_out
    if quick: return pickle.load(open(savename,'r'))
    radio = hp.read_map(datadir+'lambda_haslam408_dsds.fits')
    if coordinates=='gal':
        radio_out = hp.ud_grade(radio, nside_out)
    if coordinates=='eq':
        from astropy.coordinates import FK5
        from astropy import units as u
        # Up-sample and convert from galactic to equatorial.
        nside_up = nside_out*2
        radio_up_gal = hp.ud_grade(radio, nside_up)        
        # Find the indices in an up-sampled *galactic* map that belong to these 
        # *equatorial* coordinates.
        theta, phi = hp.pix2ang(nside_up, np.arange(hp.nside2npix(nside_up)))
        ra = phi
        dec = np.pi/2.-theta
        coord = FK5(ra=ra, dec=dec, unit=(u.rad, u.rad))
        l_gal = coord.galactic.l.rad
        b_gal = coord.galactic.b.rad
        phi = l_gal
        theta = np.pi/2.-b_gal
        ind_up = hp.ang2pix(nside_up, theta, phi)
        radio_up_eq = radio_up_gal[ind_up]
        radio_out = hp.ud_grade(radio_up_eq, nside_out)
    mask_out = np.ones_like(radio_out)
    pickle.dump((radio_out, mask_out), open(savename,'w'))
    return radio_out, mask_out


def get_sdss_mask(sample='LOWZ', nside=512, fwhm1=5.0, fwhm2=5.0, quick=True):
    # FWHM1 and FWHM2 are in degrees.
    savename = datadir+'sdss_mask_'+sample+'_%i_%i.pkl'%(fwhm1,fwhm2)
    if quick: return pickle.load(open(savename,'r'))
    nmap = np.zeros(hp.nside2npix(nside))
    for ind in ['1','2']:
        for reg in ['North','South']:
            filename = datadir+'random'+ind+'_DR10v8_'+sample+'_'+reg+'.fits'
            data = fits.open(filename)[1].data
            nmap += ra_dec_to_hpix(data.RA, data.DEC, nside=nside)
    mask = nmap>0
    mask = hp.smoothing(mask, fwhm=fwhm1*np.pi/180., verbose=False)
    mask = mask>0.05
    mask = hp.sphtfunc.smoothing(mask, fwhm=fwhm2*np.pi/180., verbose=False)
    mask -= mask.min()
    mask /= mask.max()
    pickle.dump(mask, open(savename,'w'))
    return mask

    
def study_randoms():
    for sample in ['LOWZ','CMASS']:
        for reg in ['North','South']:
            filename = datadir+'galaxy_DR10v8_'+sample+'_'+reg+'.fits'
            filename1 = datadir+'random1'+'_DR10v8_'+sample+'_'+reg+'.fits'
            filename2 = datadir+'random2'+'_DR10v8_'+sample+'_'+reg+'.fits'
            ndata = len(fits.open(filename)[1].data)
            nr1 = len(fits.open(filename1)[1].data)
            nr2 = len(fits.open(filename2)[1].data)
            print sample, reg, 1.*nr1/ndata, 1.*nr2/ndata
            
def get_sdss_map(sample='CMASS',nside=512):
    nmap = np.zeros(hp.nside2npix(nside))
    for reg in ['North','South']:
        filename = datadir+'galaxy_DR10v8_'+sample+'_'+reg+'.fits'        
        data = fits.open(filename)[1].data
        if sample=='CMASS': zrange=[0.40,0.65]
        if sample=='LOWZ': zrange=[0.10,0.40]
        whok = np.where((data.Z>zrange[0])&(data.Z<zrange[1]))[0]
        nmap += ra_dec_to_hpix(data.RA[whok], data.DEC[whok], nside=nside)
    return nmap


def get_sdss_random_map(sample='CMASS',nside=512):
    nmap = np.zeros(hp.nside2npix(nside))
    for reg in ['North','South']:
        filename = datadir+'random1_DR10v8_'+sample+'_'+reg+'.fits'
        data = fits.open(filename)[1].data
        if sample=='CMASS': zrange=[0.40,0.65]
        if sample=='LOWZ': zrange=[0.10,0.40]
        whok = np.where((data.Z>zrange[0])&(data.Z<zrange[1]))[0]
        ngrab=int(len(whok)/50.)
        ra = np.random.choice(data.RA[whok], ngrab)
        dec = np.random.choice(data.DEC[whok], ngrab)     
        nmap += ra_dec_to_hpix(ra, dec, nside=nside)
    return nmap
    
                                        
def get_mask_gal(percentage_keep=40, nside_out=512, coordinates='eq', quick=True):
    import pyfits
    from astropy.coordinates import FK5
    from astropy import units as u    
    savename = datadir+'mask_GAL0%i_%i_'%(percentage_keep,nside_out)
    savename += coordinates+'_.pkl'
    if quick: return pickle.load(open(savename,'r'))

    pp = pyfits.open(datadir+'HFI_Mask_GalPlane_2048_R1.10.fits')
    mask_gal = pp[1].data['GAL0%i'%percentage_keep]
    mask_gal = hp.reorder(mask_gal, out='RING', inp='NESTED')
    if coordinates=='gal':
        mask_out = hp.ud_grade(mask_gal, nside_out)
    if coordinates=='eq':
        nside_up = nside_out*2
        mask_gal = hp.ud_grade(mask_gal, nside_up)
        # Find the indices in an up-sampled *galactic* map that belong to these 
        # *equatorial* coordinates.
        theta, phi = hp.pix2ang(nside_up, np.arange(hp.nside2npix(nside_up)))
        ra = phi
        dec = np.pi/2.-theta
        coord = FK5(ra=ra, dec=dec, unit=(u.rad, u.rad))
        l_gal = coord.galactic.l.rad
        b_gal = coord.galactic.b.rad
        phi = l_gal
        theta = np.pi/2.-b_gal
        ind_up = hp.ang2pix(nside_up, theta, phi)
        mask_up_eq = mask_gal[ind_up]
        mask_out = hp.ud_grade(mask_up_eq, nside_out)

    pickle.dump(mask_out, open(savename,'w'))
    return mask_out
