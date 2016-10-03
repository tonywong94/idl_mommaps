PRO mom_test
;+
;   mask the cube by first smoothing to 6" resolution and applying a
;   Gaussian smoothing of FWHM 10 km/s in velocity.  The 4-sigma contour
;   of the smoothed cube is expanded to 2-sigma to generate the mask.
;-
makemom,'../example/n4654co.line.cm.fits.gz',$
    errfile='../example/n4654co.line.sen.fits.gz',$
    baseroot='n4654co',$
    smopar=[6.,10.],$
    thresh=4,edge=2,/PVMOM0
            
END

PRO MOM_TEST_ROT
;+
;   align the galaxy major axis to x-axis
;   re-generate the mom01 and mom0pv with the rotated data
;-

baseroot='n4654co'
angle=30.0

data=readfits('../example/'+baseroot+'.line.cm.fits.gz',hd)
mask=readfits(baseroot+'.mask.fits',maskhd)
mom0=readfits(baseroot+'.mom0.fits',mom0hd)
emom0=readfits(baseroot+'.emom0.fits',emom0hd)

sz=size(data,/d)
hrot3d,mom0,mom0hd,newmom0,newmom0hd,angle,(sz[0]-1)/2.0,(sz[1]-1)/2.0,$
    1,missing=!values.f_nan,/pivot
hrot3d,emom0,emom0hd,newemom0,newemom0hd,angle,(sz[0]-1)/2.0,(sz[1]-1)/2.0,$
    1,missing=!values.f_nan,/pivot
hrot3d,data,hd,newdata,newhd,angle,(sz[0]-1)/2.0,(sz[1]-1)/2.0,$
    1,missing=!values.f_nan,/pivot
hrot3d,mask,hd,newmask,newmaskhd,angle,(sz[0]-1)/2.0,(sz[1]-1)/2.0,$
    1,missing=!values.f_nan,/pivot    

MASKMOMENT_PV,newdata,newhd,newmask,mom0xv,mom0vy,mom0xvhd=mom0xvhd,mom0vyhd=mom0vyhd,$
    vrange=vrange

WRITEFITS, baseroot+'_rot.cube.fits',newdata,newhd
WRITEFITS, baseroot+'_rot.mask.fits',newmask,newmaskhd        
WRITEFITS, baseroot+'_rot.mom0xv.fits',mom0xv,mom0xvhd
WRITEFITS, baseroot+'_rot.mom0vy.fits',mom0vy,mom0vyhd
WRITEFITS, baseroot+'_rot.mom0.fits',newmom0,newmom0hd
WRITEFITS, baseroot+'_rot.emom0.fits',newemom0,newemom0hd

PLTMOM_PV, baseroot+'_rot'
PLTMOM_PV, baseroot
PLTMOM,baseroot

END
