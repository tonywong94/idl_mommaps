PRO mom_test
; mask the cube by first smoothing to 6" resolution and applying a
; Gaussian smoothing of FWHM 10 km/s in velocity.  The 4-sigma contour
; of the smoothed cube is expanded to 2-sigma to generate the mask.
makemom,'n4654co.line.cm.fits',$
    errfile='n4654co.line.sen.fits',$
    smopar=[6.,10.],$
    thresh=4,edge=2
END

