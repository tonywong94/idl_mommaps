PRO GKERNEL,  bmaj1, bmin1, bpa1, $
              bmaj2, bmin2, bpa2, $
              bmaj, bmin, bpa, $
              ifail,quiet=quiet
;+
; NAME:
;   GKERNEL
;
; PURPOSE:
;   determine the 2D kernel required for convolving a 2D Gaussian to another
;   2D Gaussian (just like MIRIAD-convol with options=final) 
;
; INPUTS:
;   bmaj1,bmin1        Major and minor FWHM of the target Gaussian
;   bpa1               Position angle of the target Gaussian, in degrees.
;   bmaj2,bmin2        Major and minor FWHM of the origianl Gaussian 
;   bpa2               Position angle of the original Gaussian, in degrees.
;
; OUTPUTS:
;   bmaj,bmin          Major and minor axes of the Kernel.
;   bpa                Position angle of the kernel radians.
;   ifail              Success status: 0  All OK.
;                                      1  Result is pretty close to a
;                                         point source.
;                                      2  Illegal result.
;
; HISTORY:
;   
;   20120514  RX       using the same formula in miriad/gaupar.for
;-

;
d2r=!dpi/180
r2d=180.0/!dpi
theta1 = bpa1*d2r
theta2 = bpa2*d2r

alpha  = (bmaj1*cos(theta1))^2 + (bmin1*sin(theta1))^2 -$
         (bmaj2*cos(theta2))^2 - (bmin2*sin(theta2))^2
beta   = (bmaj1*sin(theta1))^2 + (bmin1*cos(theta1))^2 -$
         (bmaj2*sin(theta2))^2 - (bmin2*cos(theta2))^2
gamma  = 2 * ( (bmin1^2-bmaj1^2)*sin(theta1)*cos(theta1) -$
               (bmin2^2-bmaj2^2)*sin(theta2)*cos(theta2) )

s = alpha + beta
t = sqrt((alpha-beta)^2 + gamma^2)
limit = min([bmaj1,bmin1,bmaj2,bmin2])
limit = 0.1*limit^2

if  ((alpha lt 0.) or (beta lt 0.) or (s lt t)) then begin
  bmaj = 0
  bmin = 0
  bpa = 0
  if  ((0.5*(s-t) lt limit) and (alpha gt -limit) and (beta gt -limit)) then $
    ifail = 1 else $
    ifail = 2
endif else begin
  bmaj = float(sqrt(0.5*(s+t)))
  bmin = float(sqrt(0.5*(s-t)))
  if  (abs(gamma)+abs(alpha-beta) eq 0.0) then $
    bpa = 0.0 else $
    bpa = float(0.5*atan(-gamma,alpha-beta)*r2d)
  ifail = 0
endelse
if  not keyword_set(quiet) then begin
  message,/info, "original beam: "+string(bmaj2)+','+string(bmin2)+',   ('+string(bpa2,format='(i4)')+')'
  message,/info, "target   beam: "+string(float(bmaj1))+','+string(float(bmin1))+',   ('+string(float(bpa1),format='(i4)')+')'
  message,/info, "convol kernel: "+string(bmaj)+','+string(bmin)+',   ('+string(bpa,format='(i4)')+')'
endif

END

PRO gkernel_test

gkernel,70.,10.,0.,69.,9.,45.,bmaj,bmin,bpa,ifail
print,bmaj,bmin,bpa,ifail

gkernel,70.,50.,0.,50.,9.,45.,bmaj,bmin,bpa,ifail
print,bmaj,bmin,bpa,ifail

END
