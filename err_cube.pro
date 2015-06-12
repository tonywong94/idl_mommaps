FUNCTION ERR_CUBE,  im, $
                    pattern=pattern,$
                    mask=mask,$
                    planes=planes
;+
; NAME:
;   ERR_CUBE
;
; PURPOSE:
;   generate a cube containing estimated errors for each pixel.
;   if the sensitivity pattern is not known, sigma_cube.pro might be
;   a better choice.
;
; INPUTS:
;   im        data cube
;   [pattern] sensitivity pattern (doesn't require normalization)
;             (e.g. "sen" maps from miriad/mossen or .nsen.fits from casa_mossen.pro)
;             could be a 2D image or 3D cube    
;   [planes]  line-free velocity planes for error estimations
;             default: [0,1,-1,-2] (first & last two channels)
;   [mask]    optional mask cube; zero values excluded from noise level estimation       
;
; OUTPUTS:
;   sen       error cube
;
; HISTORY:
;
;   20120710  RX  introduced
;   20130227  RX  change the name from find_sens.pro to err_cube.pro
;   20130319  RX  add options for different methods of error estimations
;                 change the keyword GAIN to PATTERN...
;                 PATTERN could be either 3D or 2D
;   20130603  RX  correctly handle the cube with partial coverage in 2D 
;                 or 3D (e.g. MAGMA)
;             
;-

; SELECT REGION AND NORMALIZE NOISE FOR ERROR ESTIMATIONS
dim=size(im,/dimension)
if  keyword_set(pattern) then begin
  img=pattern
  img[where(img le 0.0,/null)]=!VALUES.F_NAN
  img_min=min(img,/nan)
  print, "Current minimum rms value: ", img_min
  img=img/img_min
  if (size(img))[0] eq 2 then begin
    img=CMREPLICATE(img,dim[2])
  endif
endif else begin
  img=im*0.0+1.0
endelse

; missing data will be Nan in assi
assi=im*0.0
assi2d=total(im,3,/nan)*0.0+1.0    
if not keyword_set(planes) then begin
  sch=where(total(total(im eq im,1),1) ne 0)
  planes=[sch[0],sch[1],sch[-2],sch[-1]]
endif
for i=0,n_elements(planes)-1 do begin
  assi[*,*,planes[i]]=assi2d
endfor
if keyword_set(mask) $
  then assi[where(im ne im or mask eq 0.0,/null)]=!values.f_nan $
  else assi[where(im ne im,/null)]=!values.f_nan

ncube=im/img
ifail=0

; USE SELECTED REGIONS
tag=where(assi eq 1.0)
if tag[0] eq -1 then tmp=[0.0] else tmp=ncube[tag]
rms=ROBUST_SIGMA(tmp,/zero)
print, "Robustly estimated value:  ", rms
if rms eq -1 then ifail=1

sen=img*rms
print, "New min and max rms value: ", min(sen,/nan), max(sen,/nan)

;writefits, 'testx.sens.fits', sen
;writefits, 'testx.sens_assi.fits', assi
;writefits, 'testx.pattern.fits', pattern

return,sen

END
