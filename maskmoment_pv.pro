PRO MASKMOMENT_PV, data, hd, mask, $
                   mom0xv, mom0vy,$
                   mom0xvhd=mom0xvhd, mom0vyhd=mom0vyhd,$
                   keep3d=keep3d, mask0=mask0,$
                   vrange=vrange
;+
; NAME:
;   MASKMOMENT
;
; PURPOSE:
;   derive moment maps by collapsing along x/y-axis of the cube 
;
; INPUTS:
;   DATA        -- data cube
;   HD          -- data header
;   MASK        -- mask cube (1: unmasked, 0: masked, !values.f_nan: no data)
;                  if no mask is available, just set it to 1.0
;   VRANGE      -- choose the velocity range for moments calculations
;
; KEYWORDS:
;   mask0       -- mask 2D pixels where some 3D pixels along LOS are blank.
;   keep3d      -- keep the dummy axis
;
; OUTPUTS:
;   MOM0XV,  MOM0VY
;   MOMXVHD, MOM0VYHD
;   
; NOTES:
;   By default, MOMXV and MOMVY will become 2D (|v-x & |y-v). 
;   If keep3d is set, the output data will still have the same number of dimensions 
;   as the orginal data. Then, the procedures handling ppv data (e.g. rd_hd.pro) will still work.
;   
;   The collapsing is done along the pixel directions rather than coodinate directions, so the 2D 
;   header might not be very meaningful (execept proj=CAR) 
;   
; HISTORY:
;
;   20131101    RX  introduced
;   20150602    TW  adapted for mommaps package
;-

; COLLAPSE CUBE IN Y
fov=float(total((data eq data),2))
if  keyword_set(mask0) then begin
    fov[where(fov ne max(fov), /null)]=!values.f_nan
endif else begin
    fov[where(fov eq 0.0, /null)]=!values.f_nan
endelse
mom0xv=data[*,0,*]
mom0xv[*,0,*]=total(data*float(mask),2,/nan)+fov-fov

; COLLAPSE CUBE IN X
fov=float(total((data eq data),1))
if  keyword_set(mask0) then begin
    fov[where(fov ne max(fov), /null)]=!values.f_nan
endif else begin
    fov[where(fov eq 0.0, /null)]=!values.f_nan
endelse
mom0vy=data[0,*,*]
mom0vy[0,*,*]=total(data*float(mask),1,/nan)+fov-fov

; SCALE VALUES
;getrot,hd,rotang,cdelt
RADIOHEAD, hd, s = h
xsize=abs(h.cdelt[0]*60.*60.)
ysize=abs(h.cdelt[1]*60.*60.)
if  n_elements(vrange) ne 0 then begin
    tag_invrange=where(h.v ge vrange[0] and h.v le vrange[1])
endif else begin
    vrange=[min(h.v),max(h.v)]
    tag_invrange=indgen(n_elements(h.v))
endelse

mom0vy=mom0vy*ysize
mom0xv=mom0xv*xsize

bunit=sxpar(hd,'BUNIT')
mhd=hd
sxaddpar,mhd,'BUNIT',strtrim(bunit,2)+'.arcsec'
mom0xvhd=mhd
mom0vyhd=mhd

sxaddpar, mom0xvhd, 'DATAMAX', max(mom0xv,/nan)
sxaddpar, mom0xvhd, 'DATAMIN', min(mom0xv,/nan)
sxaddpar, mom0vyhd, 'DATAMAX', max(mom0vy,/nan)
sxaddpar, mom0vyhd, 'DATAMIN', min(mom0vy,/nan)

; DROP THE DUMMY DIMENSION
if  not keyword_set(keep3d) then begin
    
    keylist=['CDELT','CRPIX','CRVAL','CTYPE','NAXIS','CUNIT','CROTA']

    foreach key,keylist do begin
        tmp=sxpar(mom0xvhd,key+'3',count=ct)
        if ct ne 0 then $
            SXADDPAR,mom0xvhd,key+'2',sxpar(mhd,key+'3') $
        else $
            SXDELPAR,mom0xvhd,key+'2'
        SXDELPAR,mom0xvhd,key+'3'
    endforeach
    mom0xv=total(mom0xv[*,*,[tag_invrange]],2) ; |V-RA
    mom0xv[where(total(mask[*,*,[tag_invrange]],2,/nan) eq 0.0, /null)]=!values.f_nan
    SXADDPAR,mom0xvhd,'CRPIX2',sxpar(mom0xvhd,'CRPIX2')-min(tag_invrange)
    SXADDPAR,mom0xvhd,'NAXIS2',n_elements(tag_invrange)
    
    foreach key,keylist do begin
        tmp=sxpar(mom0vyhd,key+'3',count=ct)
        if ct ne 0 then $
            SXADDPAR,mom0vyhd,key+'1',sxpar(mhd,key+'3') $
        else $
            SXDELPAR,mom0vyhd,key+'1'
        SXDELPAR,mom0vyhd,key+'3'
    endforeach
    mom0vy=total(mom0vy[*,*,[tag_invrange]],1) ; |V-DEC
    mom0vy[where(total(mask[*,*,[tag_invrange]],1,/nan) eq 0.0, /null)]=!values.f_nan
    mom0vy=rotate(mom0vy,4)                    ; |DEC-V
    SXADDPAR,mom0vyhd,'CRPIX1',sxpar(mom0vyhd,'CRPIX1')-min(tag_invrange)
    SXADDPAR,mom0vyhd,'NAXIS1',n_elements(tag_invrange)    
    
endif

END


