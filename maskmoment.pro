PRO MASKMOMENT, data, mask, ecube, vel, $
                mom0 = mom0, emom0 = emom0, $
                mom1 = mom1, emom1 = emom1, $
                mom2 = mom2, emom2 = emom2, $
                peak = peak, snrpk=snrpk, $
                mask0 = mask0, $
                nchref=nchref, emommx=emommx
;+
; NAME:
;   MASKMOMENT
;
; PURPOSE:
;   derive moment maps
;
; INPUTS:
;   DATA    -- data cube
;   MASK    -- mask cube
;   ECUBE   -- error cube 
;   VEL     -- velocity values
;   NCHREF  -- number of channels for calculating conservative emom0
;   
; KEYWORDS:  
;   MASK0   -- mask pixels where some channels are blank.
;              this option might be useful for fits files generated
;              by CASA when the primary beam masking is not constant across
;              different channels
;
; OUTPUTS:
;   MOMX, EMOMX, PEAK, SNRPK
;
; HISTORY:
;
;   20130214  TW  introduced
;   20130214  RX  fix a bug in emom0, add optional masking for pixels where 
;                 some channels are blank.
;                 mask negative values for mom1 & mom2
;   20130410  RX  enhance exception catching in mom1 & mom2
;                 set a minial value for emom0 by determining the smallest line width 
;                 add output values emom0_min and emom0_max for deriving 
;                 a better emom0 map
;   20140821  TW  Make the calculation of emom0 outside the mask optional,
;                 for consistency with other emom's and for varying spectral coverage
;   20150503  TW  Cleanup
;   20160201  TW  Change the way the secondary emom0 estimate is calculated
;-

sz = size(data)

;if n_elements(nchmin) eq 0 then nchmin=2.0

; FOV MASK
fov=float(total((data eq data and ecube eq ecube),3))
if  keyword_set(mask0) then begin
    fov[where(fov ne sz[3], /null)]=!values.f_nan
endif else begin
    fov[where(fov eq 0.0, /null)]=!values.f_nan
endelse

; MOMENT 0 AND UNCERTAINTY
mom0  = total(data*float(mask),3,/nan)+fov-fov
mom0[where(total(mask,3,/nan) eq 0.0, /null)]=!values.f_nan

emom0 = sqrt(total(ecube^2.0*float(mask),3,/nan))+fov-fov
emom0[where(total(mask,3,/nan) eq 0.0, /null)]=!values.f_nan

; CONSERVATIVE ERROR IN MOM0 USING ALL OR SPECIFIED NO. OF CHANNELS
if  nchref gt 0.0 then begin
    print, "Number of channels assumed in emom0max: ", nchref
    emom0_1ch=sqrt(total(ecube^2.0,3,/nan)/(total(ecube eq ecube,3)>1))+fov-fov
    ;nch_min=total(float(mask),3,/nan)>nchmin
    emommx = emom0_1ch*sqrt(nchref)
    ;tag=where(emom0 lt emom0_min)
    ;if tag[0] ne -1 then emom0[tag]=emom0_min[tag]
endif else begin
    ;vmaskcub = data*0.0
    ;vmaskarr = make_array(sz[3],/float,/nozero)
    ;for i=0,sz[3]-1 do vmaskarr[i]=max(max(float(mask),/nan,dim=1),/nan,dim=1)
    ;vmaskcub[*,*,where(vmaskarr GT 0,/null)] = 1.
    emommx = sqrt(total(ecube^2.0,3,/nan))+fov-fov
endelse

; MOMENT 1 & 2, PEAK, SNRPK
mom1  = fltarr(sz[1],sz[2])+!values.f_nan
mom2  = mom1
emom1 = mom1
emom2 = mom1
peak  = mom1
snrpk = mom1

for i = 0, sz[1]-1 do begin
  for j = 0, sz[2]-1 do begin
    spectrum = data[i, j, *]
    error = ecube[i, j, *]
    ; IGNORE MASK FOR SNR PEAK IMAGE
    ind2 = where((spectrum eq spectrum) and (error eq error), ct)
    if ct gt 0 then snrpk[i,j]=max(spectrum/error,/nan)+fov[i,j]-fov[i,j]
    ; CLIP NEGATIVE PIXELS FOR MOM-1 & MOM-2 
    ind2 = where(mask[i,j,*] and (data[i,j,*] gt 0.0), ct)
    if ct lt 2 then continue
    mom = WT_MOMENT(vel[ind2], spectrum[ind2], error = error[ind2])
    peak[i, j] = max(spectrum[ind2],/nan)+fov[i,j]-fov[i,j]
    ; RESTRICT MOM-1 TO ALLOWED VELOCITY RANGE
    if mom.mean ge min(vel[ind2],/nan) and mom.mean le max(vel[ind2]) then begin 
      mom1[i, j]  = mom.mean +fov[i,j]-fov[i,j]
      emom1[i, j] = mom.errmn+fov[i,j]-fov[i,j]
    endif
    if mom.stdev gt 0.0 then begin
      mom2[i, j]  = mom.stdev+fov[i,j]-fov[i,j]
      emom2[i, j] = mom.errsd+fov[i,j]-fov[i,j]
    endif
  endfor
endfor

END
