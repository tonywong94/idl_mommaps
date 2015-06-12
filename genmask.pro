FUNCTION PADDING,im,edge

; if edge>0 pad data with some zero-value pixels in each dimension 
; if edge<0 remove the padding pixels

if  edge gt 0 then begin
    struct=size(im, /STRUCTURE)
    imdims=size(im)
    tmp=MAKE_ARRAY(struct.dimensions[0:imdims[0]-1]+2*edge, TYPE=struct.TYPE)
    if imdims[0] eq 1 then tmp[edge:(edge+imdims[1]-1)]=im
    if imdims[0] eq 2 then tmp[edge:(edge+imdims[1]-1),edge:(edge+imdims[2]-1)]=im
    if imdims[0] eq 3 then tmp[edge:(edge+imdims[1]-1),edge:(edge+imdims[2]-1),edge:(edge+imdims[3]-1)]=im
    if imdims[0] eq 4 then tmp[edge:(edge+imdims[1]-1),edge:(edge+imdims[2]-1),edge:(edge+imdims[3]-1),edge:(edge+imdims[4]-1)]=im
endif
if  edge lt 0 then begin
    imdims=size(im)
    if imdims[0] eq 1 then tmp=im[-edge:(edge+imdims[1]-1)]
    if imdims[0] eq 2 then tmp=im[-edge:(edge+imdims[1]-1),-edge:(edge+imdims[2]-1)]
    if imdims[0] eq 3 then tmp=im[-edge:(edge+imdims[1]-1),-edge:(edge+imdims[2]-1),-edge:(edge+imdims[3]-1)]
    if imdims[0] eq 4 then tmp=im[-edge:(edge+imdims[1]-1),-edge:(edge+imdims[2]-1),-edge:(edge+imdims[3]-1),-edge:(edge+imdims[4]-1)]
endif
if  edge eq 0 then tmp=im

return,tmp
END


FUNCTION GENMASK, im, err=err, hd=hd,$
                  spar=spar,$
                  sig=sig,$
                  grow=grow,$
                  guard=guard
;+
; NAME:
;   GENMASK
;
; PURPOSE:
;   generate a signal mask from a data cube using sigma clipping, with optional
;   dilated masking and pre-smoothing.
;
; INPUTS:
;   IM        data cube
;   ERR       error cube
;   HD        data header
;   [SPAR]    2-element vector
;             We will degrade the cube angular resolution to <spar[0]> arcsec, and
;             smooth spectra with a Gaussian function of fwhm=<spar[1]> km/s for 
;             signal detection.
;   [SIG]     detection threshold, in units of cube/smoothed cube sigma
;   [GROW]    threshold for expanding core mask, in units of cube/smoothed cube sigma
;   [GUARD]   thickness of guard band in pixels
;
; OUTPUTS:
;   MASK      signal mask cube (1 as signal, 0 as noise, nan as missing data)
;
; HISTORY:
;
;   20130624  RX  introduced (for replacing grmask.pro and gsmask.pro)
;   20150528  TW  removed eprop option.
;   
;-

forward_function dilate_mask

if not keyword_set(spar) then spar=[0.0,0.0]
if not keyword_set(sig)  then sig =0.0
if not keyword_set(grow) then grow=0.0

im_smo=im
sen_smo=err

; SMOOTH DATA IF REQUESTED
if spar[0] gt 0.0 then begin
    nchan=(size(im,/d))[2]
    SMOOTH3D,im/err,hd,im_smo,hd_smo,[spar[0],spar[0],0.],svel=spar[1]
    sen_smo=ERR_CUBE(im_smo,planes=indgen(nchan))
    ; calculate rms after clipping based on MIRIAD sigest.for
    mask1 = abs(im_smo) lt 2.5*sqrt(!PI/2)*meanabsdev(im_smo,/nan)
    sen_smo=ERR_CUBE(im_smo,planes=indgen(nchan),mask=mask1)
endif

; GENERATE MASK BY SIGMA CLIPPING (MIN 2 CHANS)
if sig gt 0.0 then begin
    mask = im*0.0
    mask[where(im_smo gt sen_smo*sig,/null)]=1.0
    mask = padding(mask,1)
    mask = (mask*(shift(mask,0,0,1)+shift(mask,0,0,-1)) gt 0)
    ; DILATE THE MASK IF REQUESTED
    if grow gt 0.0 then begin
        ; CONSTRAINT MASK
        growmask = im_smo gt grow*sen_smo
        growmask = padding(growmask,1)
        growmask = (growmask*(shift(growmask,0,0,1)+shift(growmask,0,0,-1)) gt 0)
        ; EXPAND MASK
        mask = float(dilate_mask(mask, constraint = growmask))
    endif
    mask = padding(mask,-1)  
    ; ADD A GUARD BAND IF REQUESTED
    if total(guard) gt 0 then begin
        mask = MASKGUARD_3D(mask, guard = guard)
    endif
endif else begin
    mask = im*0.0 + 1.0
endelse

return,mask
  
END
