PRO HEXTRACT3D,oldim,oldhd,newim,newhd,$
    xyrange,zrange=zrange
;+
; NAME:
;   HEXTRACT3D
;   
; PURPOSE:
;   extract a subregion from a cube and keep the orginal astrometry
;   without regridding
;   (similar with writechunk.pro from E.Rosolowsky, without file I/O)
;
; INPUTS:
;   xyrange     --  [x0,x1,y0,y1] 
;                   bottom left / top right corner positions
;                   (x0,y0)-(x1,y1)
;   zrange      --  [z0,z1]
;
; HISTORY:
; 
;   20130912    RX  introduced
;
;-

newim=oldim[xyrange[0]:xyrange[1],xyrange[2]:xyrange[3],*]
tmphd=oldhd

SXDELPAR,tmphd,'NAXIS3'
SXDELPAR,tmphd,'NAXIS4'
SXADDPAR,tmphd,'NAXIS',2

HEXTRACT,oldim[*, *, 0],tmphd,tmpim,newhd,$
    xyrange[0],xyrange[1],xyrange[2],xyrange[3],/silent

SXADDPAR,newhd,'NAXIS',sxpar(oldhd, 'NAXIS')
tmp=SXPAR(oldhd,'NAXIS3',count=ct)
if  ct ne 0 $
    then SXADDPAR,newhd,'NAXIS3',sxpar(oldhd, 'NAXIS3'),after='NAXIS2'
tmp=SXPAR(oldhd,'NAXIS4',count=ct)
if  ct ne 0 $
    then SXADDPAR,newhd,'NAXIS4',sxpar(oldhd, 'NAXIS4'),after='NAXIS3'

if  n_elements(zrange) eq 2 then begin
    newim=newim[*,*,zrange[0]:zrange[1]]
    ;sxaddpar,newhd,'CRPIX3',sxpar(oldhd, 'CRPIX3')-zrange[0]
    sxaddpar,newhd,'CRVAL3',$
        sxpar(oldhd, 'CRVAL3')+zrange[0]*sxpar(oldhd, 'CDELT3')
endif

END



