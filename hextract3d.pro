PRO HEXTRACT3D,oldim,oldhd,newim,newhd,$
    xyrange,vrange=vrange
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
;   vrange      --  velocity range (not implemented)
;
; HISTORY:
; 
;   20130912    RX  introduced
;
;-

newhd=oldhd
newim=oldim[xyrange[0]:xyrange[1],xyrange[2]:xyrange[3],*]


; Update header information!
tmphd=oldhd
SXDELPAR,tmphd,'NAXIS3'
SXDELPAR,tmphd,'NAXIS4'
SXADDPAR,tmphd,'NAXIS',2
HEXTRACT,oldim[*, *, 0],tmphd,tmpim,newhd,$
    xyrange[0],xyrange[1],xyrange[2],xyrange[3],/silent

SXADDPAR,newhd,'NAXIS3',sxpar(oldhd, 'NAXIS3')
SXADDPAR,newhd,'NAXIS4',sxpar(oldhd, 'NAXIS4')

cpkey=['NAXIS3','NAXIS4']
for i=0,n_elements(cpkey)-1 do begin
    tmp=SXPAR(oldhd,cpkey[i],count=ct)
    if ct ne 0  then SXADDPAR,newhd,cpkey[i],sxpar(oldhd,cpkey[i]) $
                else SXDELPAR,newhd,cpkey[i]
endfor

SXADDPAR,newhd,'NAXIS',sxpar(oldhd, 'NAXIS')

END
