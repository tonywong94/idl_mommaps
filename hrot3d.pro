PRO hrot3d, oldim, oldhd, newim, newhd, $
    angle, xc, yc, $
    int, MISSING=missing, $
    INTERP = interp, CUBIC = cubic, PIVOT = pivot,ERRMSG= errmsg
;+
; NAME:
;   HROT3D
;   
; PURPOSE:
;   Rotate an image or a cube and create new FITS header with updated astrometry
;   A wrapper for astron/hrot.pro
;
; INPUTS:
;   see hrot.pro
;
; HISTORY:
; 
;   20150912    RX  introduced (similar with cube_hrot.pro in cpropstoo)
;
;-

sz=size(oldim)

if  sz[0] gt 2 then begin
    newim=oldim
    newhd=oldhd
    for i=0,sz[3]-1 do begin
        oldhd_i=oldhd
        sxaddpar,oldhd_i,'NAXIS3',1
        hrot, oldim[*,*,i], oldhd_i, newim_i, newhd_i, $
            angle, xc, yc,$
            int, MISSING=missing, $
            INTERP = interp, CUBIC = cubic, PIVOT = pivot,ERRMSG= errmsg
        newim[0,0,i]=newim_i
        newhd=newhd_i
    endfor
endif else begin
    hrot, oldim, oldhd, newim, newhd, $
        angle, xc, yc,$
        int, MISSING=missing, $
        INTERP = interp, CUBIC = cubic, PIVOT = pivot,ERRMSG= errmsg
endelse

END
