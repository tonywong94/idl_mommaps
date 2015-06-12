pro radiohead, hd, structure = structure, ms = ms

;+
; NAME:
;    radiohead
; PURPOSE:
;    To extract vital information from the header of a FITS file and
;    return it in a structure.  Designed for Interferometer Data.
;
; INPUTS:
;   HD - A FITS header
;
; KEYWORD PARAMETERS:
;   STRUCTURE - The name of a structure in which the header
;               information is returned.
;   MS - Set this keyword to return velocities in m/s.  Defaults to km/s.
;
; OUTPUTS:
;
;
; MODIFICATION HISTORY:
;
;       20150526  TW  based on RDHD.PRO from Erik Rosolowsky.
;       20130204  RX
;       correctly read the rest frequency in the new FITS convention
;       change the file name to rd_hd.pro
;       improve the acurracy of jypb2k
;-

; Basic dimensions
  dim = SXPAR(hd, 'NAXIS')
  naxis1 = SXPAR(hd, 'NAXIS1')
  naxis2 = SXPAR(hd, 'NAXIS2')
  naxis3 = SXPAR(hd, 'NAXIS3')
  
; Beam information
  bunit  = SXPAR(hd, 'BUNIT')
  bm_maj = SXPAR(hd, 'BMAJ')*3600. > SXPAR(hd, 'BMMAX')*3600.
  bm_min = SXPAR(hd, 'BMIN')*3600. > SXPAR(hd, 'BMMIN')*3600.
  bpa    = SXPAR(hd, 'BPA')

; AIPS beam extraction bits.
  if bm_maj eq 0 then begin
    ind = where(strpos(hd, 'CLEAN') gt 0 and $
                strpos(hd, 'BMAJ') gt 0 and $
                strpos(hd, 'AIPS') gt 0,  ct)
    if ct gt 0 then begin
      str = hd[max(ind)]
      str = strsplit(str, /ext)
      bm_maj = float(str[4])*3600
      bm_min = float(str[6])*3600
      bpa = float(str[8])
    endif
  endif
  
; Rest frequency and observation date
  freq = SXPAR(hd, 'RESTFREQ') > SXPAR(hd, 'RESTFRQ') > SXPAR(hd, 'FREQ0')
  date = SXPAR(hd, 'DATE-OBS')
  EXTAST, hd, astrom

; Radio astronomy parameters
  lam_mm=!VALUES.F_NAN
  jypb2k=!VALUES.F_NAN
  k2jypb=!VALUES.F_NAN
  if freq ne 0. then begin
    lam_mm = 1.e3*!CONST.C/freq
    ; Conversion Factors for the map to change from Jy/beam to K.
    if bm_maj ne 0 then begin
      jypb2k = 13.58491*(lam_mm/(sqrt(bm_maj*bm_min)))^2
      k2jypb = 1./jypb2k
    endif
  endif

; Get typical values of RA and DEC
  if dim eq 1 then begin
    ravec = SXPAR(hd, 'CRVAL2')
    decvec = SXPAR(hd, 'CRVAL3')
  endif else begin
;    XY2AD, astrom.crpix[0], indgen(naxis2), astrom, null, decvec
;    XY2AD, indgen(naxis1), astrom.crpix[1], astrom, ravec, null
    XY2AD, naxis1/2, indgen(naxis2), astrom, null, decvec
    XY2AD, indgen(naxis1), naxis2/2, astrom, ravec, null
  endelse

; For datacubes
  if naxis3 gt 1 then begin
    astrom2 = {CD:astrom.cd, $
               CDELT:[astrom.cdelt, SXPAR(hd, 'CDELT3')], $
               CRPIX:[astrom.crpix, SXPAR(hd, 'CRPIX3')], $
               CRVAL:[astrom.crval, SXPAR(hd, 'CRVAL3')], $
               CTYPE:[astrom.ctype, SXPAR(hd, 'CTYPE3')], $
               LONGPOLE:astrom.longpole, $
               LATPOLE:astrom.latpole, $
               PV2:astrom.pv2}
    astrom = astrom2
    velvec = (findgen(naxis3)+1-astrom.crpix[2])*$
             astrom.cdelt[2]+astrom.crval[2]
  endif else begin
    velvec = 0 
    dim = 2
  endelse

; Output data structure (uses km/s for velocity by default)
  if not keyword_set(ms) then velvec = velvec/1000.

  ppbeam = abs((bm_maj*bm_min/3600.^2)/(astrom.cdelt[0]*astrom.cdelt[1])*$
               2.*!dpi/(8.*alog(2.)))
  if  SXPAR(hd,'CTYPE1') eq 'RA---CAR' then $ 
  ppbeam = abs((bm_maj*bm_min/3600.^2)/(astrom.cdelt[1]*astrom.cdelt[1])*$
    2.*!dpi/(8.*alog(2.)))                 
  specsys=SXPAR(hd, 'SPECSYS')
  velref=SXPAR(hd, 'VELREF')
  structure = {naxis1:naxis1, naxis2:naxis2, $
                 naxis3:naxis3, ra:ravec, dec:decvec, v:velvec, $
                 ctype:astrom.ctype, $
                 crpix:astrom.crpix, $
                 crval:astrom.crval, $
                 cdelt:astrom.cdelt, $
                 k2jypb:k2jypb, jypb2k:jypb2k, freq:freq, $
                 bmaj:bm_maj, bmin:bm_min, date:date, bpa:bpa, $
                 ppbeam:ppbeam,specsys:specsys,$
                 velref:velref}

  return
end


