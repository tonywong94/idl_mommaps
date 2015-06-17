FUNCTION maskguard_3d, mask, guard = guard
  
;+
;
; NAME:
;   maskguard_3d
;
; PURPOSE:
;
;   Read a binary cube mask and expand it by specified increment
;
; CALLING SEQUENCE:
;
;   maskguard_3d, mask, guard = [3, 3, 1]
;
; INPUTS:
;
;   mask : the 3D mask array.
;
;   guard=[xshift, yshift, vshift]: integers describing the number of voxels along
;   each axis that the mask should be expanded.
;
; KEYWORD PARAMETERS:
;
;
; MODIFICATION HISTORY:
;
;    Shamelessly hacked from fits2props and dilate_mask routines by Erik Rosolowsky.
;    Annie H, 1 March 2011
;    tw  20150601 - modified for use in mommaps
;    tw  20150616 - use padding to avoid wraparound
;

  if  n_elements(guard) eq 0 then begin
    guard = [0,0,0]
  endif else if  n_elements(guard) eq 1 then begin
    guard = [0,0,guard[0]]
  endif else if  n_elements(guard) eq 2 then begin
    guard = [guard[0],guard[1],0]
  endif

  mask = padding(mask, max(guard))
  expmask = mask*0.0

  print, 'Expanding mask by (x,y,v) = (',guard[0],guard[1],guard[2],') pixels'

; Grow the mask by specified number of pixels along each axis
  for i=0, guard[0] do begin
      for j = 0, guard[1] do begin 
          for k = 0, guard[2] do begin 
            ;print, i, j, k   
            expmask = (expmask+shift(mask,i,j,k)+shift(mask,-i,-j,-k)+shift(mask,i,-j,k)+shift(mask,-i,j,k)+shift(mask,i,j,-k)+shift(mask,i,-j,-k)+shift(mask,-i,j,-k)+shift(mask,-i,-j,k))
            ;print,expmask
            ;print, " "
          endfor
      endfor
  endfor 
 
  expmask = padding(expmask gt 0,-max(guard))
  return, expmask

end
