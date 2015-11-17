PRO PLTMOM,prefix,scale=scale,label=label
;+
; NAME:
;   PLTMOM
;
; PURPOSE:
;   summary plot for mom0 and mom1 maps
;
; INPUTS:
;   PREFIX    -- prefix for maps
;                e.g.'yourpath/n4254co.gm'
;   
; OUTPUTS:
;   prefix.mom01.eps
;
; HISTORY:
;
;   20130214  RX  introduced
;   20130910  RX  changed to absolute coordinate
;   20150910  RX  merged from makemom_pl.pro
;                   clean up dependency
;                   correct a bug related to BPA
;-

;   CHECK IMAGE FILES

if  ~file_test(prefix+'.mom0.fits') and $
    ~file_test(prefix+'.emom0.fits') and $
    ~file_test(prefix+'.mom1.fits') then begin
    print,'no mom0/emom0/mom1 image found'
    return
endif

set_plot, 'ps'
device, filename=prefix+'.mom01.eps', $
    bits_per_pixel=8,/encapsulated,$
    xsize=8,ysize=4.8,/inches,/col,xoffset=0,yoffset=0

!p.thick = 1.5
!x.thick = 1.5
!y.thick = 1.5
!z.thick = 1.5
!p.charsize=1.0
!p.charthick=1.5
xyouts,'!6'

eim=prefix+'.emom'+strtrim(0,2)+'.fits'
eim=READFITS(eim,eimhd)

for i=0,1 do begin
  
    im=prefix+'.mom'+strtrim(i,2)+'.fits'
    im=READFITS(im,imhd)
    
    if n_elements(scale) ne 0 and i eq 0 then im=im*scale
    pos=[0.08,0.1,0.52,0.85]+i*[0.46,0.0,0.46,0.0]
    
    ;   COLOR SCALING
    if  i eq 0 then begin
        CGLOADCT,13
        minvalue=min(im,/nan)
        maxvalue=max(im,/nan)
        if  minvalue eq maxvalue then begin
            minvalue=-1 & maxvalue=1
        endif
        CGIMAGE,im,pos=pos,stretch=1,/noe,/KEEP_ASPECT_RATIO,minvalue=minvalue,maxvalue=maxvalue
    endif
    if  i eq 1 then begin
        CGLOADCT,13
        CGIMAGE,im,pos=pos,stretch=1,/noe,/KEEP_ASPECT_RATIO
    endif
  
    RADIOHEAD,imhd,s=s
    psize=abs(s.cdelt[0])*3600
    sz=size(im,/d)

    CGLOADCT,0
    if  i eq 1 then ytickformat='(A1)' else ytickformat=''
    imcontour,im,imhd,$
        /noe,/nodata,pos=pos,$
        title=dataid,/overlay,$
        ytickformat=ytickformat,$
        subtitle=' ',$
        color='black',AXISCOLOR='black'
    imcontour,im,imhd,$
        /noe,/nodata,pos=pos,$
        title=dataid,/overlay,$
        ytickformat='(A1)',xtickformat='(A1)',$
        xtitle=' ',ytitle=' ',$
        subtitle=' ',$
        color='black',AXISCOLOR='red'
    if  (where(eim ne eim))[0] ne -1 then begin   
        bb=find_boundary(where(eim eq eim),xsize=sz[0],ysize=sz[1])
        oplot,bb[0,*],bb[1,*],color=cgcolor('yellow'),linestyle=2
    endif
    
    pos=[pos[0],pos[3],pos[2],pos[3]]+[0.0,0.06,0.0,0.10]
    cgloadct,13
    if  i eq 0 then title='Intensity ['+strtrim(sxpar(imhd,'BUNIT'),2)+']'
    if i eq 1 then title='Velocity Field ['+strtrim(sxpar(imhd,'BUNIT'),2)+']'
    crange=[min(im,/nan),max(im,/nan)]
    tickint=fix(crange[1]-crange[0])/5
    if crange[0] eq crange[1] or (where(im eq im))[0] eq -1 then crange=[-1.,1.]
    cgCOLORBAR, range=crange, POSITION=pos,title=title,tlocation='TOP',tickinterval=tickint
    
    cgloadct,0
    tvellipse,s.bmaj/2.0/psize,s.bmin/2.0/psize,$
        sz[0]/10.0,sz[1]/10.0,$
        s.bpa-90.0,$
        /data,noclip=0,color=cgcolor('cyan'),/fill
    
    if  keyword_set(label) then al_legend,label,/top,/right,textcolor='yellow',box=0
    
endfor

device, /close
set_plot,'X'


END
