PRO PLTMOM_PV,prefix
;+
; NAME:
;   PLTMOM_PV
;
; PURPOSE:
;   summary plot for mom0 and mom0-pv maps
;
; INPUTS:
;   PREFIX    -- prefix for maps
;                e.g.'yourpath/n4254co.gm'
;
; OUTPUTS:
;   prefix.mom0pv.eps
;
; HISTORY:
;
;   20130214  RX  introduced
;   20150910  RX  merged from makemompv_pl.pro
;-

;   CHECK IMAGE FILES

im=prefix+'.mom0.fits'
eim=prefix+'.emom0.fits'
imxv=prefix+'.mom0xv.fits'
imvy=prefix+'.mom0vy.fits'

if  ~file_test(im) and $
    ~file_test(eim) and $
    ~file_test(imxv) and $
    ~file_test(imvy) then begin
    print,'no mom0/emom0/mom0xy/mom0vy image found'
    return
endif

im=readfits(im,imhd)
eim=READFITS(eim,eimhd)
imxv=readfits(imxv,imxvhd)
imvy=readfits(imvy,imvyhd)

sz_xv=size(imxv,/d)
sz_vy=size(imvy,/d)
sz=size(im,/d)
scy=sz_xv[1]*1.0/sz[1]
scx=sz_vy[0]*1.0/sz[0]

nvels=(size(imvy))[1]
vels=sxpar(imvyhd,'CRVAL1')+(findgen(nvels)+1-sxpar(imvyhd,'CRPIX1'))*sxpar(imvyhd,'CDELT1')
vels=vels/1000.

set_plot, 'ps'
device, filename=prefix+'.mom0pv.eps', $
    bits_per_pixel=8,/encapsulated,$
    xsize=10,ysize=10.*(sz[1]+sz_xv[1])/(sz[0]+sz_vy[0]),/inches,/col,xoffset=0,yoffset=0

!p.thick = 1.5
!x.thick = 1.5
!y.thick = 1.5
!z.thick = 1.5
!p.charsize=1.2
!p.charthick=1.5
xyouts,'!6'

; MOM-0 (XY) PLOT 
pos=[0.1,0.1,0.1+0.88/(scx+1.0),0.1+0.88/(scy+1.0)]

cgloadct,13
cgimage,im,pos=pos,stretch=1,/noe
imcontour,im,imhd,nlevels=10,$
    /noe,pos=pos,/nodata,color='red',AXISCOLOR='red',subtitle=' '
imcontour,im,imhd,nlevels=10,$
    /noe,pos=pos,/nodata,color='black',AXISCOLOR='black',subtitle=' ',ticklen=0
if  (where(eim ne eim))[0] ne -1 then begin
    sz=size(eim,/d)
    bb=find_boundary(where(eim eq eim),xsize=sz[0],ysize=sz[1])
    oplot,bb[0,*],bb[1,*],color=cgcolor('yellow'),linestyle=2
endif
    

subpos_xv=[pos[0],pos[3]+0.01,pos[2],pos[3]+0.9/(scx+1.0)*scx]
cgloadct,13
cgimage,imxv,pos=subpos_xv,stretch=1,/noe
cgloadct,0
plot,[0,1],[min(vels),max(vels)],/nodata,/noe,pos=subpos_xv,$
    xstyle=5,ystyle=1,$
    ytitle='Velocity [km/s]',color=cgcolor('red')
plot,[0,1],[min(vels),max(vels)],/nodata,/noe,pos=subpos_xv,$
    xstyle=5,ystyle=1,$
    ytitle='Velocity [km/s]',ticklen=0.0   
imcontour,im,imhd,nlevels=10,$
    /noe,pos=subpos_xv,/nodata,color='red',AXISCOLOR='red',ystyle=5,xtickformat='(A1)',$
    subtitle=' ',xticklen=!p.ticklen*(0.6/0.25),xtitle=' ',ytitle=' '
    
subpos_yv=[pos[2]+0.01,pos[1],pos[2]+0.9/(scy+1.0)*scy,pos[3]]
cgloadct,13
cgimage,imvy,pos=subpos_yv,stretch=1,/noe;,minvalue=0.0
cgloadct,0
plot,[min(vels),max(vels)],[0,1],/nodata,/noe,pos=subpos_yv,$
    xstyle=1,ystyle=5,$
    xtitle='Velocity [km/s]',color=cgcolor('red')
plot,[min(vels),max(vels)],[0,1],/nodata,/noe,pos=subpos_yv,$
    xstyle=1,ystyle=5,$
    xtitle='Velocity [km/s]',ticklen=0.0 
imcontour,im,imhd,nlevels=10,$
    /noe,pos=subpos_yv,/nodata,color='red',AXISCOLOR='red',xstyle=5,ytickformat='(A1)',$
    subtitle=' ',yticklen=!p.ticklen*(0.6/0.25),xtitle=' ',ytitle=' '

device, /close
set_plot,'X'

END
