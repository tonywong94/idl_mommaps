;+
; NAME:
;	CONVOL3D
; PURPOSE:
;	Convolution or correlation of two 3D arrays (volumetric data),
;	or autocorrelation of a single 3D array.
;	Default is to compute using product of Fourier transforms.
;
; CALLING SEQUENCE:
;
;	imconv = convol3d( voldata, psf, FT_PSF = psf_FT )
;  or:
;	correl = convol3d( voldata1, voldata2, /CORREL )
;  or:
;	correl = convol3d( voldata, /AUTO )
;
; INPUTS:
;	voldata = 3-D array to be convolved with PSF.
;	psf = the 3-D array Point Spread Function,
;		with size < or = to size of voldata.
;
; KEYWORDS:
;
;	FT_PSF = passes out/in the Fourier transform of the PSF,
;		(so that it can be re-used the next time function is called).
;
;	FT_VOLDATA = passes out/in the Fourier transform of voldata.
;
;	/CORRELATE uses the conjugate of the Fourier transform of PSF,
;		to compute the cross-correlation of voldata and PSF,
;		(equivalent to IDL function convol() with NO rotation of PSF).
;
;	/AUTO_CORR computes the auto-correlation function of voldata using FFT.
;
;	/NO_FT overrides the use of FFT, using IDL function convol() instead.
;		(then PSF is rotated by 180 degrees to give same result)
; METHOD:
;	When using FFT, PSF is centered & expanded to size of voldata.
;	Basically a copy of function convolve modified from 2D to 3D.
; HISTORY:
;	written, Frank Varosi, NASA/GSFC 1997,
;-

function convol3d, voldata, psf, FT_PSF=psf_FT, FT_VOLDATA=vdFT, NO_FT=noft, $
			CORRELATE=correlate, AUTO_CORRELATION=auto

	sv = size( voldata )
	sp = size( psf )
	spf = size( psf_FT )
	if (spf(0) ne 3) or (spf(spf(0)+1) ne 6) then no_psf_FT = 1

	if (spf(0) eq 3) and (sv(0) eq 3) then begin
		if max( abs( spf(1:sp(0)) - sv(1:sv(0)) ) ) then no_psf_FT = 1
	   endif

	if (sv(0) ne 3) or ((sp(0) ne 3) and $
	    keyword_set( no_psf_FT ) and (NOT keyword_set( auto ))) then begin
		print,"syntax:	result = convol3d( voldata, psf )
		print,"    or:	result = convol3d( voldata, psf, FT_PSF=psf_FT )
		print,"    or:	correl = convol3d( voldata1, voldata2, /CORREL )
		print,"    or:	autocorr = convol3d( voldata, /AUTO )
		return,0
	   endif

	if keyword_set( noft ) then begin
		if keyword_set( auto ) then begin
			message,"auto-correlation available only with FFT",/INFO
			return, voldata
		  endif else if keyword_set( correlate ) then $
				return, convol( voldata, psf ) $
			else	return, convol( voldata, rotate( psf, 2 ) )
	   endif

	sc = sv/2
	npix = N_elements( voldata )
	sif = size( vdFT )

	if (sif(0) ne 3) or (sif(sif(0)+1) ne 6) or $
	   (sif(1) ne sv(1)) or (sif(2) ne sv(2)) then vdFT = FFT( voldata,-1 )

	if keyword_set( auto ) then $
	 return, shift( npix*float( FFT( vdFT*conj( vdFT ), 1 ) ), $
						sc(1), sc(2), sc(3) )

	if Keyword_Set( no_psf_FT ) then begin
		s2 = sc + (sv MOD 2)*(sp LT sv)   ;correct for odd size voldata.
		Loc = (s2 - sp/2) > 0		;center PSF in new array,
		s = (sp/2 - s2) > 0	;handle all cases: smaller or bigger
		L = (s + sv-1) < (sp-1)
		psf_FT = dcomplexarr( sv(1), sv(2), sv(3) )
		psf_FT(Loc(1),Loc(2),Loc(3)) =psf(s(1):L(1),s(2):L(2),s(3):L(3))
		psf_FT = FFT( psf_FT, -1, /OVERWRITE )
	   endif

	if keyword_set( correlate ) then $
		conv = npix * float( FFT( vdFT * conj( psf_FT ), 1 ) ) $
	  else	conv = npix * float( FFT( vdFT * psf_FT, 1 ) )

return, shift( conv, sc(1), sc(2), sc(3) )
end