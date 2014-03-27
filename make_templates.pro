; -----------------------------
; SEPARATE NORMALIZATION SCRIPT
; -----------------------------
function kellel_normalize, lam, flx, norm_fact=norm_fact, rng=rng

lamn = lam
flxn = flx 
if (n_elements(rng) ne 0) then begin
 w = where(lam ge rng(0) and lam le rng(1),cnt)
 if (cnt gt 0) then begin
  lamn = lam(w)
  flxn = flx(w)
 endif
endif

norm_fact = MEAN(flxn,/nan)
;norm_fact = int_tabulated(lamn,flxn)

if ~finite(norm_fact) then MESSAGE, 'problem with normalization'
;print, 'norm', fact

return, flx/norm_fact
end

; -----------------------------
; SCRIPT TO CREATE TEMPLATES
; -----------------------------
function kellel_template, lam, flxs, flx_uncs, norm_facts=norm_facts, flags_in=flags_in, flags_out=flags_out, flags_one=flags_one,sigma=sigma, sd=sd, chi2=chi2, maxs_templ=maxs_templ, mins_templ=mins_templ, no_renorm=no_renorm

;JUST MAKE THE TEMPLATE AND FIND OUTLIERS

nspec = n_elements(flxs(0,*))
nlam = n_elements(flxs(*,0))
norm_facts = fltarr(nspec) ;set up normalivation array
flxns = flxs * 0. ;set up array
flxn_uncs = flxs * 0.
flags_out = flags_in ;initialize flags
flags_one = intarr(nspec) ;just flags from this run
if (n_elements(sigma) eq 0) then begin 
	sigma = 2.0
	message, 'using chi^2 < 2.0',/info
endif
;if (n_elements(flag) eq 0) then flag = fltarr(nspec)*0

wselect = where(flags_in eq 0,cntselect)
if (cntselect eq 0) then message, 'No spectra available for creating template'
;print, 'wselect: ', cntselect

; first normalization of each spectrum

for j=0,nspec-1 do begin
  flxns[*,j] = kellel_normalize(lam, flxs[*,j], norm_fact=norm_fact)
  flxn_uncs[*,j] = flx_uncs[*,j] / norm_fact
  norm_facts[j] = norm_fact
  ;if j eq 0 then plot, lam, flxns(*,j) else oplot, lam, flxns(*,j)
endfor

;COMBINE THE SELECTED normalized spectra using the uncertainties
mc_meancomb2, flxns[*,wselect], template_first, datavar=flxn_uncs[*,wselect]

; repeat normalization against template
if ~keyword_set(no_renorm) then begin
	for j=0,nspec-1 do begin
		norm_facts[j] = median( flxs[*,j] / template_first );
		flxns[*,j] = flxs[*,j] / norm_facts[j]
		flxn_uncs[*,j] = flx_uncs[*,j] / norm_facts[j]
		;print, 'RENORMIALIZED TO TEMPLATE'
		;print,min(flxns[*,j]),max(flxns[*,j])
	endfor
endif

;recalculate template using re-normalized spectra
mc_meancomb2, flxns[*,wselect], template, datavar=flxn_uncs(*,wselect)

;CHARACTERIZE NEW TEMPLATE
mins_templ = MIN(flxns[*,wselect],dim=2)
maxs_templ = MAX(flxns[*,wselect],dim=2)

; standard deviation of spectra used in template at each wavelength
; unweighted
sd = fltarr(nlam)
for j = 0, nlam-1 do sd[j] = stddev(flxns[j,wselect])

; identify outliers of new template
; JUST FLAG (REJECT LATER) based on reduced chi square goodness of fit test
chi2 = fltarr(nspec)
i_NotNANs_template = where(finite(template) eq 1,n_template)
for ispec=0,nspec-1 do begin
  chi2[ispec] = total( (flxns[*,ispec]-template)^2/sd^2, /NAN ) / (n_template-1.)
  if (chi2[ispec] gt sigma) then begin
	  flags_one[ispec] = flags_one[ispec]+1
	  flags_out[ispec] = flags_in[ispec]+1
	  ;print,'rejected chi2: ', chi2[ispec]
  endif
endfor

;flag = flag<1 ; sets all flags greater than 1 to 1 and leaves 0 unchanged.

return, template
end


; -----------------------------
; MAIN ROUTINE
; -----------------------------
pro make_templates, type=type, btypes=btypes, reset=reset, sigma=sigma, ptype=ptype, batch=batch, ps=ps, interactive=interactive

; if N_params() lt 1 then begin
; 	print, n_params()
;     print, "Syntax -  print, 'kellel, type=type, [/RESET, sigma=, ptype, batch, /PS]' "
; 	print, "type = 2 = L2. default sigma=2, type=1,ptype=0"
; print, "ptype = 0 (default, all), 1 (kept+template), 2 (rejects+template)"
; print, "Do more than one spectral typo at once with /batch then type=[1,4]"
;     return
; endif
	

trix = [-1,0,1,-1]
triy = [-1,1,-1,-1]
crx = rndoff(10.*cos(findgen(41)*!pi/20.))
cry = rndoff(10.*sin(findgen(41)*!pi/20.))
crx = round(10.*cos(findgen(41)*!pi/20.))
cry = round(10.*sin(findgen(41)*!pi/20.))
sqx = [-1,-1,1,1,-1]
sqy = [-1,1,1,-1,-1]

!p.font=0
!p.thick=4
!x.thick=3
!y.thick=3
;tek_color
@colors_kc

color_rejected=[green,cyan,purple,red,dkgreen,orange,gray,magenta,pink,ltorange,ltgray,purple2,dkred,dkorange,dkgrey,dkpurple, green,cyan,purple,red,dkyellow,dkgreen,orange,gray]
color_kept=blue
colorfill_sdev = gray
colorfill_minmax = vltgray


If ~keyword_set(PS) then begin
	 device, Decomposed=0
endif

; fold_root='/Users/adam/papers/kellel/'
fold_root='/Users/kelle/Dropbox/Analysis/NIRtemplates/'
bfold = fold_root
; spfold = '/Users/adam/spectra/spex/prism/'
; filtfold = '/Users/adam/idl/filters/'
if (n_elements(sigma) eq 0) then sigma=2.
if (n_elements(type) eq 0) then type=1
if (n_elements(ptype) eq 0) then ptype=0

; -----------------------------
; BATCH PROCESS
; -----------------------------

if (keyword_set(batch)) then begin
  ;types = indgen(3)+1
  for i=0,n_elements(btypes)-1 do begin
   for j=0,2 do begin
    for k=1,2 do kellel, type=btypes(i), ptype=j, sigma=k*1. , ps=ps
   endfor
  endfor
  goto, FINISH
 endif
 
; -----------------------------
; EXAMINE VARIANCE IN SPECTRA, ATTEMPT TO DEFINE "NORMAL"
; -----------------------------

spt = 'L'+strmid(strtrim(string(type),2),0,1)
print, spt, '  chi^2: ',strn(sigma)
print, 'ptype: ', strn(ptype)
sfold = bfold+spt+'/'
spec_fold = '/Users/kelle/Dropbox/Data/nir_spectra_low/'
dfold = sfold+spt+'s/'
slist = dfold+spt+'s_in.txt'
ofold = sfold+'output_'+spt+'/'
strfile = sfold+'spectra.dat'
band = ['J','H','K']

if (file_search(sfold) eq '' or file_search(dfold) eq '') then begin
	print, 'cannot find '+sfold+' or '+dfold
	goto, FINISH
endif
if (file_search(ofold) eq '') then begin
	spawn, 'mkdir '+ofold
	print, 'created directory: ' + ofold
endif

if (file_search(strfile) eq '' or keyword_set(reset)) then begin
	MESSAGE,'re-selecting spectra',/info
	
	;use .txt file if present, otherwise use all spectra in folder
	if file_search(slist) eq '' then begin
		sfiles = file_search(dfold+'*.fits')
		spec_fold='' 
		MESSAGE,'using spectra in ' + dfold, /info
	endif else begin
		READCOL,slist,sfiles,format='A'
		sfiles=sfiles 
		MESSAGE,'using spectra listed in ' + slist,/info
	endelse
	
	snr = fltarr(n_elements(sfiles))
	for ifile=0,n_elements(sfiles)-1 do begin
	 fits = READFITS(spec_fold+sfiles(ifile),hd,/silent)
	 if (ifile eq 0) then begin
		 ;setup the arrays
	  lam = fits(*,0)
	  flx = fltarr(n_elements(lam),n_elements(sfiles))
	  flx_unc = flx*0.
	  ;wnorm = where(lam ge 1.2 and lam le 1.3)		; initial normalization
	 endif
	 ;interpolate onto wavelength of first object
	 flx(*,ifile) = interpol(fits(*,1),fits(*,0),lam)
	 ;mx = max(flx(wnorm,ifile))
	 ;flx(*,ifile) = flx(*,ifile)/mx
	 flx_unc(*,ifile) = interpol(fits[*,2],fits[*,0],lam)
	 ;snr(ifile) = max(smooth(flx(*,ifile)/ns(*,ifile),9))
	 ;plot, lam,flx(*,ifile)
	endfor

	str = create_struct($
	'files',sfiles,$
	'lam',lam,$
	'flx',flx,$
	'flx_unc',flx_unc,$
	'names', strrep(strrep(sfiles,dfold,''),'.fits',''))	
 	 
	 save, str, file=strfile

endif else begin
	restore, file=strfile
	MESSAGE,'restoring spectra from save file, not re-reading',/info
endelse

nspec=n_elements(str.files)
print,'nspec=',nspec

fbase = 'comp_'+spt+'_s'+strmid(strtrim(string(sigma),2),0,3)

if keyword_set(ps) then begin	
	set_plot, 'ps'
	device, /encapsulated, ysize=18, xsize=24, filename=ofold+fbase+'_v'+strtrim(string(ptype+1),2)+'.eps', /portrait, bits_per_pixel=8, /color
	print, 'Creating: '+ ofold+fbase+'_v'+strtrim(string(ptype+1),2)+'.eps'
endif else begin
	;window,1
endelse

; iterative scan, normalizing in individual bands
lam_norm = [[0.87,1.39],[1.41,1.89],[1.91,2.39]]
max_loops = 10

flags_in = intarr(nspec)
flags_j = intarr(nspec)
flags_h = intarr(nspec)
flags_k = intarr(nspec)
chi2_j = fltarr(max_loops)
chi2_h = fltarr(max_loops)
chi2_k = fltarr(max_loops)
chi2all = fltarr(nspec,4)
i_spec = indgen(nspec)
iloop=0
loop_stop=0
loop_stop_j=0
loop_stop_h=0
loop_stop_k=0
;Make template 10 times to iterativly reject


while loop_stop eq 0 do begin
;for iloop=0,nloops-1 do begin

	;Loop over bands
	;Reject any spectrum that has chi^2 > 2 compared to template
	for mmm=0,n_elements(band)-1 do begin

		w = where(str.lam ge lam_norm(0,mmm) and str.lam le lam_norm(1,mmm),cnt)
		
		;Always re-normalize
		no_renorm = 0
		;renorm for the first half of loops, but not the second
		;if iloop ge nloops/2.  then no_renorm = 1 else no_renorm = 0
		;print, iloop,no_renorm
		
		;takes flags_in and modifies to reflect new rejects
		template = kellel_template(str.lam(w), str.flx(w,*), str.flx_unc(w,*), norm_facts=norm_facts, flags_in=flags_in, flags_out=flags_out, flags_one=flags_one,sd=sd, sigma=sigma, chi2=chi2, mins_templ=mins_templ, maxs_templ=maxs_templ,no_renorm=no_renorm)
		
		i_keep = where(flags_out eq 0,cnt_keep)
		i_rejects = where(flags_out ge 1,cnt_rejects)
				
		;keep track of band-by-band flags
		case mmm of 
			0: BEGIN	
				flags_j = flags_j + flags_one
				chi2_j[iloop] = mean(chi2[i_keep])
				if iloop ge 1 then if chi2_j[iloop] eq chi2_j[iloop-1] then loop_stop_j=1			
			END
			1: BEGIN
				flags_h = flags_h + flags_one
				chi2_h[iloop] = mean(chi2[i_keep])	
				if iloop ge 1 then if chi2_h[iloop] eq chi2_h[iloop-1] then loop_stop_h=1			
			END
			2: BEGIN
				flags_k = flags_k + flags_one
				chi2_k[iloop] = mean(chi2[i_keep])	
				if iloop ge 1 then if chi2_k[iloop] eq chi2_k[iloop-1] then loop_stop_k=1			
			END
		endcase
			
		if loop_stop_j eq 1 and loop_stop_h eq 1 and loop_stop_k eq 1 then loop_stop =1 else loop_stop =0
			
		i_rejects_j = where(flags_j ge 1,cnt_rejects_j)
		i_rejects_h = where(flags_h ge 1,cnt_rejects_h)
		i_rejects_k = where(flags_k ge 1,cnt_rejects_k)
		
		if cnt_rejects_j ge 1 then ii_rejects_j = intarr(cnt_rejects_j)
		if cnt_rejects_h ge 1 then ii_rejects_h = intarr(cnt_rejects_h)
		if cnt_rejects_k ge 1 then ii_rejects_k = intarr(cnt_rejects_k)
		for i=0,cnt_rejects_j-1 do ii_rejects_j[i] = where(i_rejects eq i_rejects_j[i])
		for i=0,cnt_rejects_h-1 do ii_rejects_h[i] = where(i_rejects eq i_rejects_h[i])
		for i=0,cnt_rejects_k-1 do ii_rejects_k[i] = where(i_rejects eq i_rejects_k[i])
	 
	 	if mmm eq 2 then print,iloop,cnt_rejects_j,cnt_rejects_h,cnt_rejects_k
	 
		if (cnt_keep eq 0) then message, 'Warning, rejected all sources!'
		
		;print, 'loop & mean chi^2 = ' + strn(iloop) + ' ' + strn(mean(chi2[i_keep]))
		;print, 'nkept, n_rejected = ', cnt_keep, cnt_reject
		;print, 'n_rejects , j ,h,k = ', cnt_reject_j, cnt_reject_h, cnt_reject_k
		;
		;if mmm eq 2 then print, ' '
			
		; PLOT	chi^sq for each band	
		if ~keyword_Set(ps) then begin
			if mmm eq 0 then begin
				window, 0
				!p.multi=[0,3,1]
				!p.thick=2
			endif else begin
				wset, 0
				!p.multi=[3-mmm,3,1]
			endelse			
			plot, i_spec, chi2, psym=3, title=strn(mmm)	, xr=[-1,nspec],xstyle=1,symsize=2,/ylog
			if cnt_rejects ge 1 then oplot, i_spec[i_rejects], chi2[i_rejects], psym=6, symsize=3.5,color=red
			
			if cnt_rejects_j ge 1 then for j=0,cnt_rejects_j-1 do oplot, [i_spec[i_rejects_j[j]]], [chi2[i_rejects_j[j]]], psym=5, symsize=3,color=color_rejected[ii_rejects_j[j]]
			if cnt_rejects_h ge 1 then for j=0,cnt_rejects_h-1 do oplot, [i_spec[i_rejects_h[j]]], [chi2[i_rejects_h[j]]], psym=4, symsize=3,color=color_rejected[ii_rejects_h[j]]
			if cnt_rejects_k ge 1 then for j=0,cnt_rejects_k-1 do oplot, [i_spec[i_rejects_k[j]]], [chi2[i_rejects_k[j]]], psym=4, symsize=3,color=color_rejected[ii_rejects_k[j]]
			
			;oplot, i_spec[i_rejects_j], chi2[i_rejects_j], psym=5, symsize=3,color=cyan
			;if i_rejects_h[0] ne -1 then oplot, i_spec[i_rejects_h], chi2[i_rejects_h], psym=4, symsize=2,color=green
			;if cnt_rejects_k[0] ge 1 then oplot, i_spec[i_rejects_k], chi2[i_rejects_k], psym=1, symsize=2,color=magenta
		
			plots, [0,nspec-1], [mean(chi2[i_keep]), mean(chi2[i_keep])]
			plots, [0,nspec-1], [sigma,sigma], linestyle=1
			xyouts, mmm/3.0, 0, 'mean chi^2 =' + strn(mean(chi2[i_keep])), /normal
		ENDIF
				
		flags_in = flags_out
		;flags = flag+f
		chi2all(*,mmm) = chi2
		 
		if (iloop ge 1) then begin
		   
		   ;SET UP THE PLOTS
			if mmm eq 0 then begin
				   if ~keyword_set(ps) then window,1
				   !p.multi=[0,3,2]
			endif else begin
				   if ~keyword_set(ps) then wset,1
				   !p.multi=[6-mmm,3,2]
			endelse
	  	 	
			flux_mins = fltarr(nspec) & flux_maxs=fltarr(nspec)
	   		for i = 0,nspec-1 do flux_mins[i] = MIN( str.flx[w,i] / norm_facts[i] ) 
	   		for i = 0,nspec-1 do flux_maxs[i] = MAX( str.flx[w,i] / norm_facts[i] )
	    	yra = [MIN(flux_mins),MAX(flux_maxs) ] 

		   ;plot template
		   plot, str.lam(w), template, /xsty, yra=yra,/ysty, xtitle='!3Wavelength (!9m!3m)', ytitle='Normalized Flux', charsize=2, xmargin=[8,2], title=tit
		   			
		   case 1 of 
		    	ptype eq 0: begin
					;plot each object rejected from template
		     	   if (cnt_rejects gt 0) then for j=0,cnt_rejects-1 do oplot, str.lam(w), str.flx(w,i_rejects[j])/norm_facts[i_rejects[j]], $
					   color=color_rejected[j], thick=1
  		   		   ;polyfill the min-max of the spectra used in the template
  		   		   polyfill, [str.lam[w],reverse(str.lam[w])], [maxs_templ,REVERSE(mins_templ)], color=colorfill_minmax, /fill		   
				   ;polyfill the standard dev
  		   		   polyfill, [str.lam[w],reverse(str.lam[w])], [template+sd,reverse(template-sd)], color=colorfill_sdev, /fill
				   ;print out the names of rejected objects
					case mmm of
					   	0: if (cnt_rejects_j gt 0) then for j=0,cnt_rejects_j-1 do xyouts, 0.15+0.33*(mmm), 0.7-j*0.02, str.names[i_rejects_j[j]],color=color_rejected[ii_rejects_j[j]],size=0.8,/normal
						1: if (cnt_rejects_h gt 0) then for j=0,cnt_rejects_h-1 do xyouts, 0.15+0.33*(mmm), 0.7-j*0.02, str.names[i_rejects_h[j]],color=color_rejected[ii_rejects_h[j]],size=0.8,/normal
						2: if (cnt_rejects_k gt 0) then for j=0,cnt_rejects_k-1 do xyouts, 0.15+0.33*(mmm), 0.7-j*0.02, str.names[i_rejects_k[j]],color=color_rejected[ii_rejects_k[j]],size=0.8,/normal
					ENDCASE
		    	end
		    	ptype eq 1: begin ;only plot objects selected & the std dev
		     	   if (cnt_keep gt 0) then for j=0,cnt_keep-1 do oplot, str.lam[w], str.flx(w,i_keep[j])/norm_facts[i_keep[j]], color=color_kept, thick=1
				   ;polyfill the standard dev of the template
  		   		   polyfill, [str.lam[w],reverse(str.lam[w])], [template+sd,reverse(template-sd)], color=colorfill_sdev, /fill
		    	end
		    	else: begin ;plot everything
  		   		   if (cnt_keep gt 0) then for j=0,cnt_keep-1 do oplot, str.lam[w], str.flx(w,i_keep[j])/norm_facts[i_keep[j]], color=color_kept, thick=1
				   if (cnt_rejects gt 0) then for j=0,cnt_rejects-1 do oplot, str.lam(w), str.flx(w,i_rejects[j])/norm_facts(i_rejects[j]), color=color_rejected[j], thick=1
				case mmm of
				   	0: if (cnt_rejects_j gt 0) then for j=0,cnt_rejects_j-1 do xyouts, 0.15+0.33*(mmm), 0.7-j*0.02, str.names[i_rejects_j[j]],color=color_rejected[ii_rejects_j[j]],size=0.8,/normal
					1: if (cnt_rejects_h gt 0) then for j=0,cnt_rejects_h-1 do xyouts, 0.15+0.33*(mmm), 0.7-j*0.02, str.names[i_rejects_h[j]],color=color_rejected[ii_rejects_h[j]],size=0.8,/normal
					2: if (cnt_rejects_k gt 0) then for j=0,cnt_rejects_k-1 do xyouts, 0.15+0.33*(mmm), 0.7-j*0.02, str.names[i_rejects_k[j]],color=color_rejected[ii_rejects_k[j]],size=0.8,/normal
				ENDCASE
		    	end
		   endcase
		   
		   ;polyfill the min-max
		   ;polyfill, [str.lam(w),reverse(str.lam(w))], [maxs_templ,REVERSE(mins_templ)], color=colorfill_minmax, /fill
		   
		   ;polyfill the standard dev
		   ;polyfill, [str.lam(w),reverse(str.lam(w))], [template+sd,reverse(template-sd)], color=colorfill_sdev, /fill

		   oplot, str.lam(w), template, thick=2
		 		  
		   ;if (iloop eq nloops-1) then begin
			   if loop_stop eq 1 then begin
		   ; save template
		   dat = [[str.lam(w)],[template],[sd],[mins_templ],[maxs_templ]]
		   fxhmake,hdr,dat
		   output_fits = ofold+spt+'_template_'+band(mmm)+'.fits'
		   writefits, output_fits , dat, hdr
		   message, 'wrote' + output_fits, /info 
  	 	endif
	endif

		;if  (iloop eq nloops-1) then stop

	 endfor ;end loop over 3 (JHK) bands, mmm

 	if keyword_set(interactive) and iloop ge 1 then begin & print,'press any key to continue' & tmp=GET_KBRD() & endif
		iloop++
endwhile	
;endfor ;end iterative 10 loops, iloop

nloops=iloop
i_loop = indgen(nloops)

if ~keyword_set(ps) then begin
	window,3
	!p.multi=0
	plot, i_loop, i_loop, yr =[0.5,1.1], xr=[-1,nloops+1],/nodata,xstyle=1
	plots,[0,nloops],[1,1]
	oplot, i_loop, chi2_j[0:nloops-1], psym = 5, color=cyan, symsize=2
	oplot, i_loop, chi2_h[0:nloops-1], psym = 4, color=green, symsize=2
	oplot, i_loop, chi2_k[0:nloops-1], psym = 1, color=magenta, symsize=2
	;xyouts,0.5,0.5,/normal
	
ENDIF

;WRITE TXT FILES
 tb='	'
 openw, unit, ofold+fbase+'_band_rejects.txt', /get_lun
 printf, unit, '# '+strtrim(string(cnt_keep),2)+' '+spt+' dwarfs rejected from template construction because '
 printf, unit, '# reduced chi2 > '+strtrim(string(sigma),2)+' in any band (but not full spectrum)'
 printf, unit, '# Name'+tb+'J chi2'+tb+'H chi2'+tb+'K chi2'
 if (cnt_rejects gt 0) then for i=0,cnt_rejects-1 do printf, unit, str.names(i_rejects[i])+' '+tb+strtrim(string(chi2all(i_rejects[i],0)),2)+tb+strtrim(string(chi2all(i_rejects[i],1)),2)+tb+strtrim(string(chi2all(i_rejects[i],2)),2)
 close, unit
 free_lun, unit

 openw, unit, ofold+fbase+'_band_keepers.txt', /get_lun
 printf, unit, '# '+strtrim(string(cnt_keep),2)+' '+spt+' dwarfs used for template construction because'
 printf, unit, '# reduced chi2 < '+strtrim(string(sigma),2)+' in all bands (but not full spectrum)'
 printf, unit, '# Name'+tb+'J chi2'+tb+'H chi2'+tb+'K chi2'
 if (cnt_keep gt 0) then for i=0,cnt_keep-1 do printf, unit, str.names(i_keep[i])+' ' + tb+strtrim(string(chi2all(i_keep[i],0)),2)+tb+strtrim(string(chi2all(i_keep[i],1)),2)+tb+strtrim(string(chi2all[i_keep[i]],2),2)
 close, unit
 free_lun, unit

 rstr = spt + '!C'+'!9s!3 = '+strtrim(string(sigma),2)+' N_kept = '+strtrim(string(cnt_keep),2)+'!C'+ 'N_reject = '+strtrim(string(cnt_rejects),2)+'!C'+'N_spec = '+strtrim(string(cnt_keep+cnt_rejects),2)

 if i_rejects[0] eq -1 then print, 'no band-by-band rejects'; else print, "band-by-band rejects:", str.names(wnselect)

 print, " number selected = ", cnt_keep
 print, " number rejected = ", cnt_rejects

;=================== 
;plot the entire range
;=================== 
w = where(str.lam ge 0.9 and str.lam le 2.4,cnt)
template_full = kellel_template( str.lam(w), str.flx(w,*), str.flx_unc(w,*), norm_facts=norm_facts,flags_in=flags_in, flags_out=flags_out, sd=sd, sigma=sigma, chi2=chi2,mins_templ=mins_templ,maxs_templ=maxs_templ,no_renorm=no_re_norm)
;  flags_in = flags_out

; plot out results if at end of loop
;  if (i eq nloops-1) then begin
	  flags=flags_out
 ;  wnselect = where(flags ne 0,cntnselect)
  ; wselect = where(flags eq 0,cntselect)
   
   ;SET UP THE PLOTS
   if ~keyword_set(ps) then wset,1
   	!p.multi=[3,3,2]
   	;find the minimum and max to setup plot
	flux_mins = fltarr(nspec) & flux_maxs=fltarr(nspec)
	for i = 0,nspec-1 do flux_mins[i] = MIN( str.flx[w,i] / norm_facts[i] ) 
	for i = 0,nspec-1 do flux_maxs[i] = MAX( str.flx[w,i] / norm_facts[i] )
	
   yra = [MIN(flux_mins),MAX(flux_maxs) ] 
   plot, str.lam(w), template_full, /xsty, yra=yra,/ysty, xtitle='!3Wavelength (!9m!3m)', ytitle='Normalized Flux', charsize=2, xmargin=[8,-70]
  
   ;if wnselect[0] eq -1 then print, 'no full band rejects' else print, "full band rejects:", str.names(wnselect)
   ;print, " number selected = ",cntselect
   ;print, " number rejected = ",cntnselect

   case 1 of 
    ptype eq 0: begin
    	if (cnt_rejects gt 0) then begin
	 		for j=0,cnt_rejects-1 do begin 
		 		oplot, str.lam[w], str.flx(w,i_rejects[j])/norm_facts[i_rejects[j]], color=color_rejected[j], thick=1
				;xyouts, 0.95, yra(1)-0.15*j/2.5*(yra(1)-yra(0)), str.names(wnselect(j)),color=color_rejected[j],size=0.8
			endfor
		endif
	    polyfill, [str.lam[w],reverse(str.lam[w])], [maxs_templ,reverse(mins_templ)], color=colorfill_minmax, /fill
	    polyfill, [str.lam[w],reverse(str.lam[w])], [template_full+sd,reverse(template_full-sd)], color=colorfill_sdev, /fill
		
    end
    ptype eq 1: begin ;only plot keepers and the std dev
		if (cnt_keep gt 0) then for j=0,cnt_keep-1 do oplot, str.lam[w], str.flx(w,i_keep[j])/norm_facts[i_keep[j]], color=color_kept, thick=2
	    ;polyfill, [str.lam[w],reverse(str.lam[w])], [maxs_templ,reverse(mins_templ)], color=colorfill_minmax, /fill
	    polyfill, [str.lam[w],reverse(str.lam[w])], [template_full+sd,reverse(template_full-sd)], color=colorfill_sdev, /fill
		
    end
    else: begin ;only plot rejects
		if (cnt_rejects gt 0) then for j=0,cnt_rejects-1 do oplot, str.lam[w], str.flx(w,i_rejects[j])/norm_facts[i_rejects[j]], color=color_rejected[j], thick=1
	    polyfill, [str.lam[w],reverse(str.lam[w])], [maxs_templ,reverse(mins_templ)], color=colorfill_minmax, /fill
	    polyfill, [str.lam[w],reverse(str.lam[w])], [template_full+sd,reverse(template_full-sd)], color=colorfill_sdev, /fill		
    end
   endcase
   
 oplot, str.lam[w], template_full, thick=3
  
 xyouts, 2.35, yra(1)-0.15*(yra(1)-yra(0)), rstr , align=1, charsize=1.5

; save template
   dat = [[str.lam[w]],[template_full],[sd]]
   fxhmake,hdr,dat
   writefits, ofold+spt+'_template_FULL.fits', dat, hdr
  ;endif
;endfor

;print, ''

; openw, unit, sfold+fbase+'_full_rejects.txt', /get_lun
; wnselect = where(flag ne 0,cntnselect)
; if (cntnselect gt 0) then for i=0,cntnselect-1 do printf, unit, str.names(wnselect(i))+tb+strtrim(string(chi2all(wnselect(i),0)),2)+tb+strtrim(string(chi2all(wnselect(i),1)),2)+tb+strtrim(string(chi2all(wnselect(i),2)),2) else printf, unit, 'No rejects'
; close, unit
; free_lun, unit
 
!p.multi=0

if keyword_set(ps) then begin
	device, /close
	set_plot, 'x'
endif

FINISH: return
end