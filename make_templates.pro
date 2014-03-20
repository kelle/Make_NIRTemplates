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

;fact = total(flxn)
norm_fact = int_tabulated(lamn,flxn)

if ~finite(norm_fact) then print, 'problem with normalization'
;print, 'norm', fact

return, flx/norm_fact
end

; -----------------------------
; SCRIPT TO CREATE TEMPLATES
; -----------------------------
function kellel_template, lam, flxs, flx_uncs, norm_facts=norm_facts, flags_in=flags_in, flags_out=flags_out, flags_one=flags_one,sigma=sigma, sd=sd, chi2=chi2, maxs_templ=maxs_templ, mins_templ=mins_templ

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
  flxns(*,j) = kellel_normalize(lam, flxs(*,j), norm_fact=norm_fact)
  flxn_uncs(*,j) = flx_uncs(*,j) / norm_fact
  norm_facts(j) = norm_fact
  ;if j eq 0 then plot, lam, flxns(*,j) else oplot, lam, flxns(*,j)
endfor

;COMBINE THE SELECTED normalized spectra using the uncertainties
mc_meancomb2, flxns[*,wselect], template_first, datavar=flxn_uncs[*,wselect]

; repeat normalization against template
for j=0,nspec-1 do begin
	norm_facts(j) = median( flxs[*,j] / template_first )
	flxns[*,j] = flxs[*,j] / norm_facts[j]
	;print,min(flxns[*,j]),max(flxns[*,j])
endfor
  	
;recalculate template using re-normalized spectra
mc_meancomb2, flxns[*,wselect], template, datavar=flxn_uncs(*,wselect)

;CHARACTERIZE NEW TEMPLATE
mins_templ = MIN(flxns[*,wselect],dim=2)
maxs_templ = MAX(flxns[*,wselect],dim=2)

; standard deviation of spectra used in template at each wavelength
; currently unweighted
sd = fltarr(nlam)
for j = 0, nlam-1 do sd[j] = stddev(flxns[j,wselect])

; identify outliers of new template
; JUST FLAG (REJECT LATER) based on reduced chi square goodness of fit test
chi2 = fltarr(nspec)
for ispec=0,nspec-1 do begin
  chi2[ispec] = total( (flxns[*,ispec]-template)^2/sd^2) / (n_elements(template)-1.)
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
dfold = sfold+spt+'s/'
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
	sfiles = file_search(dfold+'*.fits')
	snr = fltarr(n_elements(sfiles))
	for ifile=0,n_elements(sfiles)-1 do begin
	 fits = readfits(sfiles(ifile),hd,/silent)
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
	 flx_unc(*,ifile) = interpol(fits(*,2),fits(*,0),lam)
	 ;snr(ifile) = max(smooth(flx(*,ifile)/ns(*,ifile),9))
	 ;plot, lam,flx(*,ifile)
	endfor

	str = create_struct($
	'file',sfiles,$
	'lam',lam,$
	'flx',flx,$
	'flx_unc',flx_unc,$
	;'snr',snr, $ 
	;'names',strrep(strrep(strrep(sfiles,dfold,''),'.fits',''),'spex_prism_',''))
	'names', strrep(strrep(sfiles,dfold,''),'.fits',''))	
 	 save, str, file=strfile

endif else restore, file=strfile

nspec=n_elements(sfiles)
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
nloops = 6

flags_in = intarr(nspec)
flags_j = intarr(nspec)
flags_h = intarr(nspec)
flags_k = intarr(nspec)
chi2_j = fltarr(nloops)
chi2_h = fltarr(nloops)
chi2_k = fltarr(nloops)
chi2all = fltarr(nspec,4)
i_spec = indgen(nspec)
i_loop = indgen(nloops)
 
;Make template 10 times to iterativly reject
for iloop=0,nloops-1 do begin
		
	;Loop over bands
	;Reject any spectrum that has chi^2 > 2 compared to template
	for mmm=0,n_elements(band)-1 do begin

		w = where(str.lam ge lam_norm(0,mmm) and str.lam le lam_norm(1,mmm),cnt)
		
		;takes flags_in and modifies to reflect new rejects
		template = kellel_template(str.lam(w), str.flx(w,*), str.flx_unc(w,*), norm_facts=norm_facts, flags_in=flags_in, flags_out=flags_out, flags_one=flags_one,sd=sd, sigma=sigma, chi2=chi2, mins_templ=mins_templ, maxs_templ=maxs_templ)
		
		i_keep = where(flags_out eq 0,cnt_keep)
		i_rejects = where(flags_out ge 1,cnt_reject)
				
		;keep track of band-by-band flags
		case mmm of 
			0: BEGIN	
				flags_j = flags_j + flags_one
				chi2_j[iloop] = mean(chi2[i_keep])			
			END
			1: BEGIN
				flags_h = flags_h + flags_one
				chi2_h[iloop] = mean(chi2[i_keep])	
			END
			2: BEGIN
				flags_k = flags_k + flags_one
				chi2_k[iloop] = mean(chi2[i_keep])	
			END
		endcase
			
		i_rejects_j = where(flags_j ge 1,cnt_reject_j)
		i_rejects_h = where(flags_h ge 1,cnt_reject_h)
		i_rejects_k = where(flags_k ge 1,cnt_reject_k)
	 
		if (cnt_keep eq 0) then message, 'Warning, rejected all sources!'
		
		;print, 'loop & mean chi^2 = ' + strn(iloop) + ' ' + strn(mean(chi2[i_keep]))
		;print, 'nkept, n_rejected = ', cnt_keep, cnt_reject
		;print, 'n_rejects , j ,h,k = ', cnt_reject_j, cnt_reject_h, cnt_reject_k
		;
		;if mmm eq 2 then print, ' '
			
		; PLOT	chi^sq for each band	
		if mmm eq 0 then begin
			window, 0
			!p.multi=[0,3,1]
			!p.thick=2
		endif else begin
			wset, 0
			!p.multi=[3-mmm,3,1]
		endelse			
		plot, i_spec, chi2, psym=3, title=strn(mmm)	, xr=[-1,nspec],xstyle=1,symsize=2,/ylog
		if i_rejects[0] ne -1 then oplot, i_spec[i_rejects], chi2[i_rejects], psym=6, symsize=3.5,color=red
		if i_rejects_j[0] ne -1 then oplot, i_spec[i_rejects_j], chi2[i_rejects_j], psym=5, symsize=3,color=cyan
		if i_rejects_h[0] ne -1 then oplot, i_spec[i_rejects_h], chi2[i_rejects_h], psym=4, symsize=2,color=green
		if i_rejects_k[0] ne -1 then oplot, i_spec[i_rejects_k], chi2[i_rejects_k], psym=1, symsize=2,color=magenta
		
		plots, [0,nspec-1], [mean(chi2[i_keep]), mean(chi2[i_keep])]
		plots, [0,nspec-1], [sigma,sigma], linestyle=1
		xyouts, mmm/3.0, 0, 'mean chi^2 =' + strn(mean(chi2[i_keep])), /normal
				
		flags_in = flags_out
		;flags = flag+f
		chi2all(*,mmm) = chi2
		 
		 
		 
	; plot out results if at end of loop
		if (iloop eq nloops-1) then begin
			flags=flags_out
		   	wnselect = where(flags ne 0,cntnselect)
			wnselect_j = where(flags_j ne 0,cntnselect_j)
			wnselect_h = where(flags_h ne 0,cntnselect_h)
			wnselect_k = where(flags_k ne 0,cntnselect_k)
		   	wselect = where(flags eq 0,cntselect) 
		   
		   ;SET UP THE PLOTS
		   if mmm eq 0 then begin
			   window,1
			   !p.multi=[0,3,2]
		   endif else begin
			   wset,1
			   !p.multi=[6-mmm,3,2]
		   endelse
			flux_mins = fltarr(nspec) & flux_maxs=fltarr(nspec)
	   		for i = 0,nspec-1 do flux_mins[i] = MIN( str.flx[w,i] / norm_facts[i] ) 
	   		for i = 0,nspec-1 do flux_maxs[i] = MAX( str.flx[w,i] / norm_facts[i] )
	      yra = [MIN(flux_mins),MAX(flux_maxs) ] 

		   ;plot template
		   plot, str.lam(w), template, /xsty, yra=yra,/ysty, xtitle='!3Wavelength (!9m!3m)', ytitle='Normalized Flux', charsize=2, xmargin=[8,2], title=tit
		   color_rejected=[green,cyan,purple,red,dkgreen,orange,gray,magenta,pink,ltorange,ltgray,purple2,dkred,dkorange,dkgrey,dkpurple, green,cyan,purple,red,dkyellow,dkgreen,orange,gray]
		   color_kept=blue
		   colorfill_sdev = gray
		   colorfill_minmax = vltgray
		    
		   case 1 of 
		    	ptype eq 0: begin
					;plot each object rejected from template
		     	   if (cntnselect gt 0) then for j=0,cntnselect-1 do oplot, str.lam(w), str.flx(w,wnselect(j))/norm_facts(wnselect(j)), color=color_rejected[j], thick=1
		     	   ;plot each target included in template
				   ;if (cntselect gt 0) then for j=0,cntselect-1 do oplot, str.lam(w), str.flx(w,wselect(j))/scl(wselect(j))/sclt, color=color_kept, thick=1
  		   		   ;polyfill the min-max of the spectra used in the template
  		   		   polyfill, [str.lam(w),reverse(str.lam(w))], [maxs_templ,REVERSE(mins_templ)], color=colorfill_minmax, /fill		   
				   ;polyfill the standard dev
  		   		   polyfill, [str.lam(w),reverse(str.lam(w))], [template+sd,reverse(template-sd)], color=colorfill_sdev, /fill
				   case mmm of
				   	0: if (cntnselect_j gt 0) then for j=0,cntnselect_j-1 do xyouts, 0.15+0.33*(mmm), 0.7-j*0.02, str.names(wnselect_j(j)),color=color_rejected[j],size=0.8,/normal
					1: if (cntnselect_h gt 0) then for j=0,cntnselect_h-1 do xyouts, 0.15+0.33*(mmm), 0.7-j*0.02, str.names(wnselect_h(j)),color=color_rejected[j],size=0.8,/normal
					2: if (cntnselect_k gt 0) then for j=0,cntnselect_k-1 do xyouts, 0.15+0.33*(mmm), 0.7-j*0.02, str.names(wnselect_k(j)),color=color_rejected[j],size=0.8,/normal
				ENDCASE
		    	end
		    	ptype eq 1: begin ;only plot objects selected
		     	   if (cntselect gt 0) then for j=0,cntselect-1 do oplot, str.lam(w), str.flx(w,wselect(j))/scl(wselect(j))/sclt, color=color_kept, thick=1
				   ;polyfill the standard dev of the template
  		   		   polyfill, [str.lam(w),reverse(str.lam(w))], [template+sd,reverse(template-sd)]/sclt, color=colorfill_sdev, /fill
		    	end
		    	else: begin ;only plot objects NOT selected	
		     	   if (cntnselect gt 0) then for j=0,cntnselect-1 do oplot, str.lam(w), str.flx(w,wnselect(j))/norm_facts(wnselect(j))/sclt, color=color_rejected[j], thick=1
  		   		   ;polyfill the min-max of the spectra used in the template
  		   		   polyfill, [str.lam(w),reverse(str.lam(w))], [max_templ,REVERSE(min_templ)]/sclt, color=colorfill_minmax, /fill		   
				   ;polyfill the standard dev
  		   		   polyfill, [str.lam(w),reverse(str.lam(w))], [template+sd,reverse(template-sd)]/sclt, color=colorfill_sdev, /fill
		    	end
		   endcase
		   
		   ;polyfill the min-max
		   polyfill, [str.lam(w),reverse(str.lam(w))], [maxs_templ,REVERSE(mins_templ)], color=colorfill_minmax, /fill
		   
		   ;polyfill the standard dev
		   polyfill, [str.lam(w),reverse(str.lam(w))], [template+sd,reverse(template-sd)], color=colorfill_sdev, /fill

		   oplot, str.lam(w), template, thick=2
		 		  
		   ;calculate new normalization for template
		   

		   ; save template
		   dat = [[str.lam(w)],[template],[sd],[mins_templ],[maxs_templ]]
		   fxhmake,hdr,dat
		   output_fits = ofold+spt+'_template_'+band(mmm)+'.fits'
		   writefits, output_fits , dat, hdr
		   message, 'wrote' + output_fits, /info 
  	 	endif
	
	 endfor ;end loop over 3 (JHK) bands, mmm

 	if keyword_set(inter) then begin & print,'press any key to continue' & tmp=GET_KBRD() & endif
	 ;stop ;after each loop
	 
endfor ;end iterative 10 loops, iloop

window,3
!p.multi=0
plot, i_loop, i_loop, yr =[0.5,1.1], xr=[-1,nloops+1],/nodata,xstyle=1
plots,[0,nloops],[1,1]
oplot, i_loop, chi2_j, psym = 5, color=cyan, symsize=2
oplot, i_loop, chi2_h, psym = 4, color=green, symsize=2
oplot, i_loop, chi2_k, psym = 1, color=magenta, symsize=2

;WRITE TXT FILES
 tb='	'
 openw, unit, ofold+fbase+'_band_rejects.txt', /get_lun
 ;wnselect = where(flags ne 0,cntnselect)
 printf, unit, '# '+strtrim(string(cntnselect),2)+' '+spt+' dwarfs rejected from template construction because '
 printf, unit, '# reduced chi2 > '+strtrim(string(sigma),2)+' in any band (but not full spectrum)'
 printf, unit, '# Name'+tb+'J chi2'+tb+'H chi2'+tb+'K chi2'
 if (cntnselect gt 0) then for i=0,cntnselect-1 do printf, unit, str.names(wnselect(i))+' '+tb+strtrim(string(chi2all(wnselect(i),0)),2)+tb+strtrim(string(chi2all(wnselect(i),1)),2)+tb+strtrim(string(chi2all(wnselect(i),2)),2)
 ;if (cntnselect gt 0) then for i=0,cntnselect-1 do printf, unit, str.file(wnselect(i))+tb+strtrim(string(chi2all(wnselect(i),0)),2)+tb+strtrim(string(chi2all(wnselect(i),1)),2)+tb+strtrim(string(chi2all(wnselect(i),2)),2)
 close, unit
 free_lun, unit

 openw, unit, ofold+fbase+'_band_keepers.txt', /get_lun
 ;wselect = where(flags eq 0,cntselect)
 printf, unit, '# '+strtrim(string(cntselect),2)+' '+spt+' dwarfs used for template construction because'
 printf, unit, '# reduced chi2 < '+strtrim(string(sigma),2)+' in all bands (but not full spectrum)'
 printf, unit, '# Name'+tb+'J chi2'+tb+'H chi2'+tb+'K chi2'
 if (cntselect gt 0) then for i=0,cntselect-1 do printf, unit, str.names(wselect(i))+' ' + tb+strtrim(string(chi2all(wselect(i),0)),2)+tb+strtrim(string(chi2all(wselect(i),1)),2)+tb+strtrim(string(chi2all(wselect(i),2)),2)
 ;if (cntselect gt 0) then for i=0,cntselect-1 do printf, unit, str.file(wselect(i))+tb+strtrim(string(chi2all(wselect(i),0)),2)+tb+strtrim(string(chi2all(wselect(i),1)),2)+tb+strtrim(string(chi2all(wselect(i),2)),2)
 close, unit
 free_lun, unit

 rstr = spt + '!C'+'!9s!3 = '+strtrim(string(sigma),2)+' N_kept = '+strtrim(string(cntselect),2)+'!C'+ 'N_reject = '+strtrim(string(cntnselect),2)+'!C'+'N_spec = '+strtrim(string(cntselect+cntnselect),2)

 if wnselect[0] eq -1 then print, 'no band-by-band rejects'; else print, "band-by-band rejects:", str.names(wnselect)

 print, " number selected = ", cntselect
 print, " number rejected = ", cntnselect

 ;=====================
; repeat for full band
;=====================
;flags_in(*) = 0

;for i=0,nloops-1 do begin

	; select "normal" sources
 ;	wselect = where(flags_in eq 0,cntselect)
	
 	;if (i ne nloop-1) then flags_in(*) = 0		; reset except for last loop
 ;	if (cntselect eq 0) then message, 'Warning, rejected all sources!'
 
w = where(str.lam ge 0.9 and str.lam le 2.4,cnt)
  
  ;if i eq nloops-1 then begin
;	  wset,0
;	  !p.multi=[1,4,1]
  template_full = kellel_template( str.lam(w), str.flx(w,*), str.flx_unc(w,*), norm_facts=norm_facts,flags_in=flags_in, flags_out=flags_out, sd=sd, sigma=sigma, chi2=chi2,mins_templ=mins_templ,maxs_templ=maxs_templ)
;  flags_in = flags_out

; plot out results if at end of loop
;  if (i eq nloops-1) then begin
	  flags=flags_out
 ;  wnselect = where(flags ne 0,cntnselect)
  ; wselect = where(flags eq 0,cntselect)
   
   ;SET UP THE PLOTS
   ;sclt = max(template) = 1
   	wset,1
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
     if (cntnselect gt 0) then begin
	 	for j=0,cntnselect-1 do begin 
		 	oplot, str.lam(w), str.flx(w,wnselect(j))/norm_facts(wnselect(j)), color=color_rejected[j], thick=1
			;xyouts, 0.95, yra(1)-0.15*j/2.5*(yra(1)-yra(0)), str.names(wnselect(j)),color=color_rejected[j],size=0.8
		endfor
	endif
     ;if (cntselect gt 0) then for j=0,cntselect-1 do oplot, str.lam(w), str.flx(w,wselect(j))/scl(wselect(j))/sclt, color=color_kept, thick=2
    end
    ptype eq 1: begin
     if (cntselect gt 0) then for j=0,cntselect-1 do oplot, str.lam(w), str.flx(w,wselect(j))/norm_facts(wselect(j)), color=color_kept, thick=2
    end
    else: begin
     if (cntnselect gt 0) then for j=0,cntnselect-1 do oplot, str.lam(w), str.flx(w,wnselect(j))/norm_facts(wnselect(j)), color=color_rejected[j], thick=1
    end
   endcase
   
 polyfill, [str.lam(w),reverse(str.lam(w))], [maxs_templ,reverse(mins_templ)], color=colorfill_minmax, /fill
 polyfill, [str.lam(w),reverse(str.lam(w))], [template_full+sd,reverse(template_full-sd)], color=colorfill_sdev, /fill


   oplot, str.lam(w), template_full, thick=3

   ;xyouts, 2.35, yra(1)-0.15*(yra(1)-yra(0)), '!9s!3 = '+strtrim(string(fix(sigma)),2), align=1, charsize=1.5
   xyouts, 2.35, yra(1)-0.15*(yra(1)-yra(0)), rstr , align=1, charsize=1.5

; save template
   dat = [[str.lam(w)],[template_full],[sd]]
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
