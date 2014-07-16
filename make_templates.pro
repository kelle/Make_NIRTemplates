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
pro make_templates, type=type, btypes=btypes, reset=reset, sigma=sigma, ptype=ptype, batch=batch, ps=ps, interactive=interactive, young=young, optical=optical

; if N_params() lt 1 then begin
; 	print, n_params()
;     print, "Syntax -  print, 'kellel, type=type, [/RESET, sigma=, ptype, batch, /PS]' "
; 	print, "type = 2 = L2. default sigma=2, type=1,ptype=0"
; print, "ptype = 0 (default, all), 1 (kept+template), 2 (rejects+template)"
; print, "Do more than one spectral typo at once with /batch then btypes=[1,4] or indgen(9)"
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
if (n_elements(type) eq 0 and ~keyword_set(batch)) then message, 'must specify type or btypes'
if (n_elements(ptype) eq 0) then ptype=0

; -----------------------------
; BATCH PROCESS
; -----------------------------

if (keyword_set(batch)) then begin
  ;types = indgen(3)+1
  for i=0,n_elements(btypes)-1 do begin
   ;for j=0,2 do begin
    ;for k=1,2 do make_templates, type=btypes[i], ptype=ptype, sigma=k*1. , ps=ps, reset=reset
	make_templates, type=btypes[i], ptype=ptype, sigma=sigma , ps=ps, reset=reset,young=young
   ;endfor
  endfor
  goto, FINISH
 endif
 
; -----------------------------
; EXAMINE VARIANCE IN SPECTRA, ATTEMPT TO DEFINE "NORMAL"
; -----------------------------

spt = 'L'+strmid(strtrim(string(type),2),0,1)
if keyword_set(young) then spt = spt+'lg'
print, spt, '  chi^2: ',strn(sigma)
print, 'p(lot) type: ', strn(ptype)
sfold = bfold+spt+'/'
spec_fold = '/Users/kelle/Dropbox/Data/nir_spectra_low/'
spec_fold_opt = '/Users/kelle/Dropbox/Data/optical_spectra/'
dfold = sfold+spt+'s/'
;slist = dfold+spt+'s_in.txt'

if keyword_set(young) then slist = bfold+'optNIR_LowG_Opt.txt' else slist = bfold+'spectra_in.txt'

txt_fold = bfold+'NIRSpecFigures/data/' ;sfold+'output_'+spt+'/'
ofold = bfold
temp_fold = bfold+'templates/'
strfile = bfold+'save_files/spectra'+spt+'.dat'
band = ['J','H','K']

if (file_search(strfile) eq '' or keyword_set(reset)) then begin
	MESSAGE,'re-selecting spectra',/info
	
	READCOL,slist,sfiles_all_opt,sfiles_all,stypes,format='A,A,I',DELIM = string(9b),comment='#'
	MESSAGE,'using spectra listed in ' + slist,/info
	
	sfiles=sfiles_all[where(stypes eq type+10,cnt_files)]
	sfiles_opt=sfiles_all_opt[where(stypes eq type+10)]
	Message, 'reading in '+ strn(cnt_files) +' files where type=' +strn(type+10),/info
	
	;====================
	; BOOTSTRAP GOES HERE
	;  
	; i_files_choosen = int(n_spec * randu(19))
 	; files_chosen = files[i_files_chosen]
	;
	;====================
	
	for ifile=0,cnt_files-1 do begin
	 fits = READFITS(spec_fold+sfiles[ifile],hd,/silent)
	 if sfiles_opt[ifile] eq 'include' then begin
		 fits_opt = 0.1*1.0+fltarr(n_elements(lam_opt)) 
		 message, strn(ifile) +' ' +sfiles_opt[ifile]+ ' --include file found',/info
		endif ELSE $
	 	fits_opt = KREADSPEC(spec_fold_opt+sfiles_opt[ifile],hd,/silent,/norm)
	 
	 
	 
	 if (ifile eq 0) then begin
		 ;setup the arrays
	  lam = fits(*,0)
	  flx = fltarr(n_elements(lam),cnt_files)
	  flx_unc = flx*0.
	  lam_opt = reform(fits_opt[0,*])
	  flx_opt = fltarr(n_elements(lam_opt),cnt_files)
	 endif
	 
	 ;interpolate onto wavelength of first object
	 flx(*,ifile) = interpol(fits(*,1),fits(*,0),lam)
	 flx_unc(*,ifile) = interpol(fits[*,2],fits[*,0],lam)
	 if sfiles_opt[ifile] ne 'include' then $ 
	 	flx_opt[*,ifile] = smooth(interpol(fits_opt[1,*],fits_opt[0,*],lam_opt),5) else $
		flx_opt[*,ifile] = fits_opt
	 ;print,ifile,flx_opt[0:4,ifile]
	endfor
	
	str = create_struct($
	'files',sfiles,$
	'files_opt',sfiles_opt,$
	'lam',lam,$
	'flx',flx,$
	'flx_unc',flx_unc,$
	'lam_opt',lam_opt,$
	'flx_opt',flx_opt,$
	'names', strrep(strrep(sfiles,dfold,''),'.fits',''))	
 	 
	 save, str, file=strfile
endif else begin
	restore, file=strfile
	MESSAGE,'restoring spectra from save file, not re-reading',/info
endelse

nspec=n_elements(str.files)
print,'nspec=',nspec

fbase = 'comp_'+spt+'_s'+strmid(strtrim(string(sigma),2),0,3)

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
make_plot=0
make_plot2=0

if keyword_set(ps) then set_plot,'ps'

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
		ENDCASE
			
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
		;=======
		;i_keep_j = 
		;i_keep_h = 
		;i_keep_k = 
		
		;if cnt_keep_j ge 1 then ii_keep_j = intarr(cnt_keep_j)
		;if cnt_keep_h ge 1 then ii_keep_h = intarr(cnt_keep_h)
		;if cnt_keep_k ge 1 then ii_keep_k = intarr(cnt_keep_k)
		;for i=0,cnt_keep_j-1 do ii_keep_j[i] = where(i_keep eq i_keep_j[i])
		;for i=0,cnt_keep_h-1 do ii_keep_h[i] = where(i_keep eq i_keep_h[i])
		;for i=0,cnt_keep_k-1 do ii_keep_k[i] = where(i_keep eq i_keep_k[i])	 
	 
	 	if mmm eq 2 then print,iloop,cnt_rejects_j,cnt_rejects_h,cnt_rejects_k
	 
		if (cnt_keep eq 0) then message, 'Warning, rejected all sources!'
		
		;print, 'loop & mean chi^2 = ' + strn(iloop) + ' ' + strn(mean(chi2[i_keep]))
		;print, 'nkept, n_rejected = ', cnt_keep, cnt_reject
		;print, 'n_rejects , j ,h,k = ', cnt_reject_j, cnt_reject_h, cnt_reject_k
		;
		;if mmm eq 2 then print, ' '
			
		; PLOT	chi^sq for each band	
		if ~keyword_Set(ps) then begin
			make_plot2 = 1
			if mmm eq 0 then begin
				window, 0
				!p.multi=[0,3,1]
				!p.thick=2
			endif else begin
				wset, 0
				!p.multi=[3-mmm,3,1]
			endelse
		endif 
		
		if keyword_set(ps) and mmm eq 0 then begin
			device, /encapsulated, ysize=24, xsize=18, filename=ofold+fbase+'_p'+strtrim(string(ptype),2)+'.eps', /portrait, bits_per_pixel=8, /color
			print, 'Creating: '+ ofold+fbase+'_p'+strtrim(string(ptype),2)+'.eps'
			;device, /encapsulated, ysize=18, xsize=24, filename=ofold+fbase+'_v'+strtrim(string(ptype+1),2)+'_chi2.eps', /portrait, bits_per_pixel=8, /color
			!p.multi=[0,3,3]
			!p.thick=2
			make_plot2=1
		endif
		if keyword_set(ps) and mmm ge 1 then begin
			!p.multi=[9-mmm,3,3]
			make_plot2=1
		ENDIF
		
		;plot the chi2 for each object in each band
		IF make_plot2 eq 1 then begin	
			plot, i_spec, chi2, psym=7, title=strn(mmm)	, xr=[-1,nspec],xstyle=1,symsize=0.5,/ylog,/ynozero
			if cnt_rejects ge 1 then oplot, i_spec[i_rejects], chi2[i_rejects], psym=6, symsize=2.5,color=red
			
			if cnt_rejects_j ge 1 then for j=0,cnt_rejects_j-1 do oplot, [i_spec[i_rejects_j[j]]], [chi2[i_rejects_j[j]]], psym=5, symsize=2,color=color_rejected[ii_rejects_j[j]]
			if cnt_rejects_h ge 1 then for j=0,cnt_rejects_h-1 do oplot, [i_spec[i_rejects_h[j]]], [chi2[i_rejects_h[j]]], psym=4, symsize=2,color=color_rejected[ii_rejects_h[j]]
			if cnt_rejects_k ge 1 then for j=0,cnt_rejects_k-1 do oplot, [i_spec[i_rejects_k[j]]], [chi2[i_rejects_k[j]]], psym=1, symsize=2,color=color_rejected[ii_rejects_k[j]]
			
			;oplot, i_spec[i_rejects_j], chi2[i_rejects_j], psym=5, symsize=3,color=cyan
			;if i_rejects_h[0] ne -1 then oplot, i_spec[i_rejects_h], chi2[i_rejects_h], psym=4, symsize=2,color=green
			;if cnt_rejects_k[0] ge 1 then oplot, i_spec[i_rejects_k], chi2[i_rejects_k], psym=1, symsize=2,color=magenta
		
			plots, [0,nspec-1], [mean(chi2[i_keep]), mean(chi2[i_keep])]
			plots, [0,nspec-1], [sigma,sigma], linestyle=1
			xyouts, mmm/3.0+0.1, 0.72, 'mean chi^2 =' + strn(mean(chi2[i_keep]),length=5), /normal
		ENDIF
				
		flags_in = flags_out
		;flags = flag+f
		chi2all(*,mmm) = chi2
		 
		if (iloop ge 1) and ~keyword_set(ps) then  make_plot =1 
		;if keyword_set(ps) and mmm eq 0 then begin	
		;	device, /encapsulated, ysize=18, xsize=24, filename=ofold+fbase+'_v'+strtrim(string(ptype+1),2)+'.eps', /portrait, bits_per_pixel=8, /color
		;	print, 'Creating: '+ ofold+fbase+'_v'+strtrim(string(ptype+1),2)+'.eps'
		;endif 
		if keyword_set(ps) then make_plot =1
		

		if make_plot eq 1 then begin
			
			;SET UP THE PLOTS
			if ~keyword_set(ps) then begin
			   case mmm of
			    0: begin
				   window,1
				   !p.multi=[0,3,2]
				end
				else: begin
					wset,1
					!p.multi=[6-mmm,3,2]
				end
			   ENDCASE
			endif 
			if keyword_set(ps) then begin
 			   case mmm of
 			    0: begin
 				   !p.multi=[6,3,3]
 				end
 				else: begin
 					!p.multi=[6-mmm,3,3]
 				end
			   endcase			
			endif
			
			flux_mins = fltarr(nspec) & flux_maxs=fltarr(nspec)
	   		for i = 0,nspec-1 do flux_mins[i] = MIN( str.flx[w,i] / norm_facts[i] ) 
	   		for i = 0,nspec-1 do flux_maxs[i] = MAX( str.flx[w,i] / norm_facts[i] )
	    	yra = [MIN(flux_mins),MAX(flux_maxs) ] 

		   ;plot template
		   if mmm eq 0 then plot, str.lam(w), template, /xsty, yra=yra,/ysty, ytitle='Normalized Flux', charsize=2, xmargin=[8,-2],/nodata
		   if mmm eq 1 then plot, str.lam(w), template, /xsty, yra=yra,/ysty, charsize=2, xmargin=[6,0],/nodata
		   if mmm eq 2 then plot, str.lam(w), template, /xsty, yra=yra,/ysty, charsize=2, xmargin=[4,2],/nodata
		   
		   ;set up label locations
		   name_loc_x=0.14
		   name_loc_offset_x=0.3
		   name_loc_y=0.45
		   
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
					   	0: if (cnt_rejects_j gt 0) then for j=0,cnt_rejects_j-1 do xyouts, name_loc_x+name_loc_offset_x*(mmm), name_loc_y-j*0.02, str.names[i_rejects_j[j]],color=color_rejected[ii_rejects_j[j]],size=0.8,/normal
						1: if (cnt_rejects_h gt 0) then for j=0,cnt_rejects_h-1 do xyouts, name_loc_x+name_loc_offset_x*(mmm), name_loc_y-j*0.02, str.names[i_rejects_h[j]],color=color_rejected[ii_rejects_h[j]],size=0.8,/normal
						2: if (cnt_rejects_k gt 0) then for j=0,cnt_rejects_k-1 do xyouts, name_loc_x+name_loc_offset_x*(mmm), name_loc_y-j*0.02, str.names[i_rejects_k[j]],color=color_rejected[ii_rejects_k[j]],size=0.8,/normal
					ENDCASE
		    	end
		    	ptype eq 1: begin ;only plot objects selected & the std dev
		     	   if (cnt_keep gt 0) then for j=0,cnt_keep-1 do oplot, str.lam[w], str.flx(w,i_keep[j])/norm_facts[i_keep[j]], color=color_rejected[j], thick=1
				   ;polyfill the standard dev of the template
  		   		   ;polyfill, [str.lam[w],reverse(str.lam[w])], [template+sd,reverse(template-sd)], color=colorfill_sdev, /fill
				   
				case mmm of
				   	0: if (cnt_keep gt 0) then for j=0,cnt_keep-1 do xyouts, name_loc_x+0.33*(mmm), name_loc_y-j*0.02, str.names[i_keep[j]],color=color_rejected[j],size=0.8,/normal
					1: if (cnt_keep gt 0) then for j=0,cnt_keep-1 do xyouts, name_loc_x+0.33*(mmm), name_loc_y-j*0.02, str.names[i_keep[j]],color=color_rejected[j],size=0.8,/normal
					2: if (cnt_keep gt 0) then for j=0,cnt_keep-1 do xyouts, name_loc_x+0.33*(mmm), name_loc_y-j*0.02, str.names[i_keep[j]],color=color_rejected[j],size=0.8,/normal
				ENDCASE
				   
		    	end
		    	else: begin ;plot everything
  		   		   if (cnt_keep gt 0) then for j=0,cnt_keep-1 do oplot, str.lam[w], str.flx(w,i_keep[j])/norm_facts[i_keep[j]], color=color_kept, thick=1
				   if (cnt_rejects gt 0) then for j=0,cnt_rejects-1 do oplot, str.lam(w), str.flx(w,i_rejects[j])/norm_facts(i_rejects[j]), color=color_rejected[j], thick=1
					case mmm of
					   	0: if (cnt_rejects_j gt 0) then for j=0,cnt_rejects_j-1 do xyouts, name_loc_x+0.33*(mmm), name_loc_y-j*0.02, str.names[i_rejects_j[j]],color=color_rejected[ii_rejects_j[j]],size=0.8,/normal
						1: if (cnt_rejects_h gt 0) then for j=0,cnt_rejects_h-1 do xyouts, name_loc_x+0.33*(mmm), name_loc_y-j*0.02, str.names[i_rejects_h[j]],color=color_rejected[ii_rejects_h[j]],size=0.8,/normal
						2: if (cnt_rejects_k gt 0) then for j=0,cnt_rejects_k-1 do xyouts, name_loc_x+0.33*(mmm), name_loc_y-j*0.02, str.names[i_rejects_k[j]],color=color_rejected[ii_rejects_k[j]],size=0.8,/normal
					ENDCASE
		    	end
		   endcase
		   oplot, str.lam(w), template, thick=2
   		endif ;end make_plot
		
	 	;if loop_stop eq 1 then begin
		   ; save template
		   dat = [[str.lam(w)],[template],[sd],[mins_templ],[maxs_templ]]
		   fxhmake,hdr,dat
		   output_fits = temp_fold+spt+'_template_'+band(mmm)+'.fits'
		   writefits, output_fits , dat, hdr
		   message, 'wrote' + output_fits, /info 
 	  	;endif

	endfor ;end loop over 3 (JHK) bands, mmm

 	if keyword_set(interactive) and iloop ge 1 then begin & print,'press any key to continue' & tmp=GET_KBRD() & endif
		iloop++
endwhile	
;endfor ;end iterative 10 loops, iloop

nloops=iloop
i_loop = indgen(nloops)



;WRITE TXT FILES
 tb='	'
 openw, unit, txt_fold+fbase+'_band_rejects.txt', /get_lun
 printf, unit, '# '+strtrim(string(cnt_rejects),2)+' '+spt+' dwarfs rejected from template construction because '
 printf, unit, '# reduced chi2 > '+strtrim(string(sigma),2)+' in any band (but not full spectrum)'
 printf, unit, '# Name'+tb+'J chi2'+tb+'H chi2'+tb+'K chi2'
 if (cnt_rejects gt 0) then for i=0,cnt_rejects-1 do printf, unit, str.names(i_rejects[i])+' '+tb+strtrim(string(chi2all(i_rejects[i],0)),2)+tb+strtrim(string(chi2all(i_rejects[i],1)),2)+tb+strtrim(string(chi2all(i_rejects[i],2)),2)
 close, unit
 free_lun, unit

 openw, unit, txt_fold+fbase+'_band_keepers.txt', /get_lun
 printf, unit, '# '+strtrim(string(cnt_keep),2)+' '+spt+' dwarfs used for template construction because'
 printf, unit, '# reduced chi2 < '+strtrim(string(sigma),2)+' in all bands (but not full spectrum)'
 printf, unit, '# Name'+tb+'J chi2'+tb+'H chi2'+tb+'K chi2'
 if (cnt_keep gt 0) then for i=0,cnt_keep-1 do printf, unit, str.names(i_keep[i])+' ' + tb+strtrim(string(chi2all(i_keep[i],0)),2)+tb+strtrim(string(chi2all(i_keep[i],1)),2)+tb+strtrim(string(chi2all[i_keep[i]],2),2)
 close, unit
 free_lun, unit

 rstr = spt + '!C'+'!9s!3 = '+strtrim(string(sigma),2)+' N_kept = '+strtrim(string(cnt_keep),2)+'!C'+ 'N_reject = '+strtrim(string(cnt_rejects),2)+'!C'+'N_spec = '+strtrim(string(cnt_keep+cnt_rejects),2)

 if i_rejects[0] eq -1 then print, 'no band-by-band rejects'; else print, "band-by-band rejects:", str.names(wnselect)
 
 message, 'wrote:' + txt_fold+fbase+'_band_rejects.txt',/info
 message, 'wrote:' + txt_fold+fbase+'_band_keepers.txt',/info

 print, " number selected = ", cnt_keep
 print, " number rejected = ", cnt_rejects

;=================== 
;plot the entire range
;=================== 
w = where(str.lam ge 0.9 and str.lam le 2.4,cnt)
template_full = kellel_template( str.lam(w), str.flx(w,*), str.flx_unc(w,*), norm_facts=norm_facts,flags_in=flags_in, flags_out=flags_out, sd=sd, sigma=sigma, chi2=chi2,mins_templ=mins_templ,maxs_templ=maxs_templ,no_renorm=no_re_norm)
   
   ;SET UP THE PLOTS
   if ~keyword_set(ps) then wset,1
   	!p.multi=[3,3,3]
   	;find the minimum and max to setup plot
	flux_mins = fltarr(nspec) & flux_maxs=fltarr(nspec)
	for i = 0,nspec-1 do flux_mins[i] = MIN( str.flx[w,i] / norm_facts[i] ) 
	for i = 0,nspec-1 do flux_maxs[i] = MAX( str.flx[w,i] / norm_facts[i] )
	
   yra = [MIN(flux_mins),MAX(flux_maxs) ] 
   plot, str.lam(w), template_full, /xsty, yra=yra,/ysty, xtitle='!3Wavelength (!9m!3m)', ytitle='Normalized Flux', charsize=2, xmargin=[8,-50]
 
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
		if (cnt_keep gt 0) then for j=0,cnt_keep-1 do oplot, str.lam[w], str.flx(w,i_keep[j])/norm_facts[i_keep[j]], color=color_rejected[j], thick=1
		;if (cnt_keep gt 0) then for j=0,cnt_keep-1 do oplot, str.lam[w], str.flx(w,i_keep[j])/norm_facts[i_keep[j]], color=color_kept, thick=2
	    ;polyfill, [str.lam[w],reverse(str.lam[w])], [maxs_templ,reverse(mins_templ)], color=colorfill_minmax, /fill
	    ;polyfill, [str.lam[w],reverse(str.lam[w])], [template_full+sd,reverse(template_full-sd)], color=colorfill_sdev, /fill
		
    end
    else: begin ;only plot rejects
		if (cnt_rejects gt 0) then for j=0,cnt_rejects-1 do oplot, str.lam[w], str.flx(w,i_rejects[j])/norm_facts[i_rejects[j]], color=color_rejected[j], thick=1
	    polyfill, [str.lam[w],reverse(str.lam[w])], [maxs_templ,reverse(mins_templ)], color=colorfill_minmax, /fill
	    polyfill, [str.lam[w],reverse(str.lam[w])], [template_full+sd,reverse(template_full-sd)], color=colorfill_sdev, /fill		
    end
   endcase
   
 oplot, str.lam[w], template_full, thick=3
  
 xyouts, 2.35, yra(1)-0.15*(yra(1)-yra(0)), rstr , align=1, charsize=1.3

; save template
   dat = [[str.lam[w]],[template_full],[sd]]
   fxhmake,hdr,dat
   writefits, temp_fold+spt+'_template_FULL.fits', dat, hdr
  ;endif
;endfor

;===============
;plot the optical spectra
;===============
if ~keyword_set(ps) then window, 4
!p.multi=0

;stop

;normalize optical
; make norm_facts_opt

;flux_mins_opt = fltarr(nspec) & flux_maxs_opt=fltarr(nspec)
;for i = 0,nspec-1 do flux_mins_opt[i] = MIN( str.flx_opt[w,i] / norm_facts[i] ) 
;for i = 0,nspec-1 do flux_maxs_opt[i] = MAX( str.flx_opt[w,i] / norm_facts[i] )

 ;yra = [MIN(flux_mins),MAX(flux_maxs) ] 
  plot, str.lam_opt, str.flx_opt[*,0], xr=[6500,9000],/xsty, /ysty, xtitle='!3Wavelength ()', ytitle='Normalized Flux', charsize=2,/nodata; xmargin=[8,-50]

  ;stop
  
  name_loc_x=0.2
  name_loc_y=0.8

  case 1 of 
   ptype eq 0: begin
   	if (cnt_rejects gt 0) then begin
 		for j=0,cnt_rejects-1 do begin 
	 		oplot, str.lam_opt, str.flx_opt(*,i_rejects[j]), color=color_rejected[j], thick=1
			xyouts, name_loc_x,name_loc_y+0.02,'REJECTS:',/normal
			xyouts, name_loc_x, name_loc_y-j*0.02, str.files_opt[i_rejects[j]],color=color_rejected[j],size=0.8,/normal
		endfor
	endif	
   end
   ptype eq 1: begin ;only plot keepers and the std dev
	if (cnt_keep gt 0) then begin
		for j=0,cnt_keep-1 do begin			
			oplot, str.lam_opt, str.flx_opt(*,i_keep[j]), color=color_rejected[j], thick=1	
			xyouts, name_loc_x, name_loc_y+0.02, 'KEEPERS:',/normal
			xyouts, name_loc_x, name_loc_y-j*0.02, str.files_opt[i_keep[j]],color=color_rejected[j],size=0.8,/normal
;			stop
		endfor
	endif
   end
   else: begin ;only plot rejects
	if (cnt_rejects gt 0) then begin
		for j=0,cnt_rejects-1 do begin
			oplot, str.lam_opt, str.flx_opt(*,i_rejects[j]), color=color_rejected[j], thick=1
			xyouts, name_loc_x,name_loc_y+0.02,'REJECTS:',/normal
			xyouts, name_loc_x, name_loc_y-j*0.02, str.files_opt[i_rejects[j]],color=color_rejected[j],size=0.8,/normal
		endfor
	endif
   end
  endcase
   
xyouts, 2.35, yra(1)-0.15*(yra(1)-yra(0)), rstr , align=1, charsize=1.3

;===============
;plot the chi2s
;================
if ~keyword_set(ps) then begin 
	window,3
endif else if keyword_set(ps) then begin
	device, /encapsulated, ysize=18, xsize=24, filename=ofold+fbase+'_chi2avg.eps', /portrait, bits_per_pixel=8, /color
ENDIF
!p.multi=0
plot, i_loop, i_loop, yr =[0.5,1.1], xr=[-1,nloops+1],/nodata,xstyle=1
plots,[0,nloops],[1,1]
oplot, i_loop, chi2_j[0:nloops-1], psym = 5, color=cyan, symsize=2
oplot, i_loop, chi2_h[0:nloops-1], psym = 4, color=green, symsize=2
oplot, i_loop, chi2_k[0:nloops-1], psym = 1, color=magenta, symsize=2

!p.multi=0

if keyword_set(ps) then begin
	device, /close
	set_plot, 'x'
endif

FINISH: return
end