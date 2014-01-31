; -----------------------------
; SEPARATE NORMALIZATION SCRIPT
; -----------------------------
function kellel_normalize, lam, flx, fact=fact, rng=rng

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
fact = int_tabulated(lamn,flxn)

if ~finite(fact) then print, 'problem with normalization'
;print, 'norm', fact

return, flx/fact
end

; -----------------------------
; SCRIPT TO CREATE TEMPLATES
; -----------------------------
function kellel_template, lam, flx, scl=scl, flag=flag, sigma=sigma, sd=sd, chi2=chi2

nloop = 10
nspec = n_elements(flx(0,*))
nlam = n_elements(flx(*,0))
scl = fltarr(nspec)
flxn = flx
if (n_elements(sigma) eq 0) then sigma = 3.
if (n_elements(flag) eq 0) then flag = fltarr(nspec)*0

wselect = where(flag eq 0,cntselect)
if (cntselect eq 0) then message, 'No spectra available for creating template'

; first normalization 
 for j=0,nspec-1 do begin
  flxn(*,j) = kellel_normalize(lam,flx(*,j),fact=fact)
  scl(j) = fact
  ;print,j,scl(j)
 endfor
 template = total(flxn(*,wselect),2)/(1.*cntselect)

; repeat normalization against template
 for kkk=0,4 do begin
  for j=0,nspec-1 do begin
   scl(j) = median(flx(*,j)/template)
   flxn(*,j) = flx(*,j)/scl(j)
  endfor
  template = total(flxn(*,wselect),2)/(1.*cntselect)
 endfor

; standard deviation
 sd = template*0.
 for j=0,nlam-1 do sd(j) = stddev(flxn(j,wselect))
 diff = flx*0.
 chi2 = fltarr(nspec)
; plot, lam, template, /xsty

; identify outlers
 for j=0,nspec-1 do begin

; reject based on max difference
;  diff = abs(flxn(*,j)-template)/sd
;  if (max(diff) gt sigma) then flag(j) = flag(j)+1

; reject based on chi square
  chi2(j) = total((flxn(*,j)-template)^2/sd^2)/(n_elements(template)-1.)
  if (chi2(j) gt sigma) then begin
	  flag(j) = flag(j)+1
  endif

;  oplot, lam, flxn(*,j), color=150
 endfor
; oplot, str.lam(w), template
; oplot, str.lam(w), template+sd, linestyle=1
; oplot, str.lam(w), template-sd, linestyle=1

 flag = flag<1


;if ~keyword_set(ps) then begin
;  	plot,lam,template
;  	print,"press any key to continue"
;  	gibbirish=get_kbrd(1)
;endif

return, template
end


; -----------------------------
; MAIN ROUTINE
; -----------------------------
pro kellel, type=type, btypes=btypes, reset=reset, sigma=sigma, ptype=ptype, batch=batch, ps=ps

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


tb='	'
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
print, spt, '  sigma: ',strn(sigma)
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
for i=0,n_elements(sfiles)-1 do begin
 fits = readfits(sfiles(i),hd,/silent)
 if (i eq 0) then begin
  lam = fits(*,0)
  flx = fltarr(n_elements(lam),n_elements(sfiles))
  ns = flx*0.
  wnorm = where(lam ge 1.2 and lam le 1.3)		; initial normalization
 endif
 flx(*,i) = interpol(fits(*,1),fits(*,0),lam)
 mx = max(flx(wnorm,i))
 flx(*,i) = flx(*,i)/mx
 ns(*,i) = interpol(fits(*,2),fits(*,0),lam)/mx
 snr(i) = max(smooth(flx(*,i)/ns(*,i),9))
 plot, lam,flx(*,i)
endfor

 str = create_struct($
 'file',sfiles,$
 'lam',lam,$
 'flx',flx,$
 'ns',ns,$
 'snr',snr, $
 'names',strrep(strrep(strrep(sfiles,sfold+'spectra/',''),'.fits',''),'spex_prism_',''))
 
save, str, file=strfile
endif else restore, file=strfile

; iterative scan, normalizing in individual bands
flag = intarr(n_elements(str.file))*0
lam_norm = [[0.87,1.39],[1.41,1.82],[1.91,2.39]]
nloop = 10

fbase = 'comp_'+spt+'_s'+strmid(strtrim(string(sigma),2),0,3)

if keyword_set(ps) then begin	
	set_plot, 'ps'
	device, /encapsulated, ysize=18, xsize=24, filename=ofold+fbase+'_v'+strtrim(string(ptype+1),2)+'.eps', /portrait, bits_per_pixel=8, /color
	print, 'Creating: '+ ofold+fbase+'_v'+strtrim(string(ptype+1),2)+'.eps'
endif 

!p.multi=[0,3,2]

for i=0,nloop-1 do begin

; select "normal" sources
 wselect = where(flag eq 0,cntselect)
 if (i ne nloop-1) then flag(*) = 0		; reset except for last loop
 if (cntselect eq 0) then message, 'Warning, rejected all sources!'
 chi2all = fltarr(n_elements(flag),3)
 
 ;Loop over bands
 for mmm=0,n_elements(band)-1 do begin
	print,'band  '  ,strn(mmm)
	w = where(str.lam ge lam_norm(0,mmm) and str.lam le lam_norm(1,mmm),cnt)
	template = kellel_template(str.lam(w),str.flx(w,*),scl=scl,flag=f, sd=sd, sigma=sigma, chi2=chi2)
	flag = flag+f
	chi2all(*,mmm) = chi2

; plot out results if at end of loop
  if (i eq nloop-1) then begin
   wnselect = where(flag ne 0,cntnselect)
   wselect = where(flag eq 0,cntselect)
   sclt = max(template)
   yra = [min(template-3.*sd),max(template+3.*sd)]/sclt

   print, "number selected = ",cntselect
   print, "number rejected = ",cntnselect
   
   ;plot template
   plot, str.lam(w), template/sclt, /xsty, yra=yra,/ysty, xtitle='!3Wavelength (!9m!3m)', ytitle='Normalized Flux', charsize=2, xmargin=[8,2], title=tit
   
   ; color_rejected=vltgray
   ; color_kept=blue
   ; colorfill_template=vltblue
   
   ; color_rejected=vltgray
   color_rejected=[green,cyan,purple,red,dkyellow,dkgreen,orange,gray,magenta,pink,ltorange,ltgray,purple2,dkred,dkorange,dkgrey,dkpurple, green,cyan,purple,red,dkyellow,dkgreen,orange,gray]
   color_kept=blue
   colorfill_template=vltgray

   case 1 of 
    ptype eq 0: begin
     if (cntnselect gt 0) then for j=0,cntnselect-1 do oplot, str.lam(w), str.flx(w,wnselect(j))/scl(wnselect(j))/sclt, color=color_rejected[j], thick=1
     if (cntselect gt 0) then for j=0,cntselect-1 do oplot, str.lam(w), str.flx(w,wselect(j))/scl(wselect(j))/sclt, color=color_kept, thick=1
    end
    ptype eq 1: begin
     if (cntselect gt 0) then for j=0,cntselect-1 do oplot, str.lam(w), str.flx(w,wselect(j))/scl(wselect(j))/sclt, color=color_kept, thick=1
    end
    else: begin
     if (cntnselect gt 0) then for j=0,cntnselect-1 do oplot, str.lam(w), str.flx(w,wnselect(j))/scl(wnselect(j))/sclt, color=color_rejected[j], thick=1
    end
   endcase
   
   polyfill, [str.lam(w),reverse(str.lam(w))], [template+sd,reverse(template-sd)]/sclt, color=colorfill_template, /fill
   
   oplot, str.lam(w), template/sclt, thick=2

; save template
   dat = [[str.lam(w)],[template/sclt],[sd/sclt]]
   fxhmake,hdr,dat
   writefits, ofold+'template_'+band(mmm)+'.fits', dat, hdr
  endif

 endfor
endfor

 openw, unit, ofold+fbase+'_band_rejects.txt', /get_lun
 wnselect = where(flag ne 0,cntnselect)
 printf, unit, '# '+strtrim(string(cntnselect),2)+' '+spt+' dwarfs rejected from template construction because '
 printf, unit, '# reduced chi2 > '+strtrim(string(sigma),2)+' in any band (but not full spectrum)'
 printf, unit, '# Name'+tb+'J chi2'+tb+'H chi2'+tb+'K chi2'
 if (cntnselect gt 0) then for i=0,cntnselect-1 do printf, unit, str.names(wnselect(i))+tb+strtrim(string(chi2all(wnselect(i),0)),2)+tb+strtrim(string(chi2all(wnselect(i),1)),2)+tb+strtrim(string(chi2all(wnselect(i),2)),2)
 ;if (cntnselect gt 0) then for i=0,cntnselect-1 do printf, unit, str.file(wnselect(i))+tb+strtrim(string(chi2all(wnselect(i),0)),2)+tb+strtrim(string(chi2all(wnselect(i),1)),2)+tb+strtrim(string(chi2all(wnselect(i),2)),2)
 close, unit
 free_lun, unit

 openw, unit, ofold+fbase+'_band_keepers.txt', /get_lun
 wselect = where(flag eq 0,cntselect)
 printf, unit, '# '+strtrim(string(cntselect),2)+' '+spt+' dwarfs used for template construction because'
 printf, unit, '# reduced chi2 < '+strtrim(string(sigma),2)+' in all bands (but not full spectrum)'
 printf, unit, '# Name'+tb+'J chi2'+tb+'H chi2'+tb+'K chi2'
 if (cntselect gt 0) then for i=0,cntselect-1 do printf, unit, str.names(wselect(i))+tb+strtrim(string(chi2all(wselect(i),0)),2)+tb+strtrim(string(chi2all(wselect(i),1)),2)+tb+strtrim(string(chi2all(wselect(i),2)),2)
 ;if (cntselect gt 0) then for i=0,cntselect-1 do printf, unit, str.file(wselect(i))+tb+strtrim(string(chi2all(wselect(i),0)),2)+tb+strtrim(string(chi2all(wselect(i),1)),2)+tb+strtrim(string(chi2all(wselect(i),2)),2)
 close, unit
 free_lun, unit

 rstr = '!9s!3 = '+strtrim(string(sigma),2)+' N = '+strtrim(string(cntselect),2)+'/'+strtrim(string(cntselect+cntnselect),2)

 if wnselect[0] eq -1 then print, 'no band-by-band rejects' else print, "band-by-band rejects:", str.names(wnselect)


; repeat for full band
flag(*) = 0
for i=0,nloop-1 do begin

; select "normal" sources
 wselect = where(flag eq 0,cntselect)
 if (i ne nloop-1) then flag(*) = 0		; reset except for last loop
 if (cntselect eq 0) then message, 'Warning, rejected all sources!'
 
  w = where(str.lam ge 0.9 and str.lam le 2.4,cnt)
  template = kellel_template(str.lam(w),str.flx(w,*),scl=scl,flag=f, sd=sd, sigma=sigma, chi2=chi2)
  flag = f

; plot out results if at end of loop
  if (i eq nloop-1) then begin
   wnselect = where(flag ne 0,cntnselect)
   wselect = where(flag eq 0,cntselect)
   sclt = max(template)
   yra = [min(template-3.*sd),max(template+3.*sd)]/sclt
   plot, str.lam(w), template/sclt, /xsty, yra=yra,/ysty, xtitle='!3Wavelength (!9m!3m)', ytitle='Normalized Flux', charsize=2, xmargin=[8,-70]
  
   if wnselect[0] eq -1 then print, 'no full band rejects' else print, "full band rejects:", str.names(wnselect)

   case 1 of 
    ptype eq 0: begin
     if (cntnselect gt 0) then for j=0,cntnselect-1 do oplot, str.lam(w), str.flx(w,wnselect(j))/scl(wnselect(j))/sclt, color=color_rejected[j], thick=1
     if (cntselect gt 0) then for j=0,cntselect-1 do oplot, str.lam(w), str.flx(w,wselect(j))/scl(wselect(j))/sclt, color=color_kept, thick=2
    end
    ptype eq 1: begin
     if (cntselect gt 0) then for j=0,cntselect-1 do oplot, str.lam(w), str.flx(w,wselect(j))/scl(wselect(j))/sclt, color=color_kept, thick=2
    end
    else: begin
     if (cntnselect gt 0) then for j=0,cntnselect-1 do oplot, str.lam(w), str.flx(w,wnselect(j))/scl(wnselect(j))/sclt, color=color_rejected[j], thick=1
    end
   endcase

 polyfill, [str.lam(w),reverse(str.lam(w))], [template+sd,reverse(template-sd)]/sclt, color=colorfill_template, /fill

   oplot, str.lam(w), template/sclt, thick=3

;   xyouts, 2.35, yra(1)-0.15*(yra(1)-yra(0)), '!9s!3 = '+strtrim(string(fix(sigma)),2), align=1, charsize=1.5
;   xyouts, 2.35, yra(1)-0.25*(yra(1)-yra(0)), strtrim(string(cntselect),2)+'/'+strtrim(string(cntselect+cntnselect),2)+' sources in template', align=1, charsize=1.5
   xyouts, 2.35, yra(1)-0.15*(yra(1)-yra(0)), rstr, align=1, charsize=1.5

; save template
   dat = [[str.lam(w)],[template/sclt],[sd/sclt]]
   fxhmake,hdr,dat
   writefits, ofold+'template_FULL.fits', dat, hdr
  endif
endfor

print, ''

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
