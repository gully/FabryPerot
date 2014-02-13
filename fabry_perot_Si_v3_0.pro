pro fabry_perot_Si_v3_0

;Author: gully
;date: Feb 21, 2013
;desc: This program makes a Fabry Perot model for the Si voids
; in between bonded Si pucks
;The original prototpe for this program was an Excel spreadsheet:
; Si_optics_Feb2013

pi=3.14159265D0

;Outline / pseudocode

;0) Import the data, manipulate the arrays [this should be a sub-function]
;1) Define wavelength range
;2) Set the range of temperature
;3) Determine the refractive index
;4) Calculate the Fabry-Perot Physics [This should be a sub-function]
;[4.1) Set the Si absorption per unit length]
;4.2) Set the gap widths
;5) Plot the data and overlay the models 


;0) Import the data, manipulate the arrays [this should be a sub-function]
fn1='20130607_cary5000_ascii.csv'
dir1='/Users/gully/Astronomy/Silicon/ACT/transmission_tests/20130607/'
cd, dir1
at=ascii_template(fn1) ;you must use "GROUP ALL in the third ascii-template window!!"
save, at, /verbose, filename=file_basename(fn1, '.csv')+'.sav'
restore, file_basename(fn1, '.csv')+'.sav', /verbose
d=cary_5000_csv(fn1, at)

;1) Wavelength
wl=float(d.wl)
nwls=n_elements(wl)

;---------------------
;2) Physics here:
T_imr = incoh_mult_reflect_v2_0(wl, 50.0)
;---------------------
print, 1

;---------------------
;3) Plots
;---------------------
device, decomposed=1
outname='20130607_dat01-02.eps'
psObject = Obj_New("FSC_PSConfig",/Color, /Helvetica, /Bold, $
            Filename=outname, xsize=5.0, ysize=4.5, /encapsulate)
thisDevice = !D.Name
Set_Plot, "PS"
!p.font=0
Device, _Extra=psObject->GetKeywords()

nc=10
c1=ceil(findgen(nc)/nc*255)

xtit=greek('lambda', /append_font)+' (nm)'
device, /helvetica
;ytit='T!De!N'
ytit='Transmission'

Polyfill, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=cgColor('Papaya');'Pale Goldenrod')
!p.multi=0

loadct, 13
plot, wl, reform(d.dat[2, *]/100), xtitle=xtit, ytitle=ytit, thick=3.0, charthick=2.0,$
 yrange=[0.0, 1.00], psym=10, $
 xrange=[1000.0, 2500.0], xstyle=1, ystyle=1, /nodata, /noerase
 
oplot, wl, reform(d.dat[1, *])/100.0, color=255, psym=10, thick=3.0
oplot, wl, T_imr.t_net, color=0.0, linestyle=1, thick=3.0
oplot, wl, T_imr.t_old, color=50.0, linestyle=3, thick=3.0
 
Device, /Close_File
Set_Plot, thisDevice
Obj_Destroy, psObject



print, 1


end


function incoh_mult_reflect_v2_0, wl, dgap ;input wavelength must be in nm
;author: gully
;date: 6/7/2013
;desc: figures out the multiple reflections from the front and backs of Si pucks
; this only works for incoherent multiple reflections
;input the wavelength array
;outputs the net Fresnel transmission
 

;1) Set the temperature (this doesn't matter that much for room temp)
temp1=295.0 ;K

;2) Determine the refractive index
nwls=n_elements(wl)
n1=fltarr(nwls)

for i=0, nwls-1 do begin
  n1[i]=sellmeier_si(wl[i], temp1)
endfor

;3) Silicon reflectance (Fresnel losses at 1 interface)
R0=((n1-1.0)/(n1+1.0))^2.0

;4) Coefficient of Finesse
F=4.0*R0/(1.0-R0)^2.0

delta=2.0*3.141592654*dgap/wl

T_net=2.0*n1/(1.0+2.0*n1*F*sin(delta)^2.0+n1^2.0)
T_old=1.0/(1.0+F*sin(delta)^2.0)

ret_struct={t_net:t_net, t_old:t_old}
return, ret_struct

end