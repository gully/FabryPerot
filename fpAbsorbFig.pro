pro fpAbsorbFig


;Author: gully
;date: Feb 21, 2013
;Edited: March 6, 2014
;        May 1, 2014 
;desc:  Quantify absorption in Si by plotting the 
;       Si transmission through a puck divided by transmission 
;       through a thinner Si wafer.

pi=3.14159265D0

;Outline / pseudocode

;0) Import the data.
;1) Compute the ratio of transmissions.
;2) Plot it.

;0) Import the data, manipulate the arrays [this should be a sub-function]
fn1='20130218_cary5000_Si_all.csv'
dir1='/Volumes/cambridge/Astronomy/silicon/ACT/bonding/cary5000/20130218'
cd, dir1
;at=ascii_template(fn1) ;you must use "GROUP ALL in the third ascii-template window!!"
;save, at, /verbose, filename=file_basename(fn1, '.csv')+'.sav'
restore, file_basename(fn1, '.csv')+'.sav', /verbose
d=cary_5000_csv(fn1, at)

wl=float(d.wl)
nwls=n_elements(wl)

;---------------------
;1) Ratio
;---------------------
;Make a correction for the baseline, which drifted by 2% over 1 hour.
baseline2=d.dat[14, *]
baseline1=d.dat[0, *]
driftR=d.dat[16, *]/d.dat[2, *]
driftCorr = reform(baseline2/baseline1)
;Compute the Tpuck / Twafer
Tpuck=reform(d.dat[13, *])
TpuckAlt=reform(d.dat[4, *])
Twafer=reform(d.dat[16, *])
TwaferAlt=reform(d.dat[2, *])
Tp_Tw1 = Twafer / (Tpuck / driftCorr) 
Tp_Tw2 = TwaferAlt / (TpuckAlt)

;Compute the absolute deviation
Tp_Tw2Abs=abs(Tp_Tw2-1.0)

;---------------------
;2) Plots
;---------------------
device, decomposed=1
outdir='/Users/gully/IDLWorkspace/FabryPerot/figs/'
outname='fpAbsorbfig.eps'
psObject = Obj_New("FSC_PSConfig",/Color, /Helvetica, /Bold, $
            Filename=outdir+outname, xsize=7.5, ysize=4.0, /encapsulate)
thisDevice = !D.Name
Set_Plot, "PS"
!p.font=0
!p.thick=3.0
!x.thick=3.0
!y.thick=3.0
Device, _Extra=psObject->GetKeywords()

nc=10
c1=ceil(findgen(nc)/nc*255)

xtit=greek('lambda', /append_font)+' (nm)'
device, /helvetica
;ytit='T!De!N'
ytit='Transmission'

!p.multi=[0, 2, 1, 0, 1]

Polyfill, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=cgColor('Papaya');'Pale Goldenrod')

loadct, 13
plot, wl, Tp_Tw2, xtitle=xtit, ytitle=ytit, thick=3.0, charthick=2.0,$
 yrange=[0.96, 1.04], psym=10, charsize=0.8, $
 xrange=[1000.0, 2500.0], xstyle=1, ystyle=1, /nodata, /noerase
 
oplot, wl, Tp_Tw2, color=255, psym=10, thick=3.0
ones=wl*0.0+1.0
oplot, wl, ones, linestyle=2, thick=3.0
;oplot, wl, reform(d.dat[13, *]/100.0), color=100.0, linestyle=1, thick=3.0
oplot, [1200, 1200], [1.0E-8, 1.0E6], color=0, linestyle=2, thick=3.0

leg_text=['Measured', $
          'No drift']
legend,leg_text,linestyle=[10, 2], color=[255,0], $
  position=[1300.0, 0.8], charthick=2.0

!p.multi=[1, 2, 1, 0, 1]

loadct, 13
plot, wl, Tp_Tw2Abs, xtitle=xtit, ytitle='Absolute deviation', thick=3.0, charthick=2.0,$
 yrange=[1.0E-5, 100.0], psym=10, /ylog, charsize=0.8, $
 xrange=[1000.0, 2500.0], xstyle=1, ystyle=1, /nodata, /noerase
 
oplot, wl, Tp_Tw2Abs, color=255, psym=10, thick=3.0
oplot, [1200, 1200], [1.0E-8, 1.0E6], color=0, linestyle=2, thick=3.0
oplot, [800, 3000], [0.002, 0.002], color=100, linestyle=1, thick=3.0



Device, /Close_File
Set_Plot, thisDevice
Obj_Destroy, psObject

print, 1
end
