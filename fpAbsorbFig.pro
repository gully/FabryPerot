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
ratio1 = reform(d.dat[14, *]/100.0)


;---------------------
;2) Plots
;---------------------
device, decomposed=1
outdir='/Users/gully/IDLWorkspace/FabryPerot/figs/'
outname='fpAbsorbfig.eps'
psObject = Obj_New("FSC_PSConfig",/Color, /Helvetica, /Bold, $
            Filename=outdir+outname, xsize=5.0, ysize=4.5, /encapsulate)
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

Polyfill, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=cgColor('Papaya');'Pale Goldenrod')
!p.multi=0

loadct, 13
plot, wl, reform(d.dat[14, *]/100.0), xtitle=xtit, ytitle=ytit, thick=3.0, charthick=2.0,$
 yrange=[0.95, 1.05], psym=10, $
 xrange=[1000.0, 2500.0], xstyle=1, ystyle=1, /nodata, /noerase
 
oplot, wl, ratio1, color=255, psym=10, thick=3.0
ones=wl*0.0+1.0
oplot, wl, ones, linestyle=2, thick=3.0
;oplot, wl, reform(d.dat[13, *]/100.0), color=100.0, linestyle=1, thick=3.0

leg_text=['Measured', $
          'No drift']
legend,leg_text,linestyle=[10, 2], color=[255,0], $
  position=[1300.0, 0.8], charthick=2.0

Device, /Close_File
Set_Plot, thisDevice
Obj_Destroy, psObject

end
