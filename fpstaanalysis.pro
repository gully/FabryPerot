pro fpSTAanalysis

;Author: gully
;date: May 8, 2014
;desc:  Sample Transport Accessory analysis

pi=3.14159265D0

;Outline / pseudocode

;0) Import the data.
;1) Wrangle, take means, stddev's based on positions.
;2) Plot it.

;0) Import the data, manipulate the arrays [this should be a sub-function]
fn1='20140508_STA_Si_mgsTrim.csv'
dir1='/Volumes/cambridge/Astronomy/silicon/ACT/bonding/cary5000/20140508/'
cd, dir1
;at=ascii_template(fn1) ;you must use "GROUP ALL in the third ascii-template window!!"
;save, at, /verbose, filename=file_basename(fn1, '.csv')+'.sav'
restore, file_basename(fn1, '.csv')+'.sav', /verbose
;Need a special reader for the Sample Transport Accessory data because it is so big.
; #BigData
d=cary_5000_csv_STA(fn1, at)

wl=float(d.wl)
nwls=n_elements(wl)

;---------------------
;1) Wrangle, take means, stddev's based on positions.
;---------------------
;Reform the array:
thisArr=d.dat[1:*, *]
thisArr=thisArr/100.0
outArr=reform(thisArr, 21, 13, 51) ; Position, sample, wavelength bin.  Wow!

meanWL=median(outArr, dimension=3)
sigWL=SIG_ARRAY(outArr, 3)

xpos=findgen(21)*5.0+0.0

;---------------------
;2) Plot it.
;---------------------

c1=fpObstructionFig(xpos, meanWL, sigWL)


end

function fpObstructionFig, xpos, meanWL, sigWL
device, decomposed=1
outdir='/Users/gully/IDLWorkspace/FabryPerot/figs/'
outname='obstructionSTA.eps'
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

;xtit=greek('lambda', /append_font)+' (nm)'
xtit='x (mm)'
device, /helvetica
;ytit='T!De!N'
ytit='T (%)'

!p.multi=0

Polyfill, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=cgColor('Papaya');'Pale Goldenrod')

loadct, 13

id=0
thisY=meanWL[*, id]
thisErr=sigWL[*, id]
plot, xpos, thisY, xtitle=xtit, ytitle=ytit, thick=3.0, charthick=2.0,$
 yrange=[0.999, 1.001], psym=10, charsize=0.8, $
 xrange=[-10.0, 110.0], xstyle=1, ystyle=1, /nodata, /noerase

;oplot, xpos, thisY, psym=10, color=0, thick=3.0
oploterror, xpos, thisY, thisErr,thisErr, thick=3.0, errthick=3.0, errcolor=0

id=1
thisY=meanWL[*, id]
thisErr=sigWL[*, id]
;oplot, xpos, thisY, psym=10, color=255, thick=3.0
oploterror, xpos, thisY, thisErr,thisErr, thick=3.0, psym=10, color=255, errthick=3.0, errcolor=255
 
badx1=[4.0, 16.0]
badx2=[79.0, 91.0]
yband1=[100.0, 100.0]
yband2=[-100.0, -100.0]

oband,badx1,yband1,yband2,color=220
oband,badx2,yband1,yband2,color=220 
 

leg_text=['No sample holder', $
          'With sample holder']
legend,leg_text,linestyle=[10, 10], color=[0,255], $
  position=[20, 1.0005], charthick=2.0

XYOUTS, 7.0, 0.9995, 'blocked by sample holder', ORIENTATION=75.0
XYOUTS, 81.0, 0.9995, 'blocked by sample holder', ORIENTATION=75.0

Device, /Close_File
Set_Plot, thisDevice
Obj_Destroy, psObject


spawn, 'open /Users/gully/IDLWorkspace/FabryPerot/figs/obstructionSTA.eps'
return, 1

end

