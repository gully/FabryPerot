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

;Calculate the predicted Fresnel spectrum (with Mult-Incoh-Reflects):
struct1=incoh_mult_reflect_v2_0(wl, 0.0)
t_net=struct1.t_net
meanFresWL=median(t_net)
sigFresWL=stddev(t_net)

;Compute the ratio of target 3 and target 5.
rgi=[5,6,7,8,9, 10, 11, 12] ;these have VG05 and VG06 not differing by more than 0.001
rat1=reform(outArr[rgi, 2, *]/outArr[rgi,4, *])

;---------------------
;2) Plot it.
;---------------------

;c1=fpObstructionFig(xpos, meanWL, sigWL)
goodi= where(xpos lt 5.0 or (xpos gt 15.0 and xpos lt 80.0) or xpos gt 95.0 )
;c2=fpDriftFig(outArr, goodi)
;c3=fpFresCheckFig(xpos, meanWL, sigWL)
c4 = fpAbs1250Fig( wl, rat1)
c5 = fpGapCheckFig(xpos, meanWL, sigWL)
c6 = fpGapCheckAltFig(xpos, meanWL, sigWL)
c7 = fpCherryPicFig(xpos, meanWL, sigWL, outArr)

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
ytit='Transmission'

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


function fpDriftFig, outArr, goodi
device, decomposed=1
outdir='/Users/gully/IDLWorkspace/FabryPerot/figs/'
outname='driftSTA.eps'
psObject = Obj_New("FSC_PSConfig",/Color, /Helvetica, /Bold, $
            Filename=outdir+outname, xsize=4.0, ysize=4.0, /encapsulate)
thisDevice = !D.Name
Set_Plot, "PS"
!p.font=0
!p.thick=3.0
!x.thick=3.0
!y.thick=3.0
Device, _Extra=psObject->GetKeywords()

nc=10
c1=ceil(findgen(nc)/nc*255)

xtit='Transmission'
device, /helvetica
ytit='N'

!p.multi=0

Polyfill, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=cgColor('Papaya');'Pale Goldenrod')

loadct, 13
bs=0.00002

thisY=outArr[goodi, 1, *]
h1 = HISTOGRAM( thisY, BINSIZE=bs, LOCATIONS=L1, MAX=1.001, MIN=0.999)
thisY=outArr[goodi, 3, *]
h2 = HISTOGRAM( thisY, BINSIZE=bs, LOCATIONS=L2, MAX=1.001, MIN=0.999)
thisY=outArr[goodi, 5, *]
h3 = HISTOGRAM( thisY, BINSIZE=bs, LOCATIONS=L3, MAX=1.001, MIN=0.999)
thisY=outArr[goodi, 7, *]
h4 = HISTOGRAM( thisY, BINSIZE=bs, LOCATIONS=L4, MAX=1.001, MIN=0.999)
thisY=outArr[goodi, 12, *]
h5 = HISTOGRAM( thisY, BINSIZE=bs, LOCATIONS=L5, MAX=1.001, MIN=0.999)

maxN=max(h5)

plot, L1, h1, xtitle=xtit, ytitle=ytit, thick=3.0, charthick=2.0,$
 xrange=[0.999, 1.001], psym=10, charsize=0.8, $
 yrange=[0, maxN*1.5], xstyle=1, ystyle=1, /nodata, /noerase

oplot, L1, H1, psym=10, color=0, thick=3.0
oplot, L1, H2, psym=10, color=50, thick=3.0
oplot, L1, H3, psym=10, color=100, thick=3.0
oplot, L1, H4, psym=10, color=150, thick=3.0
oplot, L1, H5, psym=10, color=255, thick=3.0

oplot, [1.0, 1.0], [-10, 1000], linestyle=1, thick=3.0

Device, /Close_File
Set_Plot, thisDevice
Obj_Destroy, psObject


spawn, 'open /Users/gully/IDLWorkspace/FabryPerot/figs/driftSTA.eps'
return, 1

end


function fpFresCheckFig, xpos, meanWL, sigWL
device, decomposed=1
outdir='/Users/gully/IDLWorkspace/FabryPerot/figs/'
outname='FresCheckSTA.eps'
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
ytit='Transmission'

!p.multi=0

Polyfill, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=cgColor('Papaya');'Pale Goldenrod')

loadct, 13

id=2
thisY=meanWL[*, id]
thisErr=sigWL[*, id]
plot, xpos, thisY, xtitle=xtit, ytitle=ytit, thick=3.0, charthick=2.0,$
 yrange=[0.500, 0.530], psym=10, charsize=0.8, $
 xrange=[-10.0, 110.0], xstyle=1, ystyle=1, /nodata, /noerase

;oplot, xpos, thisY, psym=10, color=0, thick=3.0
oploterror, xpos, thisY, thisErr,thisErr, psym=5, thick=3.0, errthick=3.0, errcolor=0

id=4
thisY=meanWL[*, id]
thisErr=sigWL[*, id]
;oplot, xpos, thisY, psym=10, color=255, thick=3.0
;oploterror, xpos, thisY, thisErr,thisErr, thick=3.0, psym=6, color=255, errthick=3.0, errcolor=255
 
badx1=[4.0, 16.0]
badx2=[79.0, 91.0]
badx3=[16, 22.0]
badx4=[64, 79.0]
yband1=[100.0, 100.0]
yband2=[-100.0, -100.0]

;Re-arrange so that the bands overplot some of the points, but not others.
oband,badx3,yband1,yband2,color=190 
oband,badx4,yband1,yband2,color=190  
oploterror, xpos, thisY, thisErr,thisErr, thick=3.0, psym=6, color=255, errthick=3.0, errcolor=255
oband,badx1,yband1,yband2,color=220
oband,badx2,yband1,yband2,color=220

leg_text=['VG05', $
          'VG06']
legend,leg_text,psym=[5, 6], color=[0,255], $
  position=[24, 0.508], charthick=2.0
  
loadct,0
mT=0.5275
sigT=0.0005482


xband1=[mT-sigT, mT-sigT]  
xband2=[mT+sigT, mT+sigT]
oband, [-1000, 1000], xband1, xband2, color=200
oplot, [-1000, 1000], [mT, mT], linestyle=2, color=0

XYOUTS, 7.0, 0.5070, 'blocked by sample holder', ORIENTATION=75.0
XYOUTS, 81.0, 0.5070, 'blocked by sample holder', ORIENTATION=75.0
XYOUTS, 18.0, 0.514, 'off of VG05 edge', ORIENTATION=80.0
XYOUTS, 69.0, 0.514, 'off of VG05 edge', ORIENTATION=80.0
XYOUTS, 30.0, mt+sigT*1.3, 'Fresnel Prediction', ORIENTATION=0.0

Device, /Close_File
Set_Plot, thisDevice
Obj_Destroy, psObject

spawn, 'open /Users/gully/IDLWorkspace/FabryPerot/figs/FresCheckSTA.eps'
return, 1

end


function fpAbs1250Fig, wl, rat1

device, decomposed=1
outdir='/Users/gully/IDLWorkspace/FabryPerot/figs/'
outname='fpAbs1250fig.eps'
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
ytit='Transmission ratio'

!p.multi=[0, 2, 1, 0, 1]

Polyfill, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=cgColor('Papaya');'Pale Goldenrod')

loadct, 13

plot, wl, rat1[0, *], xtitle=xtit, ytitle=ytit, thick=3.0, charthick=2.0,$
 yrange=[0.998, 1.002], psym=10, charsize=0.8, $
 xrange=[1250.0, 1350.0], xstyle=1, ystyle=1, /nodata, /noerase
 
oplot, wl, rat1[0, *], color=0, psym=10, thick=3.0
oplot, wl, rat1[1, *], color=40, psym=10, thick=3.0
oplot, wl, rat1[2, *], color=80, psym=10, thick=3.0
oplot, wl, rat1[3, *], color=120, psym=10, thick=3.0
oplot, wl, rat1[4, *], color=160, psym=10, thick=3.0
oplot, wl, rat1[5, *], color=200, psym=10, thick=3.0
oplot, wl, rat1[6, *], color=220, psym=10, thick=3.0
oplot, wl, rat1[7, *], color=255, psym=10, thick=3.0
ones=wl*0.0+1.0
oplot, wl, ones, linestyle=2, thick=3.0

leg_text=['Measured', $
          'No drift']
legend,leg_text,linestyle=[10, 2], color=[255,0], $
  position=[1300.0, 0.8], charthick=2.0

!p.multi=[1, 2, 1, 0, 1]

loadct, 13
plot, wl, abs(rat1[0, *]-1.0), xtitle=xtit, ytitle='Absolute deviation', thick=3.0, charthick=2.0,$
 yrange=[1.0E-5, 1.0E-2], psym=10, /ylog, charsize=0.8, $
 xrange=[1250, 1350], xstyle=1, ystyle=1, /nodata, /noerase
 
oplot, wl, abs(rat1[0, *]-1.0), color=0, psym=10, thick=3.0
oplot, wl, abs(rat1[1, *]-1.0), color=40, psym=10, thick=3.0
oplot, wl, abs(rat1[2, *]-1.0), color=80, psym=10, thick=3.0
oplot, wl, abs(rat1[3, *]-1.0), color=120, psym=10, thick=3.0
oplot, wl, abs(rat1[4, *]-1.0), color=160, psym=10, thick=3.0
oplot, wl, abs(rat1[5, *]-1.0), color=200, psym=10, thick=3.0
oplot, wl, abs(rat1[6, *]-1.0), color=220, psym=10, thick=3.0
oplot, wl, abs(rat1[7, *]-1.0), color=255, psym=10, thick=3.0
oplot, [800, 3000], [0.002, 0.002], color=100, linestyle=1, thick=3.0

Device, /Close_File
Set_Plot, thisDevice
Obj_Destroy, psObject
!P.multi=0

spawn, 'open /Users/gully/IDLWorkspace/FabryPerot/figs/fpAbs1250fig.eps'
return, 1

end


function fpGapCheckFig, xpos, meanWL, sigWL
device, decomposed=1
outdir='/Users/gully/IDLWorkspace/FabryPerot/figs/'
outname='GapCheckSTA.eps'
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
ytit='Transmission'

!p.multi=0

Polyfill, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=cgColor('Papaya');'Pale Goldenrod')

loadct, 13

id=6
thisY=meanWL[*, id]
thisErr=sigWL[*, id]
plot, xpos, thisY, xtitle=xtit, ytitle=ytit, thick=3.0, charthick=2.0,$
 yrange=[0.0, 0.60], psym=10, charsize=0.8, $
 xrange=[-10.0, 110.0], xstyle=1, ystyle=1, /nodata, /noerase

;oplot, xpos, thisY, psym=10, color=0, thick=3.0
oploterror, xpos, thisY, thisErr,thisErr, psym=5, thick=3.0, color=150, errthick=3.0, errcolor=150

id=8
thisY=meanWL[*, id]
thisErr=sigWL[*, id]
;oplot, xpos, thisY, psym=10, color=255, thick=3.0
oploterror, xpos, thisY, thisErr,thisErr, thick=3.0, psym=6, color=255, errthick=3.0, errcolor=255

id=9
thisY=meanWL[*, id]
thisErr=sigWL[*, id]
;oplot, xpos, thisY, psym=10, color=255, thick=3.0
oploterror, xpos, thisY, thisErr,thisErr, thick=3.0, psym=4, color=100, errthick=3.0, errcolor=100

id=10
thisY=meanWL[*, id]
thisErr=sigWL[*, id]
;oplot, xpos, thisY, psym=10, color=255, thick=3.0
oploterror, xpos, thisY, thisErr,thisErr, thick=3.0, psym=2, color=50, errthick=3.0, errcolor=50

id=11
rgi=[4, 5, 6, 7, 8, 9, 10, 11]
thisY=meanWL[rgi, id]
thisErr=sigWL[rgi, id]
;oplot, xpos, thisY, psym=10, color=255, thick=3.0
oploterror, xpos[rgi], thisY, thisErr,thisErr, thick=3.0, psym=10, color=0, errthick=3.0, errcolor=0

 
badx1=[4.0, 16.0]
badx2=[79.0, 91.0]
yband1=[100.0, 100.0]
yband2=[-100.0, -100.0]

;Re-arrange so that the bands overplot some of the points, but not others.
oband,badx1,yband1,yband2,color=220
oband,badx2,yband1,yband2,color=220

leg_text=['VG02', $
          'VG04', $
          'MIS01', $
          'VG08', $
          'VG05']
legend,leg_text,psym=[5, 6, 4, 2, -3], color=[150,255,100,50,0], $
  position=[24, 0.20], charthick=2.0
  
loadct,0
mT=0.5275
sigT=0.0005482


xband1=[mT-sigT, mT-sigT]  
xband2=[mT+sigT, mT+sigT]
oband, [-1000, 1000], xband1, xband2, color=200
oplot, [-1000, 1000], [mT, mT], linestyle=2, color=0

XYOUTS, 7.0, 0.10, 'blocked by sample holder', ORIENTATION=75.0
XYOUTS, 81.0, 0.10, 'blocked by sample holder', ORIENTATION=75.0
XYOUTS, 30.0, mt+sigT*1.3, 'Fresnel Prediction', ORIENTATION=0.0

Device, /Close_File
Set_Plot, thisDevice
Obj_Destroy, psObject

spawn, 'open /Users/gully/IDLWorkspace/FabryPerot/figs/GapCheckSTA.eps'
return, 1

end



function fpGapCheckAltFig, xpos, meanWL, sigWL
device, decomposed=1
outdir='/Users/gully/IDLWorkspace/FabryPerot/figs/'
outname='GapCheckAltSTA.eps'
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
ytit='Transmission'

!p.multi=0

Polyfill, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=cgColor('Papaya');'Pale Goldenrod')

loadct, 13

id=6
thisY=meanWL[*, id]
thisErr=sigWL[*, id]
plot, xpos, thisY, xtitle=xtit, ytitle=ytit, thick=3.0, charthick=2.0,$
 yrange=[0.495, 0.53], psym=10, charsize=0.8, $
 xrange=[-10.0, 110.0], xstyle=1, ystyle=1, /nodata, /noerase

;oplot, xpos, thisY, psym=10, color=0, thick=3.0
oploterror, xpos, thisY, thisErr,thisErr, psym=5, thick=3.0, color=150, errthick=3.0, errcolor=150

id=8
thisY=meanWL[*, id]
thisErr=sigWL[*, id]
;oplot, xpos, thisY, psym=10, color=255, thick=3.0
oploterror, xpos, thisY, thisErr,thisErr, thick=3.0, psym=6, color=255, errthick=3.0, errcolor=255

id=9
thisY=meanWL[*, id]
thisErr=sigWL[*, id]
;oplot, xpos, thisY, psym=10, color=255, thick=3.0
oploterror, xpos, thisY, thisErr,thisErr, thick=3.0, psym=4, color=100, errthick=3.0, errcolor=100

id=10
thisY=meanWL[*, id]
thisErr=sigWL[*, id]
;oplot, xpos, thisY, psym=10, color=255, thick=3.0
oploterror, xpos, thisY, thisErr,thisErr, thick=3.0, psym=2, color=50, errthick=3.0, errcolor=50

id=11
rgi=[4, 5, 6, 7, 8, 9, 10, 11]
thisY=meanWL[rgi, id]
thisErr=sigWL[rgi, id]
;oplot, xpos, thisY, psym=10, color=255, thick=3.0
oploterror, xpos[rgi], thisY, thisErr,thisErr, thick=3.0, psym=10, color=0, errthick=3.0, errcolor=0

 
badx1=[4.0, 16.0]
badx2=[79.0, 91.0]
yband1=[100.0, 100.0]
yband2=[-100.0, -100.0]

;Re-arrange so that the bands overplot some of the points, but not others.
oband,badx1,yband1,yband2,color=220
oband,badx2,yband1,yband2,color=220

leg_text=['VG02', $
          'VG04', $
          'MIS01', $
          'VG08', $
          'VG05']
legend,leg_text,psym=[5, 6, 4, 2, -3], color=[150,255,100,50,0], $
  position=[24, 0.507], charthick=2.0
  
loadct,0
mT=0.5275
sigT=0.0005482


xband1=[mT-sigT, mT-sigT]  
xband2=[mT+sigT, mT+sigT]
oband, [-1000, 1000], xband1, xband2, color=200
oplot, [-1000, 1000], [mT, mT], linestyle=2, color=0

XYOUTS, 7.0, 0.505, 'blocked by sample holder', ORIENTATION=75.0
XYOUTS, 81.0, 0.505, 'blocked by sample holder', ORIENTATION=75.0
XYOUTS, 30.0, mt+sigT*1.3, 'Fresnel Prediction', ORIENTATION=0.0

Device, /Close_File
Set_Plot, thisDevice
Obj_Destroy, psObject

spawn, 'open /Users/gully/IDLWorkspace/FabryPerot/figs/GapCheckAltSTA.eps'
return, 1

end



function fpCherryPicFig, xpos, meanWL, sigWL, outArr
device, decomposed=1
outdir='/Users/gully/IDLWorkspace/FabryPerot/figs/'
outname='CherryPicSTA.eps'
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
ytit='Transmission'

!p.multi=[0, 2, 1, 0, 1]

Polyfill, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=cgColor('Papaya');'Pale Goldenrod')

loadct, 13

id=6
thisY=meanWL[*, id]
thisErr=sigWL[*, id]
plot, xpos, thisY, xtitle=xtit, ytitle=ytit, thick=3.0, charthick=2.0,$
 yrange=[0.495, 0.53], psym=10, charsize=0.8, $
 xrange=[-10.0, 110.0], xstyle=1, ystyle=1, /nodata, /noerase

;oplot, xpos, thisY, psym=10, color=0, thick=3.0
oploterror, xpos, thisY, thisErr,thisErr, psym=5, thick=3.0, color=150, errthick=3.0, errcolor=150

id=8
thisY=meanWL[*, id]
thisErr=sigWL[*, id]
;oplot, xpos, thisY, psym=10, color=255, thick=3.0
oploterror, xpos, thisY, thisErr,thisErr, thick=3.0, psym=6, color=255, errthick=3.0, errcolor=255

id=9
thisY=meanWL[*, id]
thisErr=sigWL[*, id]
;oplot, xpos, thisY, psym=10, color=255, thick=3.0
oploterror, xpos, thisY, thisErr,thisErr, thick=3.0, psym=4, color=100, errthick=3.0, errcolor=100

id=10
thisY=meanWL[*, id]
thisErr=sigWL[*, id]
;oplot, xpos, thisY, psym=10, color=255, thick=3.0
oploterror, xpos, thisY, thisErr,thisErr, thick=3.0, psym=2, color=50, errthick=3.0, errcolor=50

id=11
rgi=[4, 5, 6, 7, 8, 9, 10, 11]
thisY=meanWL[rgi, id]
thisErr=sigWL[rgi, id]
;oplot, xpos, thisY, psym=10, color=255, thick=3.0
oploterror, xpos[rgi], thisY, thisErr,thisErr, thick=3.0, psym=10, color=0, errthick=3.0, errcolor=0

 
badx1=[4.0, 16.0]
badx2=[79.0, 91.0]
yband1=[100.0, 100.0]
yband2=[-100.0, -100.0]

;Re-arrange so that the bands overplot some of the points, but not others.
oband,badx1,yband1,yband2,color=220
oband,badx2,yband1,yband2,color=220

leg_text=['VG02', $
          'VG04', $
          'MIS01', $
          'VG08', $
          'VG05']
legend,leg_text,psym=[5, 6, 4, 2, -3], color=[150,255,100,50,0], $
  position=[24, 0.507], charthick=2.0
  
loadct,0
mT=0.5275
sigT=0.0005482


xband1=[mT-sigT, mT-sigT]  
xband2=[mT+sigT, mT+sigT]
oband, [-1000, 1000], xband1, xband2, color=200
oplot, [-1000, 1000], [mT, mT], linestyle=2, color=0

XYOUTS, 7.0, 0.505, 'blocked by sample holder', ORIENTATION=75.0
XYOUTS, 81.0, 0.505, 'blocked by sample holder', ORIENTATION=75.0
XYOUTS, 30.0, mt+sigT*1.3, 'Fresnel Prediction', ORIENTATION=0.0

Device, /Close_File
Set_Plot, thisDevice
Obj_Destroy, psObject

spawn, 'open /Users/gully/IDLWorkspace/FabryPerot/figs/GapCheckAltSTA.eps'
return, 1

end
