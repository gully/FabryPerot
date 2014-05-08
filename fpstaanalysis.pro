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
;Make a correction for the baseline, which drifted by 2% over 1 hour.
baseline2=d.dat[0, *]

print, 1

end
