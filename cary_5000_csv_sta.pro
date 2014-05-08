function cary_5000_csv_STA, in_fn, asc_template
;This function takes in the Cary 5000 formatted csv file, and an ascii_template
;The function outputs a structure with wavelength and data array

fn1=in_fn
at=asc_template

d=read_ascii(fn1, template=at)

;Manipulate the array.

;a) The first two lines are header info.
hdr=d.field001[*, 0]  
gi=where(hdr ne '')
ti=where(hdr eq '')
names=hdr[gi]

;b) Round up the numbers in d.field01
;find what the data sampling is:

wl_start=1350.0
wl_stop=1250.0
wl_interval=2.0

n_wls= (wl_start-wl_stop) / wl_interval
dat=float(d.field001[ti, 2:2+n_wls])

;c) While we're at it let's figure out the collection times:


ids=string(indgen(n_elements(names)))
symb=replicate('&', n_elements(names))
symb2=replicate('\\', n_elements(names))
print, '------------------------------------'
print, '----------Collection times----------'
print, transpose([[ids],[symb],[names], [symb2] ])
print, '------------------------------------'
print, '------------------------------------'

;d) We only need one wavelength array:
wl=d.field001[0, 2:2+n_wls]

;e) Transpose, so wavelength is increasing
wl=reform(wl)
wl=reverse(wl)
dat2=reverse(dat, 2)

ret_struct={wl:wl, dat:dat2}

return, ret_struct

end