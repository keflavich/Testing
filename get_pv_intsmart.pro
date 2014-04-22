;PRO get_pv_int, hops=hops, evla=evla, snake=snake, fil1=fil1
;PRO get_pv_intsmart, cfile, mfile, rfile, namestr, 
;I'm trying to create a merge conflict
;test!
; this is a different edit

sortx = 1

hops = 1
evla = 0
grsco = 0
higal_fil = 0
fil1 = 0
xsig = 10 ; number of sigma of levels to plot
;always set xfwhm to be the minor axis, parallel with line
xfwhm = 3. ;pixels
yfwhm = 6.
snake = 1
g13=0
hotseg1 = 0
hotseg2 = 0
refpoints = 1
image_contours=1

distype = 3
;distype!  1 = arcseconds, 2 = arcminutes, 
; 3 = physical pc, assuming near kinematic 
; distance, default is arcseconds
if keyword_set(grsco) then begin

    radec = 0
    buffer = 10
    
    cfile = 'grs/g32.0+0.0.fits'
    mfile = 'grs/grs_intint_88_102.fits'

    fits_read, 'grs/grs_nh2_88_102.fits', other, otherhdr
;    fits_read, 'grs/large_grsnh2_bigpix.fits'

    image_contours=1

    lvls = [4e21, 7.92e21, 1.18e22, 1.58e22, $
                1.97e22, 2.36e22, 2.75e22]
    units = '!3K'
    
    ch1 = 439
    ch2 = 510
    
    lf1 = 1
    lf2 = 40

    rfile = 'grs/g32_grs_cut1.txt'
    namestr = 'grs/g32_grs_cut1_intsmart'
    
endif


if keyword_set(higal_fil) then begin
    radec = 0
    buffer = 5
    
    cfile = 'higal/grs_g30.21-0.18.fits'
;    mfile = 'higal/grs_intint_97_112.fits'
;    fits_read, 'higal/grs_nh2_97_112.fits', other, otherhdr

    image_contours=1

    lvls = [1e22, 1.5e22, 2e22, 2.5e22, 3e22, 3.5e22]
    units = '!3K'
   
    lf1 = 1
    lf2 = 40

    if keyword_set(fil1) then begin
        
        rfile = 'higal/grs_g30_fil1.txt'
        namestr = 'higal/g30_grs_fil1_intsmart'
        mfile = 'higal/grs_intint_100_108.fits'
        fits_read, 'higal/grs_nh2_100_108.fits', other, otherhdr
        
        ch1 = 495 ;100 km/s
        ch2 = 533 ;108 km/s
        
   
    endif else begin
    
        rfile = 'higal/grs_g30_fil2.txt'
        namestr = 'higal/g30_grs_fil2_intsmart'
        mfile = 'higal/grs_intint_97_112.fits'
        fits_read, 'higal/grs_nh2_97_112.fits', other, otherhdr

        ch1 = 481 ;97 km/s
        ch2 = 552 ;112 km/s
   
    endelse

   
endif



if keyword_set(hops) then begin
    radec = 0
    buffer = 10

    cfile = 'HOPS_NH3_11_cube_l_10.5_13.5.fits'
    mfile = 'HOPS_NH3_11_mom0_l_10.5_13.5.fits'
    fits_read, mfile, other, otherhdr	 
    units = '!3K'
   
    ;line free channels
    lf1 = 0
    lf2 = 400

    if snake eq 1 then begin
        rfile = 'g11.13-0.12_filline.txt' ;snake
        ch1 = 530               ; 19.2 km/s
        ch2 = 543               ; 27.7 km/s
        namestr = 'G11.13-0.12_intsmart' ;snake
        minr = -0.25
        maxr = 1.25
    endif else begin
        rfile = 'g12.9-0.27_filline.txt'
        ch1 = 540               ; 31.5 km/s
        ch2 = 558               ; 39.2 km/s
        namestr = 'G12.9-0.27_intsmart'
        minr = -0.4
        maxr = 1.75
    endelse
    
    if g13 eq 1 then begin
        rfile = 'g13.28-0.32_filline.txt' 
        ch1 = 550
        ch2 = 568
        namestr = 'G13.28-0.32_intsmart' 
    endif


endif

if keyword_set(evla) then begin
    radec = 1
    buffer = 20

    cfile = 'evla/coldclump_11.cube_r0.5.image.fits'
;    mfile = 'evla/coldclump_11_integrated_masked.fits'
    mfile = 'evla/coldclump_11.cube_r0.5.image.moment.integrated.fits'
    
    if keyword_set(image_contours) then begin
        ;cold
        fits_read, '/orion/data/cara/filaments/pv_analysis/evla/coldclump_11_integrated_masked.fits', $
          other, otherhdr
        lvls = [0.01, 0.02167, 0.03333, 0.045, 0.0567, 0.0683, 0.08]
       

if keyword_set(hotseg1) or keyword_set(hotseg2) then begin 
        ;hot
        fits_read, '/orion/data/cara/filaments/pv_analysis/evla/hotclump_11_integrated_masked.fits', $
          other, otherhdr
        lvls = [0.03, 0.05833, 0.08667, 0.115, 0.1433, 0.1717, 0.2]
endif        

    endif

    units = '!3Jy/beam'
    
    ch1 = 120
    ch2 = 145
;    ch1 = 50 ;to see hyperfine lines
;    ch2 = 250


    ;line free channels
    lf1 = 20
    lf2 = 70

     if keyword_set(fil1) then begin
         minr = -0.01
         maxr = 0.02
         namestr = 'evla/cold_11_fil1_intsmart'
;         rfile = 'evla/cold_11_mom0_fil1.txt'
         rfile = 'evla/cold_11_mom0_fil1_try2.txt'
     endif else begin
         minr = -0.01
         maxr = 0.02
         namestr = 'evla/cold_11_fil2_intsmart'
         ;rfile = 'evla/cold_11_mom0_fil2.txt'
         rfile = 'evla/cold_11_mom0_fil2_try2.txt'
     endelse

     center1 = 0
     if keyword_set(center1) then begin
         minr = -0.01
         maxr = 0.02
         namestr = 'evla/cold_11_center1_intsmart'
         rfile = 'evla/cold_11_mom0_center1.txt'
     endif


     if keyword_set(hotseg1) then begin
         cfile = 'evla/hotclump_11.cube_r0.5_rerun.image.fits'
         mfile = 'evla/hotclump_11.cube_r0.5_rerun.image.moment.integrated.fits'
         ch1 = 115
         ch2 = 139
         namestr = 'evla/hot_11_seg1_intsmart'
         rfile = 'evla/hot_11_mom0_seg1.txt'
         sortx = 0
     endif

    if keyword_set(hotseg2) then begin
         cfile = 'evla/hotclump_11.cube_r0.5_rerun.image.fits'
         mfile = 'evla/hotclump_11.cube_r0.5_rerun.image.moment.integrated.fits'
         ch1 = 130
         ch2 = 150
         namestr = 'evla/hot_11_seg2_intsmart'
         rfile = 'evla/hot_11_mom0_seg2.txt'
         sortx = 0
     endif

endif

fits_read, cfile, data, hdr
fits_read, mfile, mdata, mhdr

sz = size(data)


;region file
; in ds9 draw line (lines) that follow filament and save regions
; in 'image' coordinate
; region file now has format
; hdr
; line(x1, y1, x2, y2)
readcol, rfile, x1, y1, x2, y2, format="(F,F,F,F)" 
;F is plenty good, don't need D

;subtract one pixel, ds9 starts at 1, idl starts at 0
 x1 -= 1 & x2 -= 1
 y1 -= 1 & y2 -= 1 

if sortx eq 1 then oo = sort(x1) else oo = sort(y1)
x1 = x1(oo) & x2 = x2(oo)
y1 = y1(oo) & y2 = y2(oo)

get_line_xy, x1, x2, y1, y2, x, y, npix

;sort x and y
nchan = ch2-ch1 +1
if sortx eq 1 then o = sort(x) else o = sort(y)
x = x(o) & y = y(o)



;for each point make an elliptical gaussian
; major axis orthogonal to line, minor parallel to line
;
; sum over the gaussian

nn = total(npix)
nline = n_elements(npix)

;size of gaussian kernel, 5*fwhm ?
npixker = [5.*xfwhm, 5.*yfwhm]

rotationang = fltarr(nline)
convdata = fltarr(sz(1), sz(2), sz(3))
pv = fltarr(nn, nchan)
interpdata=fltarr(nn, nchan)

tot = 0
pixmax = 0

;for each line
for ii = 0, nline-1 do begin

    nel = npix(ii)
    
    ;determine the angle that is perpendicular to
    ; the defined line
    totprev = tot
    tot = totprev + npix(ii)-1

    delx = float(x(tot) - x(totprev))
    dely = float(y(tot) - y(totprev))

    rad2deg = 360. / (2.*!pi)
    rotationang[ii] = atan( -1. * dely / delx)*rad2deg ;degrees
    ; -1 determined from trial and error.  going clockwise...
    ;symmetric about 0, breaks apart at 90deg.

    ;make a xfwhm, yfwhm 2d gaussian rotated to be perpendicular to
    ; the defined line
    g = make_any_2dgaussian(npixel=npixker, fwhm=[xfwhm, yfwhm], $
                            rotang=rotationang(ii), /normalize)

    ;for each channel, convolve rotated gaussian with data
    ; then interpolate to mid-pixel x, y points for "PV" diagram
    for kk = ch1, ch2 do begin
        convdata[*,*, kk] = convolve_fft(data(*,*,kk), g, /ignore_nan)
        chindex = kk-ch1
        interpdata[*,chindex] = interpolate(convdata(*,*,kk), x, y)
    endfor

  
endfor

;PV is just the interpolated convoluted data vs. channel
pv = interpdata


;get sigma
; in the line free channels at x, y points
;Would be better to do this with a noise map ... in the future ... 
signalfree = fltarr(nn,(lf2-lf1)+1)
for kk = lf1, lf2 do begin
    chindex = kk-lf1
    signalfree[*, chindex] = interpolate(data(*,*,kk), x, y)
endfor
sigma = stddev(signalfree, /nan)


;get velocity and location information
cv3 = sxpar(hdr, 'CRVAL3')
cp3 = sxpar(hdr, 'CRPIX3')
cd3 = sxpar(hdr, 'CDELT3')

chan = findgen(sz(3))
vel = ( (chan - cp3)*cd3 + cv3 ) / 1000. ;km/s

chanpv = chan(ch1:ch2)
velpv = vel(ch1:ch2)

cv1 = sxpar(hdr, 'CRVAL1')
cp1 = sxpar(hdr, 'CRPIX1')
cd1 = sxpar(hdr, 'CDELT1')
glon = ( (x - cp1)*cd1 + cv1)

cv2 = sxpar(hdr, 'CRVAL2')
cp2 = sxpar(hdr, 'CRPIX2')
cd2 = sxpar(hdr, 'CDELT2')
glat = ( (y - cp2)*cd2 + cv2)

;box size for plotting
;buffer = (30/3600.)/abs(cd1) ;30"
boxx1 = min(x)-buffer > 0 & boxx2 = max(x)+buffer < sz(1)-1
boxy1 = min(y)-buffer > 0 & boxy2 = max(y)+buffer < sz(2)-1

;needs to be an integer to plot correctly and such
boxx1 = round(boxx1) & boxx2 = round(boxx2)
boxy1 = round(boxy1) & boxy2 = round(boxy2)

; plus ones are because the x and y ranges are the actual ranges
; plus one because IDL starts at zero and left of the pixel is the
; label.
;   ooooor not?!?!  why not??
; because the line follows a similar point procedure!
glonx1 = (boxx1 - cp1)*cd1 + cv1
glonx2 = (boxx2 - cp1)*cd1 + cv1
glaty1 = (boxy1 - cp2)*cd2 + cv2
glaty2 = (boxy2 - cp2)*cd2 + cv2

;if keyword_set(radec) then euler, glon, glat, glon, glat, 1


;calculate distance along filament
; in arcseconds, arcminutes
; OR physical size based on l, b, v and near assumption
if keyword_set(radec) then begin
    ra1 = glon(0) & dec1 = glat(0)
    raall = glon & decall = glat
endif else begin
    euler, glon, glat, raall, decall, 2
    ra1 = raall(0) & dec1 = decall(0)
endelse
;u of 2 means start with RAx and DCx in degrees, dis in arcseconds
gcirc, 2, ra1, dec1, raall, decall, dis
;get distance (dis) in arcseconds

case distype of
    1: begin
        dis = dis ;arcseconds
        disunit = '!3arcseconds'
    end
    2: begin
        dis = dis/60. ;arcminutes
        disunit = '!3arcminutes'
    end
    3: begin
        if keyword_set(radec) then $
          euler, mean(glon), mean(glat), mglon, mglat, 1 else $
          mglon=mean(glon) & mglat=mean(glat)
        distkpc = kdist(mglon, mglat, (vel(ch2)-vel(ch1))/2.+vel(ch1), /near)/1000. ;dist in kpc
        print, 'distance in kpc is ~', strn(distkpc)
        dis = distkpc*1000. * (dis/3600.) *((2.*!pi)/360.)
        disunit = '!3pc'
    end
    else: begin
        dis = dis ;arcseconds
        disunit = '!3arcseconds' 
    end

endcase


;add points to plot evenly spaced along the filament
if keyword_set(refpoints) then begin
    npoint = 5 ;makes a total of npoint+1 points
    xpoint = fltarr(npoint+1)
    ypoint = fltarr(npoint+1)
    glonpoint = fltarr(npoint+1)
    glatpoint = fltarr(npoint+1)
    dispoint = fltarr(npoint+1)
    for ii = 0, npoint do begin
        dispoint[ii] = (max(dis)/npoint) * ii
        closepoint = min(abs(dispoint(ii)-dis), wclose)
        xpoint[ii] = x(wclose)
        ypoint[ii] = y(wclose)
        glonpoint[ii] = glon(wclose)
        glatpoint[ii] = glat(wclose)
    endfor
endif


if keyword_set(radec) then lnthx = 7 else lnthx = 5
if keyword_set(radec) then lnthy = 6 else lnthy = 5

;make strings for start position and end position
glon1 = strn(glon(0), length=lnthx) & glat1 = strn(glat(0), length=lnthy)
glon2 = strn(glon(nn-1), length=lnthx) & glat2 = strn(glat(nn-1), length=lnthy)
glonmid = strn(glon(nn/2.), length=lnthx) & glatmid = strn(glat(nn/2.), length=lnthy)

if keyword_set(radec) then begin
    pos1 = '!3['+greek('alpha',/ps)+','+greek('delta',/ps)+'] = ' + '[' + glon1 + '!Eo!n,' + glat1 + '!Eo!n]'
    pos2 = '!3['+greek('alpha',/ps)+','+greek('delta',/ps)+'] = ' + '[' + glon2 + '!Eo!n,' + glat2 + '!Eo!n]'
    xunit = '!3Right Ascension [degrees]'
    yunit = '!3Declination [degrees]'
endif else begin
    pos1 = '!3[l,b] = ' + '[' + glon1 + '!Eo!n,' + glat1 + '!Eo!n]'
    pos2 = '!3[l,b] = ' + '[' + glon2 + '!Eo!n,' + glat2 + '!Eo!n]'
    xunit = '!3Galactic Longitude [degrees]'
    yunit = '!3Galactic Latitude [degrees]'
endelse

;Save PV to fits
sxaddpar, pvhdr, 'NAXIS1', nn
sxaddpar, pvhdr, 'NAXIS2', nchan
sxaddpar, pvhdr, 'CDELT2', cd3
sxaddpar, pvhdr, 'CDELT1', cd2
fits_write, namestr+'_pv.fits', pv, pvhdr

ps=1
if keyword_set(ps) then begin
    psfile = namestr+'_pv.eps'
    set_plot, 'ps' 
    device, file=psfile, bits_per_pixel=8, /color, xsize = 24, ysize = 28, $
      /times, font_index=3, /isolatin1
endif

!p.multi=[0,1,2]
!p.font = 0
letter = string(109B)
mu = '!9' + letter +'!X'

!p.thick = 5
!x.thick = 5
!y.thick = 5
!p.charsize=1.3

minr = min(pv);-2.*sigma
maxr = max(pv);-1.*sigma

loadct, 13
;cgimage, pv, minvalue=minr, maxvalue=maxr, /axis,  $
cgimage, pv, minvalue=min(pv), maxvalue=max(pv), /axis,  $
  yrange=[vel(ch1), vel(ch2)], $
  ytitle = '!3Velocity [km/s]', xrange=[min(dis), max(dis)], $
  xtitle = '!3Distance Along Filament ['+disunit+']', $
  background='white', color='black', $;/keep_aspect_ratio, $
  multimargin=[2.5, 8, 10, 2]
;to plot vertical dotted lines to denote the positions for reference
; plotted below.
if keyword_set(refpoints) then begin
    yarr = fill_array(100, vel(ch1), vel(ch2))
    for ii = 0, npoint do begin
        oplot, (fltarr(100)+dispoint(ii)), yarr, $
          color=fsc_color('white'), linestyle=1
        xyouts, dispoint(ii), vel(ch1), $
          strn(ii), color=fsc_color('white'), /data
    endfor
endif

;this and associated programs not on orion for some reason?
;added ALL coyote programs to ~/idllibs/coyote/
cgcolorbar, range=[minr, maxr], title= units, $
  position=[0.15, 0.87, 0.95, 0.93]

;could do contours with smoothpv if i wanted ...
;smoothpv = convolve_res(., 0.4, pv, pvhdr)
if keyword_set(xsig) then begin
    cgcontour, pv, levels=xsig*sigma, /onimage, label=0, $
      color=cgcolor('snow'), thick = 2, charsize=1.2
endif


xyouts, 0.1, 0.81, pos1, $
  color=cgColor("white"), charsize=1.5, /normal
xyouts, 0.69, 0.81, pos2, $
  color=cgColor("white"), charsize=1.5, /normal


;now image too
mdatabox = mdata[boxx1:boxx2, boxy1:boxy2]

loadct, 13
axis_format = {xticks:nxtick};, XTickname:['Cow', 'Pig', 'Dog']}
cgimage, mdatabox, $
  minvalue=min(mdatabox), maxvalue=max(mdatabox), /axis, $
  xrange=[glonx1, glonx2], yrange=[glaty1, glaty2], $
  xtitle = xunit, ytitle = yunit, $
  multimargin=[5.5,10,1,2], /save, /keep_aspect_ratio, $
  axkeywords=axis_format, /axes; xtick = 3 ;ADDED xtick=3 on Nov.21,2013
oplot, glon, glat, color=fsc_color('black'), thick=5

if keyword_set(image_contours) then begin
;image is same size,woot!
    cgcontour, other[boxx1:boxx2, boxy1:boxy2], $
      levels = lvls, $
      color=cgcolor('black'), thick=2, charsize=1.2, $
      label=0, /onimage
endif

if keyword_set(refpoints) then begin
; to plot crosses for reference along the plot
    oplot, glonpoint, glatpoint, color=fsc_color('white'), psym=1, symsize=0.5
    for ii = 0, npoint do $
      xyouts, glonpoint(ii), glatpoint(ii), $
      strn(ii), color=fsc_color('white'), /data, charsize=1.6
endif else begin
    oplot, [glon(0),glon(0)], [glat(0), glat(0)], $
      color=fsc_color('white'), psym=1, symsize=1
    oplot, [glon(nn-1), glon(nn-1)], [glat(nn-1), glat(nn-1)], $
      color=fsc_color('white'), psym=1, symsize=1
endelse

;put the coords in the middle somewhere so more generalized?
; xyouts, 0.15, 0.1, pos1, $
;   color=cgColor("white"), charsize=1.5, /normal
; xyouts, 0.69, 0.1, pos2, $
;   color=cgColor("white"), charsize=1.5, /normal


!p.multi=0
if keyword_set(ps) then begin
    device, /close
    set_plot, 'x'
endif
cgimage, mdata[boxx1:boxx2, boxy1:boxy2], $
  minvalue=min(mdatabox), maxvalue=max(mdatabox), /axis, $
  xrange=[boxx1, boxx2], yrange=[boxy1, boxy2], $
  xtitle = xunit, ytitle = yunit
oplot, x, y, color=fsc_color('black'), thick=5
oplot, xpoint, ypoint, color=fsc_color('white'), psym=1
for ii = 0, npoint do xyouts, xpoint(ii), ypoint(ii), strn(ii), color=fsc_color('white'), /data


save, data, mdata, pv, x, y, glon, glat, $
  glonpoint, glatpoint, xpoint, ypoint, dispoint, npoint, $
  boxx1, boxx2, boxy1, boxy2, glonx1, glonx2, $
  glaty1, glaty2, pos1, pos2, sigma, xsig, minr, maxr, $
  vel, ch1, ch2, mdatabox, xunit, yunit, disunit, units, $
  dis, nn, filename=namestr+'.sav'



END
