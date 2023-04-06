# ################################################
# DJEF 01/01/2009
# Need call to DB for star file if dont have one...
# Star file needs to have at least 3 cols: name, RA, DEC
# Delete stars from tvmark that arent in the frame???
#
# Updated x0_E2V and y0_E2V from COSMOS2.13 release/ DJEF 2009.06.22
# ################################################

procedure iwcs (root, starfile, mkwcs)

   string    root          {prompt="Root name of image (e.g. iff0114)"}
   string    starfile      {prompt="Star file (e.g. 020000-300600.coo)"}
   bool      mkwcs         {prompt="Create WCS CDx_y matrix?"}
   bool      ckwcs         {prompt="Check WCS CDx_y matrix?"}

begin

   int     i, j, k, nstar
   int     xb, yb, crp1, crp2, nxm, nym, crp1_big, crp2_big, offx, offy
   int     nchips = 8
   int     xz[8] = 0,0,0,0,2049,2049,2049,2049
   int     yz[8] = 4097,4097,4097,4097,0,0,0,0
   real    mscale
   real    iroa_d = -46.25
   real    sxz[8] = 1.,1.,1.,1.,-1.,-1.,-1.,-1.
   real    syz[8] = -1.,-1.,-1.,-1.,1.,1.,1.,1.
   real    x0[8], y0[8]
   real    x0_SITe[8] = -61.555,-30.070,1.435,32.910,-30.115,-61.610,32.900,1.362
   real    y0_SITe[8] = 0.250,0.168,0.185,0.135,-61.912,-61.565,-62.070,-62.020
# Updated x0_E2V and y0_E2V from COSMOS2.13 release/ DJEF 2009.06.22
   real    x0_E2V[8] = -61.524,-29.681,2.190,34.043,-29.694,-61.545,34.041,2.159
   real    y0_E2V[8] = 2.014,1.993,2.0,1.964,-59.982,-59.965,-60.018,-59.978
#    real    x0_E2V[8] = -61.524,-29.681,2.190,34.043,-29.694,-61.545,34.041,2.159
#    real    y0_E2V[8] = 2.014,1.993,2.0,1.964,-59.982,-59.965,-60.018,-59.978
   real    x0_oldSITe = -61.555,-30.070,1.435,32.910,-30.115,-61.58,32.90,1.362
   real    y0_oldSITe = 0.250,0.168,0.185,0.135,-61.912,-61.865,-61.088,-62.020
   real    x1i, y1i, z1i, key1i, x1ii, y1ii, x1c, y1c, z1c, key1c
   real    x2i, y2i, z2i, key2i, x2ii, y2ii, x2c, y2c, z2c, key2c
   real    dx, dy, theta1, theta2, theta, theta_d, dxt, dyt
   real    cv01, cv02, cp01, cp02, sz, lim, exptime
   real    xstar[5000,8], ystar[5000,8]
   real    roffset, f, rzp, pixscale, mangle, iangle
   real    beta, droff, dew, dns, temp
   real    ra[5000], dec[5000]
   real    cd1_1, cd1_2, cd2_1, cd2_2
   real    ra1, dec1, rotan_d, dor, pa, crv1, crv2, north
   real    val1, val2, val3
   string  sval1
   string  name[5000]
   string  im, imm, im1
   string  qtmp1, qtmp2, qtmp3, qtmp4, qtmp5, qtmp6
   string  dew_mod, camera, fs, outfile, coof
   bool    mkw, OK, ckw


# Load packages:
  stsdas
  toolbox
  imgtools
  noao
  imred
  irred
  digiphot
  apphot

# Get parameters
   fs = starfile
   im = root
   mkw = mkwcs
   ckw = ckwcs

   im1 = im//"c1.fits"
   if (access(im1) != YES)
      error (0, "Image %s does not exist", im1)

# Check if a graphics window is open:
   iferr { display (im1, 1, >& "dev$null") }
          error (0, "No display server runnng") 

# READ OUT ALL ESSENTIAL HEADER INFO: RA, DEC, ROTANGLE, BINNING, CAMERA

# tmp file for stdout
  qtmp1 = mktemp ("tmpq")
  qtmp3 = mktemp ("tmpq")

  hselect (im1, "ra-d,dec-d,rotangle,dewarori", "yes") | scan(ra1,dec1,rotan_d,dew_mod)
  delete (qtmp3, veri-)
  printf ("RA = %s; DEC = %s; rotan = %f; dew_mod = %s \n", ra1,dec1,rotan_d,dew_mod)


# Binning: 

  if (access("bin") == yes) 
    delete ("bin", ver-)
  if (access("xbin") == yes) 
    delete ("xbin", ver-)
  if (access("ybin") == yes) 
    delete ("ybin", ver-)
  hselect (im1, "binning", "yes", > "bin")
  ! cat bin | awk '{x = substr($1, 1, 1); print x}' > xbin
  type ("xbin") | scan (xb)
  delete ("xbin", ver-)
  ! cat bin | awk '{y = substr($1, 3, 1); print y}' > ybin
  type ("ybin") | scan (yb)
  delete ("ybin", ver-)
  delete ("bin", ver-)
  printf ("Binning: %f, &f \n",xb, yb)

# Dewar Orientation:

  if (access("dew") == yes) 
    delete ("dew", ver-)
  if (access("xdew") == yes) 
    delete ("xdew", ver-)
  hselect (im1, "DEWARORI", "yes", > "dew")
  ! cat dew | awk '{x = substr($1, 1, 3); print x}' > xdew
  type ("xdew") | scan (dew_mod)
  delete ("dew", ver-)
  delete ("xdew", ver-)
  if (dew_mod == 'Nod') {
    print ('Dewar in Nod & Shuffle mode')
    dor = -90.
  } else if (dew_mod== 'Nor') {
    print ('Dewar in Normal mode')
    dor = 0.
  } else {
    error (0,'Dewar mode %s unknown', dew_mod)
  }

# Camera mode:

  if (access("cam") == yes) 
    delete ("cam", ver-)
  if (access("xcam") == yes) 
    delete ("xcam", ver-)
  hselect (im1, "INSTRUME", "yes", > "cam")
  ! cat cam | awk '{x = substr($2, 1, 1); print x}' > xcam
  type ("xcam") | scan (camera)
  delete ("cam", ver-)
  delete ("xcam", ver-)

  if (camera == 'S') {
    printf ("\nf/2 Camera\n")
    mscale = 67.18
    pixscale = 0.202*66.66667/mscale
    sz = 31.
    north = 0.0 - (rotan_d - (90.0+iroa_d)) + dor
    #east = north + 90.
  } else if (camera == 'L') {
    printf ("\nf/4 Camera\n")
    pixscale = 0.111
    mscale = 66.66667
    sz = 16.
    north = 180.0 + (rotan_d - (90.0+iroa_d)) + dor
    #east = north - 90.
  }

  while (north >= 360.)
    north -=360.
  while (north <= -360.)
    north +=360.
  pa = north * 3.14159 / 180.
  printf ("\nPA = %f\n", pa)

  for (j=1; j<=nchips; j+=1) {
    imm = im//"c"//str(j)//".fits"
    print (imm)
    if (camera == 'L') {
      cd1_1 = sxz[j] * xb * pixscale * cos(pa)/3600.
      cd2_2 = syz[j] * yb * pixscale * cos(pa)/3600.
      cd1_2 =-sxz[j] * xb * pixscale * sin(pa)/3600.
      cd2_1 = syz[j] * yb * pixscale * sin(pa)/3600.
      x0[j] = x0_SITe[j]
      y0[j] = y0_SITe[j]
    } else {
      cd1_1 = sxz[j] * xb * pixscale * cos(pa)/3600.
      cd2_2 =-syz[j] * yb * pixscale * cos(pa)/3600.
      cd1_2 =-sxz[j] * xb * pixscale * sin(pa)/3600.
      cd2_1 =-syz[j] * yb * pixscale * sin(pa)/3600.
      x0[j] = x0_E2V[j]
      y0[j] = y0_E2V[j]
    }
    nxm = 8800/xb
    nym = 8800/yb
    offx = 4240
    offy = 4240
    crp1_big=nxm/2
    crp2_big=nym/2
    crp1 = (crp1_big*xb-offx-x0[j]*mscale-xz[j])/(xb*sxz[j])
    crp2 = (crp2_big*yb-offy-y0[j]*mscale-yz[j])/(yb*syz[j])

# Now update hdrs in all chips IF MKW is set

    if (mkw) {
      printf ("CREATING WCS for chip %i...\n", j)
      hedit (imm, 'CTYPE1', 'RA---TAN', add+, ver-)
      hedit (imm, 'CTYPE2', 'DEC--TAN', add+, ver-)
      hedit (imm, 'crval1', ra1,  up+, add+, ver-)
      hedit (imm, 'crval2', dec1, up+, add+, ver-)
      hedit (imm, 'crpix1', crp1, up+, add+, ver-)
      hedit (imm, 'crpix2', crp2, up+, add+, ver-)
      hedit (imm, 'cd1_1', cd1_1, up+, add+, ver-)
      hedit (imm, 'cd1_2', cd1_2, up+, add+, ver-)
      hedit (imm, 'cd2_1', cd2_1, up+, add+, ver-)
      hedit (imm, 'cd2_2', cd2_2, up+, add+, ver-)

    } else {
      printf ("WCS info already in header...\n")
      hselect (im1, "cd1_1,cd1_2,cd2_1,cd2_2", "yes") | scan (cd1_1,cd1_2,cd2_1,cd2_2)
      hselect (im1, "crval1,crval2", "yes") | scan (ra1,dec1)
      printf ("CONTINUING...\n")
    }

    printf ("cd1_1 = %f \n",cd1_1)
    printf ("cd1_2 = %f \n",cd1_2)
    printf ("cd2_1 = %f \n",cd2_1)
    printf ("cd2_2 = %f \n",cd2_2)
  }


# Read star coordinate file
coof = mktemp("out")
if (access(coof))
    delete (coof, veri-)
if (access("runme_usnob"))
    delete "runme_usnob"

if (access(fs) == YES) {
     fields (fs, fields="2,3,1", print-, > coof)
} else {
    hsel (im1, "exptime", "yes") | scan (exptime)
    lim = 15.0 + 2.5 * log10 (exptime)
    if (lim >=18.0) 
      lim = 18.0
    fs = "usnob.cat"
    printf ("Attempting to download star list from USNO...\n")
    printf ("USNOB1_catalog %s %s %f %f %f %s \n", ra1, dec1, sz, sz, lim, fs)
    printf ("USNOB1_catalog %s %s %f %f %f %s \n", ra1, dec1, sz, sz, lim, fs, > "runme_usnob")
    ! source runme_usnob
# USNOB1_catalog ra1 dec1 sz sz lim fs
    fields (fs, fields="1,2,3", print-, > coof)
}


dxt = 0.0
dyt = 0.0

for (j=1; j<=nchips; j+=1) {
  imm = im//"c"//str(j)//".fits"
  outfile = imm//".coo"
  printf ("%s %s \n",imm, outfile)
  if (access(outfile) == YES)
      delete (outfile, veri-)

# Edit header based on the previous frames edit...
  if (j >=2 && ckw) {
    print int(dxt*sxz[j])
    print int(dyt*syz[j])
    hedit (imm,"CRPIX1","(($)+"//(dxt*sxz[j])//")",veri-)
    hedit (imm,"CRPIX2","(($)+"//(dyt*syz[j])//")",veri-)
  }

# Now, display and tvmark...
  display (imm, j)
  sky2xy(imm, "@"//coof) | fields (fields="5,6,3,7", quit-, > outfile)
#  ! sed '/off/d' outfile > tmp ; mv tmp outfile
  flpr >> qtmp1
  tvmark (j, outfile, color=205, mark="point", pointsize=15)
  fclean >> qtmp1
  gflush

# Do we want to check the WCS???
  if (ckw) {
     print ("OK?")
     scan (OK)
  } else {
    OK = YES
  }

# Correct WCS: Click on a star and its catalog position...

  while (!OK) {
    print('Use any key to select a star on image')
    =fscan(imcur, x1i, y1i, z1i, key1i)
    qtmp2 = mktemp("tmpq")
    print(x1i, y1i, >> qtmp2)
    center(imm,coords=qtmp2,calgori="centroid",cbox=5,interac-,verify-,verbose+)
    center(imm,coords=qtmp2,calgori="centroid",cbox=5,interac-,verify-,verbose+) | fields (fields="2,3") | scan (x1ii,y1ii)
    delete (qtmp2, veri-)
    print(' ')
    print(' ')
    print('Use any key to select catalogue position for star')
    =fscan(imcur, x1c, y1c, z1c, key1c)
# shall we refine this by selecting the closest star in the list???
# ???
    print('got ',x1ii,y1ii,x1c,y1c)
    dx = (x1ii - x1c)
    dy = (y1ii - y1c)
    dxt = dxt + dx*sxz[j]
    dyt = dyt + dy*syz[j]

# # second star:
#   print(' ')
#   print(' ')
#   print(' ')
#   print('Use any key to select 2nd star on image -- far from 1st!!!')
#   =fscan(imcur, x2i, y2i, z2i, key2i)
#   qtmp2 = mktemp("tmpq")
#   print(x2i, y2i, >> qtmp2)
#   center(imm,coords=qtmp2,calgori="centroid",cbox=5,verif-,interac-,verbose+)
#   center(imm,coords=qtmp2,calgori="centroid",cbox=5,verif-,interac-,verbose+) | fields (fields="2,3") | scan (x2ii,y2ii)
#   delete (qtmp2, veri-)
#   print(' ')
#   print(' ')
#   print('Use any key to select catalogue position for star')
#   =fscan(imcur, x2c, y2c, z2c, key2c)
#   print('got ',x2ii,y2ii,x2c,y2c)
#   dx = ((x1ii - x1c) + (x2ii - x2c)) / 2.
#   dy = ((y1ii - y1c) + (y2ii - y2c)) / 2.
#   printf ("dx = %f ; dy = %f \n", dx, dy)
#   theta1=atan((y2ii-y1ii),(x2ii-x1ii))  
#   theta2=atan((y2c-y1c),(x2c-x1c))  
#   theta=theta1-theta2
#   theta_d = theta * 180./3.14159
#   printf ("Theta = %f  -- NOT APPLIED YET!!!\n", theta_d)

# Update header:

    hedit (imm,"CRPIX1","(($)+"//int(dx)//")",veri-)
    hedit (imm,"CRPIX2","(($)+"//int(dy)//")",veri-)
  
    if (access(outfile) == YES)
       delete (outfile, veri-)
    display (imm, j)
    sky2xy(imm, "@"//coof) | fields (fields="5,6,3,7", quit-, > outfile)
#    ! sed '/off/d' outfile > tmp ; mv tmp outfile
    flpr >> qtmp1
    tvmark (j, outfile, color=205, mark="point", pointsize=15)
    fclean >> qtmp1
    gflush

    print ("OK?")
    scan (OK)
    }
  }

if (access(qtmp1))
   delete (qtmp1, veri-)
if (access(coof))
   delete (coof, veri-)

end
