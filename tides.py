#Tidal correction according to Longman 
#Program converted from TIDE-ACD.BAS and made into a module
#
#Example:
#tides.tides(longitude, latitude, timezone, datetime="2014-10-14 12:01:01")
#

import math
import numpy
import datetime

datalaikas= str(datetime.datetime.now())[:-7]
def tides(lng=-21.79, lamda = 55.3, tz=0, dt=datalaikas):
#    lng = -21.79 # longitude  WEST is positive!
#    lamda = 55.3 # latitude
    h = 0 # elevation, cm; tides are VERY insensitive to elevation changes
#    #a$ = DATE$ #  date$ function returns string of form: mm-dd-yyyy
    dt = datetime.datetime.strptime(dt,"%Y-%m-%d %H:%M:%S")
    dt = dt - datetime.timedelta(hours=tz)
    byear = dt.year #2014#VAL(RIGHT$(a$, 4)) # beginning year
    bmonth = dt.month #10#VAL(LEFT$(a$, 2)) # beginning month
    bday = dt.day #20#VAL(MID$(a$, 4, 2)) # beginning day, UTC
    bhour = dt.hour #11 # beginning hour, UTC
    bminute = dt.minute
    bsecond = dt.second

    floathour = bhour + 1./60. * bminute + 1./3600.*bsecond
#    print " Tidal values will be calculated every hour and stored in the file"
#    print " TIDAL.DAT as pairs of numbers consisting of hour, tidal acceleration"
#    print " (microgals)."

    # Constants

    pi = 3.1415927
    mu = 6.67E-08
    m = 7.3537E+25
    s = 1.993E+33
    il = .08979719
    omega = .4093146162
    ml = .074804
    el = .0549
    cl1 = 1.495E+13
    cl = 3.84402E+10
    al = 6.37827E+08
    #    Love Numbers
    h2 = .59
    k2 = .27


    LoveFactor = (1 + h2 - 1.5 * k2) # w/h2=0.59 & k2=0.27, LoveFactor=1.185

    minc = [0,31,59,90,120,151,181,212,243,273,304,334]

    minit = 0
    day = bday
    month = bmonth
    year = byear
    # algorithm doesn#t work for the first two months of 1900
    #for hrgmt in numpy.arange(xb,xe+1,hrinc):
    hrgmt = floathour
    dday = day + hrgmt / 24.
    tl0 = hrgmt + minit / 60.
    nleap = int((year - 1900) / 4)
    if (year % 4 == 0) and (month < 3):
        nleap = nleap - 1
    xm = minc[month-1]
    tdays = .5 + (year - 1900) * 365 + nleap + xm + (day - 1) + tl0 / 24

    t = tdays / 36525
    n = 4.523601612 - 33.75715303 * t + .0000367488 * t * t + .0000000387 * t * t * t
    el1 = .01675104 - .0000418 * t + .000000126 * t * t
    sl = 4.720023438 + 8399.7093 * t + .0000440695 * t * t + .0000000329 * t * t * t
    pl = 5.835124721 + 71.01800935999999 * t - .0001805446 * t * t - .0000002181 * t * t * t
    hl = 4.881627934 + 628.3319508 * t + .0000052796 * t * t
    pl1 = 4.908229467 + .0300052641 * t + 7.902400000000001E-06 * t * t + .0000000581 * t * t * t
    i = math.acos(.9136975738000001 - .0356895353 * math.cos(n))
    nu = math.asin(.0896765581 * math.sin(n) / math.sin(i))
    L = lng * .0174532925
    tl = (15 * (tl0 - 12) - lng) * .0174532925 #  ?????
    chi = tl + hl - nu
    chi1 = tl + hl
    ll1 = hl + 2 * el1 * math.sin(hl - pl1)
    cosalf = math.cos(n) * math.cos(nu) + math.sin(n) * math.sin(nu) * .9173938078
    sinalf = .3979806546 * math.sin(n) / math.sin(i)
    alf = 2 * math.atan(sinalf / (1 + cosalf))
    xi = n - alf
    sigma = sl - xi
    ll = sigma + .1098 * math.sin(sl - pl) + .0037675125 * math.sin(2 * (sl - pl)) + .0154002735 * math.sin(sl - 2 * hl + pl) + .0076940028 * math.sin(2 * (sl - hl))
    lm = lamda * .0174532925
    costht = math.sin(lm) * math.sin(i) * math.sin(ll) + math.cos(lm) * (((math.cos(.5 * i)) ** 2) * math.cos(ll - chi) + ((math.sin(.5 * i)) ** 2) * math.cos(ll + chi))
    cosphi = math.sin(lm) * .3979806546 * math.sin(ll1) + math.cos(lm) * (.9586969039 * math.cos(ll1 - chi1) + .0413030961 * math.cos(ll1 + chi1))
    c = 1 / math.sqrt(1 + .006738 * (math.sin(lm) ** 2))
    rl = 6.37827E+08 * c + h
    ap = 2.60930776E-11
    ap1 = 1 / (1.495E+13 * (1 - el1 * el1))
    dl = 1 / (1 / cl + ap * el * math.cos(sl - pl) + ap * el * el * math.cos(2 * (sl - pl)) + 1.875 * ap * ml * el * math.cos(sl - 2 * hl + pl) + ap * ml * ml * math.cos(2 * (sl - hl)))
    D = 1 / (1 / cl1 + ap1 * el1 * math.cos(hl - pl1))
    gm = mu * m * rl * (3 * (costht ** 2) - 1) / (dl * dl * dl) + 1.5 * mu * m * rl * rl * (5 * (costht ** 3) - 3 * costht) / (dl ** 4)
    gs = mu * s * rl * (3 * (cosphi ** 2) - 1) / (D * D * D)
    g0 = (gm + gs) * LoveFactor

    return(str(dt)+' '+str(g0 * 1000))

#tides(-21.79, 55.3, 0, dt="2014-10-14 12:00:00")
