#!/bin/csh

# Process in the outer .fid directory
# 
# Optimize phase, etc by running first on your favorite directory:
# ./HSQC.com favorite.fid
# Then run on all the files:
# ./HSQC.com *.fid

#mkdir processed

foreach i ( $* )

set filename = `echo $i | awk -F '.' '{print $1}'`
echo $filename


var2pipe -in ./$i/fid \
 -noaswap  \
  -xN              2048  -yN               256  \
  -xT              1024  -yT               128  \
  -xMODE        Complex  -yMODE      Rance-Kay  \
  -xSW        11261.261  -ySW         2106.962  \
  -xOBS         799.7142045  -yOBS          81.0435228  \
  -xCAR           4.6723160  -yCAR         118.0544875  \
  -xLAB              HN  -yLAB             N15  \
  -ndim               2  -aq2D          States  \
  -out ./$i/test.fid -verb -ov


nmrPipe -in ./$i/test.fid                                  \
| nmrPipe -fn SOL      \
| nmrPipe -fn SP -off 0.45 -end 0.98 -c 1 -pow 2  \
| nmrPipe -fn ZF -auto      \
| nmrPipe -fn FT          \
| nmrPipe -fn PS -p0 0.0  -p1 0 -di    \
#| nmrPipe -fn EXT -left -sw \
| nmrPipe -fn EXT -x1 6.4ppm -xn 9.6ppm -sw    \
| nmrPipe -fn TP     \
#| nmrPipe -fn LP -fb      \
| nmrPipe -fn SP -off 0.45 -end 0.98 -c 1 -pow 2   \
| nmrPipe -fn ZF -auto  \
| nmrPipe -fn FT         \
| nmrPipe -fn PS -p0 0 -p1 180  -di   \
#| nmrPipe -fn EXT -x1 106ppm -xn 133ppm -sw    \
| nmrPipe -fn POLY -auto \
| nmrPipe -fn TP \
| nmrPipe -fn POLY -auto \
-out ./$i/test.ft2 -ov -verb 

pipe2ucsf ./$i/test.ft2   $filename.ucsf

end
