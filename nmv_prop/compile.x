#!/bin/csh

set file='rho.den'
set rotfile='rho.den_rho'
set engfile='rho.den_eng'
set esqfile='rho.den_esq'

cat {$file}000 > $file
cat {$file}000_rho > $rotfile
cat {$file}000_eng > $engfile
cat {$file}000_esq > $esqfile

set ith=1

while ($ith <= 180)
  if ($ith < 10) then
    set fix='00'
  else if ($ith < 100) then
    set fix='0'
  else
    set fix=''
  endif
  echo $fix$ith
  cat {$file}$fix$ith >> $file
  cat {$file}{$fix$ith}_rho >> $rotfile
  cat {$file}{$fix$ith}_eng >> $engfile
  cat {$file}{$fix$ith}_esq >> $esqfile
  set ith=`echo "$ith + 1"|bc -l`
end
