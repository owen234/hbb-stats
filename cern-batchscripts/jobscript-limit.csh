#!/bin/tcsh

#-- contents of /afs/cern.ch/sw/lcg/external/gcc/4.6.3/x86_64-slc5/setup.csh ----------------
set gcc_config_version = 4.6.3
set mpfr_config_version = 2.4.2
set gmp_config_version=4.3.2
set LCGPLAT = x86_64-slc5
set LCG_lib_name = lib64
set LCG_arch = x86_64

set LCG_contdir = /afs/cern.ch/sw/lcg/contrib
set LCG_gcc_home = ${LCG_contdir}/gcc/${gcc_config_version}/${LCGPLAT}
set LCG_mpfr_home = ${LCG_contdir}/mpfr/${mpfr_config_version}/${LCGPLAT}
set LCG_gmp_home=${LCG_contdir}/gmp/${gmp_config_version}/${LCGPLAT}

setenv PATH ${LCG_gcc_home}/bin:${PATH}
setenv COMPILER_PATH ${LCG_gcc_home}/lib/gcc/${LCG_arch}-unknown-linux-gnu/${gcc_config_version}

if ($?LD_LIBRARY_PATH) then
setenv LD_LIBRARY_PATH ${LCG_gcc_home}/${LCG_lib_name}:${LCG_mpfr_home}/lib:${LCG_gmp_home}/lib:${LD_LIBRARY_PATH}
else
setenv LD_LIBRARY_PATH ${LCG_gcc_home}/${LCG_lib_name}:${LCG_mpfr_home}/lib:${LCG_gmp_home}/lib
endif
#---------------------------------------------------------------------------------------------


#-- from contents of /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.10/x86_64-slc5-gcc46-opt/root/bin/thisroot.csh --------------

setenv ROOTSYS /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.10/x86_64-slc5-gcc46-opt/root

set path = ($ROOTSYS/bin $path)

if ($?LD_LIBRARY_PATH) then
   setenv LD_LIBRARY_PATH $ROOTSYS/lib:$LD_LIBRARY_PATH      # Linux, ELF HP-UX
else
   setenv LD_LIBRARY_PATH $ROOTSYS/lib
endif

if ($?DYLD_LIBRARY_PATH) then
   setenv DYLD_LIBRARY_PATH $ROOTSYS/lib:$DYLD_LIBRARY_PATH  # Mac OS X
else
   setenv DYLD_LIBRARY_PATH $ROOTSYS/lib
endif

if ($?SHLIB_PATH) then
   setenv SHLIB_PATH $ROOTSYS/lib:$SHLIB_PATH                # legacy HP-UX
else
   setenv SHLIB_PATH $ROOTSYS/lib
endif

if ($?LIBPATH) then
   setenv LIBPATH $ROOTSYS/lib:$LIBPATH                      # AIX
else
   setenv LIBPATH $ROOTSYS/lib
endif

if ($?PYTHONPATH) then
   setenv PYTHONPATH $ROOTSYS/lib:$PYTHONPATH
else
   setenv PYTHONPATH $ROOTSYS/lib
endif
if ($?MANPATH) then
   setenv MANPATH `dirname $ROOTSYS/man/man1`:$MANPATH
else
   setenv MANPATH `dirname $ROOTSYS/man/man1`:$default_manpath
endif


##---------------------------------------------------------------------------------------------

 set sig_mass = $1
 set wsrootfile = $2
 set ntoy = $3
 set minmu = $4
 set maxmu = $5
 set npoints = $6
 set outputdir = $7
 set outfilebase = $8

 set finaldir = /afs/cern.ch/work/o/owen/private/$outputdir
 set fulloutfilebase = /afs/cern.ch/work/o/owen/private/$outputdir/$outfilebase

 mkdir -p $finaldir

 #set logfile = `printf "limit-sigmass-%d.log" $sig_mass`
 set logfile = `printf "%s-job.log" $outfilebase`
 echo logfile is $logfile

 which root

 pwd

 set batchworkdir = `pwd`
 echo batchworkdir is $batchworkdir

 cd /afs/cern.ch/user/o/owen/analysis-code/hbb-stats/

 echo current directory
 pwd

 echo running root
 echo  rundemo.c\(\"$wsrootfile\",0,$maxmu,$minmu,$npoints,$ntoy,\"$fulloutfilebase\"\) 
 root -b -q rundemo.c\(\"$wsrootfile\",0,$maxmu,$minmu,$npoints,$ntoy,\"$fulloutfilebase\"\) >& $finaldir/$logfile

 echo Done
