
  April 10, 2014

   In everything below, lines that start with ">" are executed from the
   unix/linux command line

 ----  For Mac OSX

   Tested with root_v5.34.14 and root_v5.34.18 built from source on a mac (OSX 10.9.2) configured
   and built with

     > ./configure macosx64 --enable-roofit --disable-xrootd --build=debug
     > make

   For some reason, I have a runtime crash if I use root built without the --build=debug flag.

 ----  At CERN

   I have also tested this at CERN on lxplus5 with root 5.34.10.  I have
   problems getting things to work on lxplus6 (SLC 6).  I also can't get
   the more recent versions of root to work at CERN on either lxplus5 or
   lxplus6.  To run at cern, log in to lxplus5 and set up root like this

     > source /afs/cern.ch/sw/lcg/external/gcc/4.6.3/x86_64-slc5/setup.csh
     > source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.10/x86_64-slc5-gcc46-opt/root/bin/thisroot.csh


===========

   Check out the code and go to the directory

      > git clone https://github.com/owen234/hbb-stats
      > cd hbb-stats

   These two will generate text input files for the likelihood builder from
   reduced trees.  They have the location hardwired to my directory (owen), so
   this will not work for you.  You can skip this step.

      > source macros/setup-input-files3
      > source macros/setup-workspace-files3

   I put the MC text input files for the likelihood builder into a directory
   in the git repository (input-files1/mc-files) so that you don't need to do
   the previous step.  This will build the RooStats workspace from the text
   intput file and save the workspace in a root file.

      > source macros/setup-workspace-files3-from-git-dir

   This does a simple fit and profile likelihood scan using the workspace
   root file as input.  The input files are for the MC test samples.

      > source macros/seutp-fitqual-files3

   This creates the workspace file for the unblind data and does the simple
   fit plus profile likelihood scan for the data.

      > source macros/process-data


===========

   To run the full CLs toys, I use the CERN batch system.  Have a look in the
   cern-batchscripts directory to see what I do.  These will not work for
   you straight out of the box, since you need to change the path to files,
   username, etc...




