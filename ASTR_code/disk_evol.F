c     ==============================================================
c        PROGRAM FOR TIME-DEPENDENT EVOLUTION OF PROTOSTELLAR / 
c               PROTOPLANETARY ACCRETION DISK WITH RADIAL
c                        MIXING OF DUST SPECIES 
c                (C) July 2004  by  C.P. Dullemond
c
c     **** BUG?? 28-03-06 **** I do not switch to gas surface
c     density when switching to gas opacity...
c     
c     ==============================================================

#include "configure.h"

#define THEEPS 1.d0

c     --------------------------------------------------------------
c                             MAIN ROUTINE
c     --------------------------------------------------------------
      program diskevol
      implicit none
      integer ir,iz,ith,iformat,idust,nr,itemp,idum,ip,i,j,nroot
      integer do_diskevol,do_traj,read_temp,do_tau,write_chem,trajmode
      integer do_cdh
      doubleprecision tdum,dummy,dum
      logical fex,fex2
      character*15 filename
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
c
c     Message
c
      write(*,*) 'VERSION WITH UNFIXED BUG IN GAS OPTICAL DEPTH'
      write(*,*) 'SEE MARKER **** BUG?? 28-03-06 **** IN CODE.'
c
c     Get natural constants.
c
      call set_natconst()
c
c     Default
c
      global_gtd        = 100. ! The default gas-to-dust ratio (not important!)
      global_flang      = 0.05 ! Fixed flaring angle
      global_mugas      = 2.3  ! Fixed value for mugas
      global_prandtlinv = 1.d0 ! Default value for the inverse Prandtl number
      isw_qvisc         = 1    ! Switch on the viscous heating
      isw_qirrad        = 1    ! Switch on the irradiative heating
      isw_growstar      = 0    ! Switch on the star mass growth
      isw_incmdisk      = 0    ! Switch on the inclusion of mass of disk in mstar
      ml_implicit       = 0    ! Switch for mass loading in implicit way
      global_ginst_plaw = 7.d0 ! The default p value for the grav-inst (1-Q)^p 
      cloud_ismooth     = 0    ! =1 ---> Smooth the end of the infall phase
      cloud_smoothparam1= 5d-1 ! Smoothing parameter
      cloud_smoothparam2= 100.*au ! Smoothing parameter
      global_visheat_mode=0    ! How to compute T_c from T_eff for active accr.
      cloud_idistr      = 0    ! Default: Ulrich way of distributing matter
      do_diskevol       = 1    ! Default: do the disk evolution run
      do_traj           = 0    ! Default: do not do trajectories
      isw_lumbl         = 0    ! Default: no irradiation by boundary layer
      isw_selfirr       = 0    ! Default: no self-irradiation
      global_stlumgrow  = 0.d0 ! If >0, then grow Lstar with Mstar
      ml_angmom_mode    = 0    ! 0=No correction to angular momentum eq
c                              !   for subkepler infall onto the disk.
c                              !   (NOTE: cloud_idistr.eq.1 mimics a solution
c                              !          to this problem; see Hueso & Guillot)
c                              ! -1=Correct using Ruud's empirical method
c                              !  1=Correct using true v0 correction assuming
c                              !    that infalling matter instantly mixes with
c                              !    prior material and remains Kepler.
c
c     Read the input file
c
      open(unit=1,file='diskevol.inp')
      read(1,*) iformat
      read(1,*) global_mstar
      global_mdisk = 0d0
      read(1,*) global_rstar
      read(1,*) global_tstar
      global_lstar = 4d0*pi*ss*global_tstar**4*global_rstar**2
      if(iformat.ge.7) then
          read(1,*) global_stlumgrow
      endif
      read(1,*) global_sigma_min
      read(1,*) global_temp_min
      read(1,*) global_flang
      read(1,*) time_dtrel
      read(1,*) time_dtmin
      read(1,*) time_end
      read(1,*) tdum
      if(tdum.lt.0.d0) read(1,*) time_tsave
      read(1,*) time_maxstep
      read(1,*) time_isave
      read(1,*) tempbg
      read(1,*) global_gravinst
      if(iformat.ge.4) then
          read(1,*) global_ginst_plaw
      endif
      read(1,*) dpr_mode
      if((dpr_mode.ge.1).and.(dpr_mode.lt.100)) then
          read(1,*) dpr_tempcryst
      endif
      if(dpr_mode.ge.100) then
          read(1,*) dpr_nu
          read(1,*) dpr_tempcryst
      endif
      if(iformat.ge.3) then
          read(1,*) idust_evap
          read(1,*) dust_tevap
      endif
      if(iformat.ge.5) then
          read(1,*) global_visheat_mode
      endif
      if(iformat.ge.2) then
          read(1,*) cloud_mass
          read(1,*) cloud_asound
          read(1,*) cloud_omega
      endif
      if(iformat.ge.6) then
          read(1,*) cloud_idistr
          if(cloud_idistr.eq.3) then
              ml_angmom_mode = -1       ! For backward compatibility with Ruud:
c                                       ! For idistr==3 we set angmom_mode=-1
c                                       ! by default but this can be modified
c                                       ! below if iformat.ge.8.
              read(1,*) do_diskevol
              if(do_diskevol.lt.0.or.do_diskevol.gt.1) then
                  write(*,*) 'Invalid do_diskevol:',do_diskevol
                  stop
              endif
              read(1,*) do_traj
              if(do_traj.lt.0.or.do_traj.gt.1) then
                  write(*,*) 'Invalid do_traj:',do_traj
                  stop
              endif
              read(1,*) write_vertden
              if(write_vertden.lt.0.or.write_vertden.gt.1) then
                  write(*,*) 'Invalid write_vertden:',write_vertden
                  stop
              endif
              if(do_traj.eq.1) then
                  read(1,*) trajmode
                  read(1,*) read_temp
                  read(1,*) do_tau
                  read(1,*) do_cdh
                  read(1,*) nparc
                  read(1,*) write_traj
                  read(1,*) write_chem
                  if(read_temp.lt.0.or.read_temp.gt.1) then
                      write(*,*) 'Invalid read_temp:',read_temp
                      stop
                  endif
                  if(do_tau.ne.0.and.do_tau.ne.2) then
                      write(*,*) 'Invalid do_tau:',do_tau
                      stop
                  endif
                  if(do_tau.eq.2.and.read_temp.ne.1)
     $                                                              then
                      write(*,*) 'do_tau is 2, so read_temp must be 1'
                      stop
                  endif
                  if(do_cdh.lt.0.or.do_cdh.gt.1) then
                      write(*,*) 'Invalid do_cdh:',do_cdh
                      stop
                  endif
                  if(do_cdh.eq.1.and.read_temp.ne.1)
     $                                                              then
                      write(*,*) 'do_cdh is 1, so read_temp must be 1'
                      stop
                  endif
                  if(write_traj.lt.0.or.write_traj.gt.1) then
                      write(*,*) 'Invalid write_traj:',write_traj
                      stop
                  endif
                  if(write_chem.lt.0.or.write_chem.gt.1) then
                      write(*,*) 'Invalid write_chem:',write_chem
                      stop
                  endif
              endif
          endif
          read(1,*) isw_lumbl
          read(1,*) isw_selfirr
      endif
      if(iformat.ge.8) then
         read(1,*) ml_angmom_mode
         if(ml_angmom_mode.lt.-1.or.ml_angmom_mode.gt.1) then
             write(*,*) 'Invalid ml_angmom_mode:',ml_angmom_mode
             stop
         endif
      endif
      read(1,*,end=110) isw_qvisc
      read(1,*,end=110) isw_qirrad
      read(1,*,end=110) global_prandtlinv
      read(1,*,end=110) isw_growstar
      read(1,*,end=110) isw_incmdisk
      read(1,*,end=110) ml_implicit
      read(1,*,end=110) cloud_ismooth
      read(1,*,end=110) cloud_smoothparam1
      read(1,*,end=110) cloud_smoothparam2
 110  continue
      close(1)
c
c     Warnings
c
      if(global_visheat_mode.ne.0) then
          write(*,*) '--------------------------------------'
          write(*,*) 'WARNING: Non-standard viscous heating!'
          write(*,*) '--------------------------------------'
      endif
      if((ml_angmom_mode.eq.0).and.(cloud_idistr.eq.3)) then
         write(*,*) 'WARNING: If you use Ruud Visser (2008/9) '
         write(*,*) '         way of distributing infalling matter'
         write(*,*) '         then it is advised to use a '
         write(*,*) '         non-zero angmom_mode to fix the'
         write(*,*) '         problem with sub-Kepler infall.'
      endif
c
c     Errors
c
      if((ml_angmom_mode.ne.0).and.(cloud_idistr.eq.1)) then
         write(*,*) 'ERROR: If you use Hueso & Guillot (2005) '
         write(*,*) '       way of distributing infalling matter'
         write(*,*) '       then using a non-zero angmom_mode '
         write(*,*) '       would be inconsistent.'
         stop
      endif
c
c     Read the radial grid.
c
      open(unit=1,file='rgrid.inp')
      read(1,*) grid_nr
      if(grid_nr.gt.FRSIZE_R) then
          write(*,*) 'Radial grid size too big'
          stop 13
      endif
      do ir=1,grid_nr
          read(1,*) grid_r(ir)
          grid_lnr(ir) = log(grid_r(ir))
      enddo
      do ir=2,grid_nr-1
          grid_dr(ir) = 0.5d0 * (grid_r(ir+1)-grid_r(ir-1))
      enddo
      close(1)
      grid_dr(1) = grid_dr(2)**2 / grid_dr(3)
      grid_dr(grid_nr) = grid_dr(grid_nr-1)**2 / grid_dr(grid_nr-2)
c
c     Read the vertical grid.
c
      dovert = 0
      inquire(file='zgrid.inp',exist=fex)
      if(fex) then
          dovert = 1
          open(unit=1,file='zgrid.inp')
          read(1,*) grid_nz
          if(grid_nz.gt.FRSIZE_Z) then
              write(*,*) 'Vertical grid size too big'
              stop 13
          endif
          do iz=1,grid_nz
              read(1,*) grid_z(iz)
          enddo
          close(1)
          do iz=2,grid_nz-1
              grid_dz(iz) = 0.5d0 * (grid_z(iz+1)-grid_z(iz-1))
          enddo
          grid_dz(1) = grid_dz(2)**2 / grid_dz(3)
          grid_dz(grid_nz) = grid_dz(grid_nz-1)**2 / grid_dz(grid_nz-2)
      endif
      if(cloud_idistr.eq.3.and.dovert.eq.0) then
          write(*,*) 'Must have zgrid.inp for idistr=3.'
          stop
      endif
c
c     Read the surface density sigma, both the total one and the
c     dust species sigma
c
      open(unit=1,file='sigma.inp')
      read(1,*) nr,global_ndust
      if(nr.ne.grid_nr) then
          write(*,*) 'Grid size sigma.inp not consistent with radius'
          stop 13
      endif
      do ir=1,grid_nr
          read(1,*) dummy
          disk_tprsigma(ir) = 2.d0*pi*grid_r(ir)*dummy
          if(disk_tprsigma(ir).lt.2.d0*pi*grid_r(ir)*
     %         global_sigma_min) then
              disk_tprsigma(ir) = 2.d0*pi*grid_r(ir)*
     %             global_sigma_min
          endif
      enddo
      if(global_ndust.gt.0) then
          do idust=1,global_ndust
              do ir=1,grid_nr
                  read(1,*) dummy
                  disk_tprsigmadust(ir,idust) = 
     %                 2.d0*pi*grid_r(ir)*dummy
              enddo
          enddo
      endif
      close(1)
c
c     Put the subkepler factor to 1., i.e. the infall is by default assumed
c     to have v_phi/v_Kepler=1.d0. In particular infall models this can be 
c     a smaller factor.
c
      do ir=1,grid_nr
         ml_subkep(ir) = 1.d0
      enddo
c
c     Test the gas-to-dust ratio. In particular for the current version
c     of the dust crystallization algorithm, the gas-to-dust ratio is
c     a global fixed number. This is because then we can simply fix the
c     value of the abundance of crystalline silicates, and this can be
c     done within the mixing routine. That makes it much faster.
c     In a later version, where I do the mixing and processing together
c     in a single Henyey algorithm, I can let go of this condition.
c
      if(global_ndust.gt.0) then
          do ir=1,grid_nr
              dummy = 0.d0
              do idust=1,global_ndust
                  dummy = dummy + disk_tprsigmadust(ir,idust)
              enddo
              if(abs(disk_tprsigma(ir)-dummy*global_gtd)/
     %             disk_tprsigma(ir).gt.1d-6) then
                  write(*,*) 'SORRY, gas-to-dust ratio is not correct'
                  stop
              endif
          enddo
      endif
c
c     Read the viscous alpha
c
      open(unit=1,file='alpha.inp')
      read(1,*) nr
      if(grid_nr.ne.nr) then
          write(*,*) 'Grid size of alpha.inp not equal'
          stop 13
      endif
      do ir=1,grid_nr
          read(1,*) global_alpha(ir)
      enddo
      close(1)
c
c     Read the Rosseland mean opacity table for dust 
c
      inquire(file='rossmean_dust.inp',exist=fex)
      if(fex) then
c
c         New way
c
          inquire(file='rossmean.inp',exist=fex)
          if(fex) then
              write(*,*) 'ERROR: Either rossmean.inp or ',
     %                   'rossmean_dust.inp'
              stop
          endif
          open(unit=1,file='rossmean_dust.inp')
          read(1,*) ross_ntemp
          if(ross_ntemp.gt.FRSIZE_TEMP) then
              write(*,*) 'ERROR: Nr of temperatures in the rossmean.inp'
              write(*,*) '       greater than FRSIZE_TEMP...'
              stop
          endif
          do itemp=1,ross_ntemp
              read(1,*) ross_temp(itemp),ross_kappa(itemp,1)
          enddo
          close(1)
      else
c
c         Old way
c
          inquire(file='rossmean_dust.inp',exist=fex)
          if(fex) then
              write(*,*) 'ERROR: Either rossmean.inp or ',
     %                   'rossmean_dust.inp'
              stop
          endif
          write(*,*) 'BACKWARD COMPATIBILITY MODE rossmean.inp...'
          open(unit=1,file='rossmean.inp')
          read(1,*) ross_ntemp
          if(ross_ntemp.gt.FRSIZE_TEMP) then
              write(*,*) 'ERROR: Nr of temperatures in the rossmean.inp'
              write(*,*) '       greater than FRSIZE_TEMP...'
              stop
          endif
          do itemp=1,ross_ntemp
              read(1,*) ross_temp(itemp),ross_kappa(itemp,1)
          enddo
          close(1)
      endif
c
c     Check if rossmean_gas.inp exists. If not, and if dust evaporation
c     is included, then protest...
c
      if(idust_evap.ne.0) then 
          inquire(file='rossmean_gas.inp',exist=fex)
          if(.not.fex) then 
              write(*,*) 'ERROR: Dust evaporation cannot be included'
              write(*,*) ' if you dont have a file rossmean_gas.inp.'
              stop
          endif
          open(unit=1,file='rossmean_gas.inp')
          read(1,*) idum
          if(idum.ne.ross_ntemp) then
              write(*,*) 'ERROR: Nr of temperatures in the '
              write(*,*) '       rossmean_gas.inp not equal '
              write(*,*) '       to rosssmean.inp'
              stop
          endif
          do itemp=1,ross_ntemp
              read(1,*) dum,ross_kappa(itemp,2)
              if(abs(dum/ross_temp(itemp)-1.d0).gt.1d-6) then
                  write(*,*) 'ERROR: Temperature grid of '
                  write(*,*) '       rossmean_gas.inp not equal '
                  write(*,*) '       to that of rossmean.inp'
                  stop
              endif
          enddo
          close(1)
      endif
c
c     Set save time stuff
c
      if(tdum.gt.0.d0) then 
          time_tsave = tdum
          varsdt = .false.
      else
          time_inctsave = abs(tdum)
          varsdt = .true.
      endif
c
c     Install the minimal temperature
c
      do ir=1,grid_nr
          temp_min(ir) = global_temp_min
      enddo
c
c     Fill the arrays from Palla & Stahler (1991) for the stellar radius
c     as a function of stellar mass.
c
      data mps/0d0,0.3d0,1d0,2.8d0,3d0,3.1d0,3.2d0,3.5d0,4d0,4.5d0,5d0,
     $         6d0,6.5d0,7d0,8d0/
      data rps/2.5d0,2.5d0,5d0,5d0,5.5d0,7d0,10d0,10.3d0,10d0,9.5d0,
     $         8.5d0,4.5d0,3.6d0,3.3d0,3.1d0/
      do i=1,NPS
          mps(i) = mps(i) * MS
          rps(i) = rps(i) * RS
      enddo
c
c     Read in the data from D'Antona & Mazzitelli (1994) for the
c     photospheric luminosity and stellar temperature as a function of
c     time and stellar mass.
c
      if(cloud_idistr.eq.3) then
          fex = .false.
          inquire(file='dm94.dat',exist=fex)
          if(fex) then
              open(unit=1,file='dm94.dat')
              do i=1,NDMM
                  read(1,10) dmmass(i)
                  dmmass(i) = dmmass(i) * MS
                  do j=1,NDMT
                      if(i.eq.1) then
                          read(1,11) dmtime(j),dmlum(i,j),dmtemp(i,j)
                          dmtime(j) = dmtime(j) * year
                      else
                          read(1,12) dmlum(i,j),dmtemp(i,j)
                      endif
                      dmlum(i,j) = (10d0**dmlum(i,j)) * LS
                      dmtemp(i,j) = 10d0**dmtemp(i,j)
                  enddo
                  if(i.ne.NDMM) read(1,13)
              enddo
   10         format(3X,F5.2)
   11         format(E8.2,F8.3,F8.3)
   12         format(8X,F8.3,F8.3)
   13         format(1X)
              close(1)
              dmpassmax = 0
          else
              write(*,*) 'File dm94.dat is absent.'
              stop
          endif
      endif
c
c     Read in the data from Terebey, Shu & Cassen (1984) for the
c     velocity and density as a function of time and position. The
c     values here are extracted from a high-resolution numerical
c     solution of the TSC84 equations. With proper integration in linear
c     or logarithmic space, these values provide the same accuracy. The
c     major advantage is that the arrays are much smaller, and thus
c     faster to work with.
c
      if(cloud_idistr.eq.3) then
          fex = .false.
          inquire(file='tsc84.dat',exist=fex)
          if(fex) then
              open(unit=1,file='tsc84.dat')
              read(1,13)
              do i=1,NTSCLOG
                  read(1,20) tscxlog(i),tscv0(i),tscwqlog(i),
     $                       tsca0log(i),tscamlog(i),tscvmlog(i),
     $                       tscvqlog(i),tscaqlog(i)
              enddo
   20         format(8(E14.6))
              read(1,13)
              do i=1,NTSCLIN
                  read(1,21) tscxlin(i),tscvmlin(i),tscvqlin(i),
     $                       tscaqlin(i)
              enddo
   21         format(4(E14.6))
              close(1)
          else
              write(*,*) 'File tsc84.dat is absent.'
              stop
          endif
      endif
c
c     Initialize etastar and etadisk for the luminosity and set the
c     outflow density. These values for etastar and etadisk are formally
c     only appropriate for a 1 MS cloud, but we can use them for a 0.5
c     MS cloud as well.
c
      etastar = 0.5d0
      etadisk = 0.75d0
      mflow = 2d-3 * cloud_mass
      rcloud = GG * cloud_mass / ( 2*cloud_asound**2 )
c
c     Take the inverse of some parameters.
c
      cloud_asound_inv = 1d0 / cloud_asound
      mpeffinv = 1d0 / (mp*global_mugas)
c
c     Destroy old save files, unless we don't want to redo the disk
c     evolution, but only want to compute infall trajectories.
c
      if(cloud_idistr.eq.3.and.do_diskevol.eq.0) goto 99
c
      open(unit=1,file='sigma.dat',status='unknown')
      write(1,*) grid_nr,global_ndust
      close(1)
      if(dovert.eq.1) then
          open(unit=1,file='midden.dat',status='unknown')
          write(1,*) grid_nr
          close(1)
          open(unit=1,file='midprs.dat',status='unknown')
          write(1,*) grid_nr
          close(1)
          open(unit=1,file='vertden.dat',status='unknown')
          write(1,*) grid_nr,grid_nz,dovert
          close(1)
          open(unit=1,file='vertprs.dat',status='unknown')
          write(1,*) grid_nr,grid_nz,dovert
          close(1)
      endif
      open(unit=1,file='velo.dat',status='unknown')
      write(1,*) grid_nr,global_ndust
      close(1)
c
c     Added by Daniel: Dec 13, 2010
c     Envelope velocities in terms of ur, uth, uphi
c     
      if(dovert.eq.1) then
          open(unit=1,file='env_vr.dat',status='unknown')
          write(1,*) grid_nr,grid_nz
          close(1)
          open(unit=1,file='env_vth.dat',status='unknown')
          write(1,*) grid_nr,grid_nz
          close(1)
          open(unit=1,file='env_vphi.dat',status='unknown')
          write(1,*) grid_nr,grid_nz
          close(1)
      endif
      open(unit=1,file='v0.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file='temperature.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file='dcoef.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file='mgrav.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file='visc.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file='alphagrav.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file='mdot.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file='mstar.dat',status='unknown')
      write(1,*)
      close(1)
      open(unit=1,file='mdisk.dat',status='unknown')
      write(1,*)
      close(1)
      open(unit=1,file='msyst.dat',status='unknown')
      write(1,*)
      close(1)
      open(unit=1,file='rstar.dat',status='unknown')
      write(1,*)
      close(1)
      open(unit=1,file='lstar.dat',status='unknown')
      write(1,*)
      close(1)
      open(unit=1,file='tstar.dat',status='unknown')
      write(1,*)
      close(1)
      open(unit=1,file='surface.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file='psi.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file='subkep.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file='mu0.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file='j_outflow.dat',status='unknown')
      write(1,*)
      close(1)
      open(unit=1,file='infall.dat',status='unknown')
      write(1,*)
      close(1)
      open(unit=1,file='infall2.dat',status='unknown')
      write(1,*)
      close(1)
      open(unit=1,file='time.dat',status='unknown')
      write(1,*)
      close(1)
c
c     Do the run now
c
      call do_run()
c
   99 continue
c
c     If cloud_idistr.ne.3 or if we don't want to compute trajectories,
c     skip to the end.
c
      if(cloud_idistr.ne.3.or.do_traj.eq.0) goto 98
c
c     ==================================================================
c              NEW MATERIAL BELOW, NOVEMBER 2008 BY RUUD VISSER
c     So far, the program has computed the collapse of a rotating cloud
c     to form a star and a disk. Now we want to follow parcels of
c     material from a given starting position into the disk and get
c     their density and temperature as function of time.
c     ==================================================================
c
c     Read in the collapse results and get some more initial parameters
c     in addition to what we already have. The parameter for
c     initialize_parcels controls whether we do or do not read in the
c     temperature results from RADMC.
c
      call read_data()
      call initialize_parcels(read_temp,write_chem,do_tau,do_cdh)
c
c     Loop over all parcels and calculate the trajectory for each one.
c
      open(unit=10,file='trajectories.dat',status='unknown')
      write(10,*) 'nparc:',nparc
      open(unit=12,file='diskentry.dat',status='unknown')
      write(12,*) 'nparc:',nparc
      do ip=1,nparc
          if(nparc.lt.100000) write(10,*) 'ip:',ip
          if(mod(ip,500).eq.0) write(*,*) 'ip:',ip
!           if(mod(ip,50).eq.0) write(*,*) 'ip:',ip
!           write(*,*) 'ip: ', ip  
          if(write_traj.eq.1) then
              write(filename,119) 'traj.',ip-1,'.dat'
  119         format(A5,I6.6,A4)
              open(unit=20,file=filename,status='unknown')
!               write(20,*) ip
          endif
!           if(ip.gt.40) stop
c
c  New Material Sep 9: DH
c  If ip > 10000: -> no need for trajectories, these parcels are 
c  for the abundance structure of the model
c
          if (ip.le.50000) then
            if(trajmode.eq.1) then
                call make_trajectory(ip,write_chem,do_tau,do_cdh)
            elseif(trajmode.eq.2) then
                call make_trajvode(ip,write_chem,do_tau,do_cdh)
            elseif(trajmode.eq.3) then
                call make_trajvodeback(ip,write_chem,do_tau,do_cdh)
            else
                write(*,*) 'Invalid trajmode:',trajmode
                stop
            endif
!             stop
          else
            call make_evolvepts(ip,write_chem,do_tau,do_cdh)
          endif
          if(write_traj.eq.1) then
              close(20)
          endif
      enddo
      close(10)
      write(12,*) -1,0d0,0d0,0d0
      close(12)
      open(unit=1,file='nncoll.dat',status='unknown')
      write(1,*) nncoll+2
      close(1)
c
   98 continue
      end



c     ------------------------------------------------------------------
c     This subroutine sets some natural constants.
c     ------------------------------------------------------------------
      subroutine set_natconst()
c
      implicit none
c
#include "common_natconst.h"
c
c     Set natural constants.
c
      pi = 3.14159265358979324d0
      mp = 1.6726d-24                 !proton mass [g]
      GG = 6.672d-8                   !gravitational constant [cm^3 g^-1 s^-2]
      MS = 1.9891d33                  !solar mass [g]
      RS = 6.96d10                    !solar radius [cm]
      LS = 3.827d33                   !solar luminosity [erg s^-1]
      au = 1.4959787d13               !astronomical unit [cm AU^-1]
      year = 31556926d0               !year [s yr^-1]
      kk = 1.38065d-16                !Boltzmann constant [erg K^-1]
      ss = 5.670399d-5                !Stefan-Boltzmann constant [g K^-4 s^-3]
      m0 = 0.975d0                    !reduced mass (S77,TSC84) [-]
c
c     Also include some inverse values (1/x), for faster computations
c     later (multiplication is faster than division).
c
      piinv = 1d0 / pi
      GGinv = 1d0 / GG
      yearinv = 1d0 / year
c
c     And some multiples.
c
      year7d4 = 7d4 * year
      m0twosq = sqrt(2d0*m0)
      m0halfsq = sqrt(0.5d0*m0)
      m0half3p5o3 = (0.5d0*m0)**3 * sqrt(0.5d0*m0) / 3d0
c
c     Precompute 1/3, which we use a lot later on.
c
      f1o3 = 1d0/3d0
c
      return
      end
c
c
c
c     ------------------------------------------------------------------
c     This subroutine sets the initial parameters for the collapse and
c     disk formation. If the readtemp flag is active, temperatures are
c     read in from RADMC.
c     ------------------------------------------------------------------
      subroutine initialize_parcels(readtemp,write_chem,do_tau,do_cdh)
c
      implicit none
      integer readtemp,write_chem,do_tau,do_cdh,i,j,k,ndum,ir,iz,it,i1
      doubleprecision dum
      logical fex
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
#include "common_chem.h"
c
c     Set model parameters.
c
      mdot = m0*cloud_asound**3/GG !mass accretion rate [g s^-1]
      rcloud = GG * cloud_mass / ( 2*cloud_asound**2 ) !envelope radius [cm]
      tcloudm = (rcloud/cloud_asound)*(2./m0) !collapse time [s]
c      tcloudm = (GG*cloud_mass) / 
c     %        (cloud_asound*cloud_asound*cloud_asound*m0)
      zerotemp = 10d0           !initial temperature [K]
      precoll = 3d5 * year      !stationary time before collapse begins [s]
      tmax = tcloudm    !maximum time for trajectory calculations [s]
c      nparc = 12000             !number of parcels [-]
      warning = 0
      nncoll = 0
c
c     Read in the initial positions for our parcels.
c
      fex = .false.
      inquire(file='r0th0.dat',exist=fex)
      if(fex) then
          open(unit=1,file='r0th0.dat')
          read(1,*) ndum
          if(ndum.ne.nparc) then
              write(*,*) 'Number of parcels in r0th0.dat is incorrect.'
              stop
          endif
          do i=1,nparc
              read(1,*) r0(i),th0(i),t0(i)
          enddo
          close(1)
      else
          write(*,*) 'File r0th0.dat is absent.'
          stop
      endif
c
c     Where does tmax fall in the diskevol time steps?
c
      imax = 1
      call hunt(de_time,de_ntime,tmax,imax)
      if(imax.le.0.or.imax.ge.de_ntime) stop 8390
c
      if(readtemp.eq.1) then
c
c     Read in the time steps for the temperature results.
c
          fex = .false.
          inquire(file='radmctime.dat',exist=fex)
          if(fex) then
              open(unit=1,file='radmctime.dat')
              read(1,31) nr_time_tdust
   31         format(I4)
              if(nr_time_tdust.gt.NTDUSTMAX) stop 2056
              do i=1,nr_time_tdust
                  read(1,30) time_tdust(i)
c                  write(*,*) time_tdust(i)
              enddo
   30         format(E16.8)
              close(1)
          else
              write(*,*) 'File radmctime.dat is absent.'
              stop
          endif
          dtime_tdust_inv(nr_time_tdust) = 0d0
          do i=1,nr_time_tdust-1
              dtime_tdust_inv(i) = 1d0 / (time_tdust(i+1)-time_tdust(i))
          enddo
c
c     Read in RADMC's temperature results. Note that these data are
c     actually the 10-logs of the temperature.
c
          fex = .false.
          inquire(file='radmctemp.dat',exist=fex)
          if(fex) then
              write(*,*) "Reading in RADMC's temperature profiles..."
              open(unit=1,file='radmctemp.dat')
              do i=1,nr_time_tdust
                  read(1,50) ndum
   50             format(I4)
                  if(ndum.ne.i) then
                      write(*,*) i,ndum
                      stop 2053
                  endif
                  do k=1,grid_nr
                      do j=1,grid_nz
c                          write(*,*) j,k
                          read(1,52) dust_temp(k,j,i)
                      enddo
                  enddo
              enddo
   52         format(F5.3)
              close(1)
              write(*,*) "Done with RADMC's temperature profiles..."
          else
              write(*,*) 'File radmctemp.dat is absent.'
              stop
          endif
      else
c
c     If not reading in temperatures, set a flag so the program knows
c     this later on.
c
          time_tdust(1) = -1d0
          write(*,*) "NOT reading in any temperature profiles."
      endif
c
c     Optionally read in the 2D radiation field from RADMC. Note that
c     these data are actually the 10-logs of the radiation field.
c
      if(write_chem.eq.1.and.do_tau.eq.2) then
          fex = .false.
          inquire(file='radmcg0.dat',exist=fex)
          if(fex) then
              write(*,*) "Reading in RADMC's 2D radiation fields..."
              open(unit=1,file='radmcg0.dat')
              do i=1,nr_time_tdust
                  read(1,50) ndum
                  if(ndum.ne.i) then
                      write(*,*) i,ndum
                      stop 2054
                  endif
                  do k=1,grid_nr
                      do j=1,grid_nz
                          read(1,53) rad2d(k,j,i)
                      enddo
                  enddo
              enddo
   53         format(E23.15)
c      Format 53 added by Maria Drozdovskaya July 2013
              close(1)
              write(*,*) "Done with RADMC's radiation fields..."
          else
              write(*,*) 'File radmcg0.dat is absent.'
              stop
          endif
      endif
c
c     Optionally read in the 2D H2 column densities from RADMC. Note
c     that these data are actually the 10-logs of the column densities.
c
      if(write_chem.eq.1.and.do_cdh.eq.1) then
          fex = .false.
          inquire(file='radmccdh.dat',exist=fex)
          if(fex) then
              write(*,*) "Reading in RADMC's 2D H2 column densities..."
              open(unit=1,file='radmccdh.dat')
              do i=1,nr_time_tdust
                  read(1,50) ndum
                  if(ndum.ne.i) then
                      write(*,*) i,ndum
                      stop 2055
                  endif
                  do k=1,grid_nr
                      do j=1,grid_nz
                          read(1,52) cdh(k,j,i)
                      enddo
                  enddo
              enddo
              close(1)
              write(*,*) "Done with RADMC's H2 column densities..."
          else
              write(*,*) 'File radmccdh.dat is absent.'
              stop
          endif
      endif
c
      return
      end
c
c
c
c     --------------------------------------------------------------
c               FUNCTION: FIND ROSSELAND MEAN OPACITY
c     --------------------------------------------------------------
      function rossmean(temp,icmp)
      implicit none
      doubleprecision rossmean,temp
      integer icmp
c
#include "common_disk.h"
c
      doubleprecision eps
      integer jlo
c
c     Find the position of temp in the rosseland mean opac table
c
      call hunt(ross_temp,ross_ntemp,temp,jlo)
c
c     Check if out of bound
c     If so, then return the closest kappa
c
      if(jlo.le.0) then
          rossmean = ross_kappa(1,icmp)
          return
      endif
      if(jlo.ge.ross_ntemp) then
          rossmean = ross_kappa(ross_ntemp,icmp)
          return
      endif
c
c     Return the linear interpol of the kappa table
c
      eps    = (temp-ross_temp(jlo))/
     %         (ross_temp(jlo+1)-ross_temp(jlo))
      if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 83745
      rossmean  = (1.d0-eps)*ross_kappa(jlo,icmp)+
     %            eps*ross_kappa(jlo+1,icmp)
c
      return
      end


c     --------------------------------------------------------------
c               HELPER FUNCTION FOR SOLVE TEMPERATURE
c                    FOR PURE ACCRETIONAL HEATING
c     --------------------------------------------------------------
      function helptemp(temp,r,w)
      implicit none
      integer w
      doubleprecision helptemp,temp,r
c
#include "common_natconst.h"
#include "common_disk.h"
c
      doubleprecision rossmean,kappa,qmin,qtoomre,qcrit
      parameter(qcrit=1.d0)
c
c     Find the rosseland mean opacity for this temperature
c
      if(ist_icmp.le.0) stop 15543
      if(ist_icmp.gt.2) stop 15544
      kappa    = rossmean(temp,ist_icmp)
c
c     Find the Toomre parameter
c
      if(st_fact_toomre.ne.0.d0) then
          qtoomre       = st_fact_toomre * sqrt(temp)
          if(qtoomre.lt.qcrit) then
              if(qtoomre.lt.0.d0) stop 13
              st_alpha_grav = (1.d0-qtoomre)**global_ginst_plaw
          else
              st_alpha_grav = 0.d0
          endif
      else
          st_alpha_grav = 0.d0
      endif
c
c     Find the qplus for the accretion
c
      st_qplus = st_fact_qpl * temp * ( st_alpha_mri + st_alpha_grav )
c
c     Find the qmin 
c
      if(global_visheat_mode.eq.0) then
c
c         The usual purely radiative cooling through an optically thick
c         disk
c
c         **** BUG?? 28-03-06 ****
c         It seems that the st_sigdust must be replaced by st_sigmagas
c         or so when the gas opacity is used. This is not yet fixed...
c
          qmin     = 8.*ss*temp**4 / (3.*st_sigdust*kappa)
      elseif(global_visheat_mode.eq.-1) then
c
c         Absolute extreme case: T_c = T_eff. So ignore the
c         blanket-effect
c
          qmin     = 2*ss*temp**4
      else
          stop 33543
      endif
c
c     Return the differences of the logs
c
      helptemp = log(st_qplus) - log(qmin)
      return
      end


c     --------------------------------------------------------------
c                   FUNCTION: COMPUTE SELF-RADIATION
c
c     A very simple self-irradiation recipe. We assume that the disk 
c     is optically thick. 
c     --------------------------------------------------------------
      function find_temp_selfirr(ir)
      implicit none
c
#include "common_natconst.h"
#include "common_disk.h"
c
      doubleprecision flux,find_temp_selfirr
      integer ir,ir2
c
      flux = 0.d0
      if(ir.gt.2) then 
          do ir2=2,ir-1
              flux = flux + ss *
     %           (global_flang*(1.d0-(grid_r(ir2)/grid_r(ir))))**2 *
     %           0.5 * ( disk_teff(ir2-1)**4 + disk_teff(ir2)**4 ) *
     %           ( grid_r(ir2)**2 - grid_r(ir2-1)**2 ) / 
     %           grid_r(ir)**2
          enddo
      endif
      find_temp_selfirr = (flux/ss)**0.25
c
      end


c     --------------------------------------------------------------
c                      FUNCTION: SOLVE TEMPERATURE 
c
c     This routine solves the local disk temperature including the
c     heating via viscous dissipation as well as the heating by
c     irradiation by the central star. As an important by-producy
c     this routine also computes the viscosity coefficient alpha
c     which may have been modified by gravitational instability.
c     --------------------------------------------------------------
      function solvetemp(ir,r,tprsigma,temp0,w)
      implicit none
      doubleprecision solvetemp,tprsigma,temp0,mu,r
      integer ir,w
c
      doubleprecision sigma,zbrent,helptemp,temp_visc,lumbl,rdx
      doubleprecision temp_irrad,temp_selfirr,temp_blirr
      doubleprecision eps,dummy,dm1,dm2,find_temp_selfirr
      external helptemp
c
#include "common_natconst.h"
#include "common_disk.h"
c
      mu             = global_mugas
c
c     First check if the ir indeed belongs to the r
c     Just an internal consistency check.
c     
      if((ir.lt.1).or.(ir.gt.grid_nr)) stop 88324
      if(ir.lt.grid_nr) then
          if((r.lt.grid_r(ir)).or.(r.ge.grid_r(ir+1))) stop 88325
      endif
c
c     Compute the basic alpha at location r by linear interpolation
c     Note: the real alpha (including grav inst effects) is computed
c     on the fly during the temperature solving.
c
      if(ir.ne.grid_nr) then
          eps    = (r-grid_r(ir))/(grid_r(ir+1)-grid_r(ir))
          st_alpha_mri = (1.d0-eps)*global_alpha(ir)+
     %                   eps*global_alpha(ir+1)
      else
          st_alpha_mri = global_alpha(grid_nr)
      endif
c
c     Make sure to use the right Rosseland mean table for the 
c     temperature determination
c
      ist_icmp = 1
c
c     First the temperature for accretion only
c
c      if(w.eq.1.and.r.gt.0.2*au) then
c          write(*,*) 'isw_qvisc,idust_evap',isw_qvisc,idust_evap
c      endif
      if(isw_qvisc.ne.0) then
c
c         ...prepare the helper function for zbrent
c
          sigma          = (tprsigma/(2*pi*r))
          st_fact_qpl    = (9./4.)*sigma*
     %                     (kk/(mu*mp))*sqrt(GG*disk_mgrav(ir)/r**3)
          if(global_gravinst.ne.0) then
              st_fact_toomre = sqrt(kk/(mu*mp))*
     %                         sqrt(GG*global_mstar/r**3)/
     %                         (pi*GG*sigma)
          else
              st_fact_toomre = 0.d0
          endif
          st_sigdust     = sigma/global_gtd
c
c         ...call zbrent to find the temperature
c
c         NOTE: Here we'll have to build in the use of the old
c               temperature as initial guess in order to be able to
c               simulate the limit-cycle instability in the disk.
c
          if(helptemp(1d0,0d0,0).gt.0.d0) then 
c#############################
c#############################
              dm1 = helptemp(1d0,0d0,0)
              dm2 = helptemp(1d10,0d0,0)
              if(dm2.gt.0.d0) then
                  temp_visc = 1d10
              else
                  temp_visc = zbrent(helptemp,1d0,1d10,1d-3,w,r)
              endif
c
c             Include dust evaporation here, if requested
c
              if(idust_evap.ne.0) then
                  if(temp_visc.gt.dust_tevap) then 
c
c                     Recompute the temperature with only the gas
c                     
                      ist_icmp = 2
                      dm1 = helptemp(1d0,0d0,0)
                      dm2 = helptemp(1d10,0d0,0)
                      if(dm2.gt.0.d0) then
                          dummy = 1d10
                      else
                          dummy = zbrent(helptemp,1d0,1d10,1d-3,0,0d0)
                      endif
c
c                     If this new temperature is BELOW the t_evap, then
c                     we must avoid going into flip-flop. In fact we
c                     will get a thermostat. Else use the new
c                     temperature.
c                     
                      if(dummy.lt.dust_tevap) then
                          temp_visc = dust_tevap
                      else
                          temp_visc = dummy
                      endif
                  endif
              endif
          else
              temp_visc = 0.d0
          endif
c
c         One more call to the temperature routine, to make sure that
c         we get the current values, both for the self-irradiation 
c         (the st_qplus) and for the gravitational instability alpha.
c         (24-11-05)
c
          if(temp_visc.gt.0.d0) then
              dummy      = helptemp(temp_visc,r,w)
          else
              st_qplus   = 0.d0
          endif
c
      else
          temp_visc = 0.d0
          st_qplus  = 0.d0
      endif
c
c     Now the temperature from the irradiation
c     
c     ...First the irradiation by the star
c
      if(isw_qirrad.ne.0) then 
c
c         Check if we wish the star luminosity to grow with star mass
c
          if((global_stlumgrow.gt.0.d0).and.
     %       (cloud_mass.ne.0.d0)) then
              rdx = (global_mstar/cloud_mass)**global_stlumgrow
              if(rdx.gt.1.d0) rdx=1.d0
          else
              rdx = 1.d0
          endif
c
c         Compute the irradiation temperature
c         
          temp_irrad = global_tstar*(0.5*global_flang*rdx*
     %                 (global_rstar/r)**2)**0.25
      else
          temp_irrad = 0.d0
      endif
c
c     ...Then the irradiation by the disk itself (self-irradiation)
c        NEWSTUFF 24-11-05
c
      if(isw_selfirr.ne.0) then
          temp_selfirr = find_temp_selfirr(ir)
      else
          temp_selfirr = 0.d0
      endif
c
c     ...Then the irradiation by the boundary layer of the disk
c        (or in other words the `accretion shock' on the star)
c        NEWSTUFF 24-11-05
c
      if(isw_lumbl.ne.0) then
          if(isw_lumbl.eq.1) then
c             
c             Use Calvet & Gullbring magnetospheric accretion model
c
              lumbl      = (1.d0-global_rstar/grid_r(1)) * 
     %                     ( GG * abs(disk_mdot(1)) * global_mstar / 
     %                       global_rstar )
              temp_blirr = (0.5*global_flang*lumbl/
     %                      (4*pi*r*r*ss))**0.25
          else
c
c             Simple boundary layer (not yet done)
c
              stop 8264
          endif
      else
          temp_blirr = 0.d0
      endif
c
c     Return the 'sum' of the two 
c     BUGFIX 13-07-05: Add the background temperature as well
c
      solvetemp      = (temp_visc**4+temp_irrad**4+tempbg**4+
     %                  temp_blirr**4+temp_selfirr**4)**0.25
c
c     Compute (for use in self-irradiation) the effective temperature
c     assuming a blackbody disk.
c     NOTE: I omit the background and the self-irradiation here
c     NEWSTUFF 24-11-05
c
      st_tempeff     = ( (st_qplus/(2*ss)) + temp_irrad**4 )**0.25
c
c     If gravitational instability mode active, then compute new alpha
c        NOTE: The renewed call to helptemp() here is only for
c        determining the alpha, and this is not affected by ist_icmp.
c
      if(global_gravinst.ne.0) then
          st_alpha   = st_alpha_mri + st_alpha_grav 
      else
          st_alpha   = st_alpha_mri
      endif
c
      return
      end


c     --------------------------------------------------------------
c                    FUNCTION: SIGMA * NU * SQRT(R)
c
c     This function contains all the microphysics needed for the
c     Shakura-Sunyaev disk. It computes the Sigma * nu * sqrt(R)
c     value, given the radius R (r) and the Sigma*R (tprsigma). The
c     arguments temp and nu are meant for returning the values
c     of temperature and viscosity, in case these are needed. 
c        A crucial property of this function is that the resulting
c     value only depends on r and tprsigma! This is because the 
c     function will be called by the functions that compute the
c     double-log derivative coefficients for the linearized 
c     equations. The 'ir' argument is only there to make it
c     easier to locate where we are on the grid. 
c
c     ARGUMENTS:
c        ir                Index
c        r                 Radius in cm
c        tprsigma          Sigma*2*pi*R in g/cm
c        disk_temp         Initial guess of temperature (usually
c                          this is not necessary, but may become
c                          important for disk instabilities due to
c                          the multiple temperature solutions).
c
c     RETURNS:
c        sigma_nu_sqr      The value of Sigma*nu*sqrt(R)
c        temp              The midplane temperature of the disk
c        nu                The viscosity coefficient
c        global_alpnew     The new value of the alpha with grav inst
c     --------------------------------------------------------------
      function sigma_nu_sqrtr(ir,r,tprsigma,temp,nu,w)
      implicit none
      integer ir,w
      doubleprecision sigma_nu_sqrtr,r,tprsigma,temp,nu
c
#include "common_natconst.h"
#include "common_disk.h"
c
      doubleprecision mu,zbrent,solvetemp
      integer number_invalid
c
      mu    = global_mugas
c
c     Compute the midplane temperature, and as a by-product compute
c     the modified alpha (modified by gravitational instability)
c
      temp  = solvetemp(ir,r,tprsigma,disk_temp(ir),w)
c
c     The viscosity coefficient, using the alpha value returned by 
c     the solvetemp() routine (=st_alpha)
c
      nu    = 3.1957972d11 * st_alpha * r**1.5 * temp /
     %              ( sqrt(disk_mgrav(ir)) * mu )
c#########################################
      if(number_invalid(nu).ne.0) then
          write(*,*) nu,ir,r,temp,disk_mgrav(ir),global_mstar
          stop 801
      endif
c#########################################
c
c     Do a check
c
      if(nu.le.0.d0) then 
          write(*,*) 'ERROR: nu zero or negative...'
          stop
      endif
c
c     Return the final value
c
      sigma_nu_sqrtr = (tprsigma/(2*pi*r)) * nu * sqrt(r)
c
      return
      end



c     --------------------------------------------------------------
c           DERIVATIVE OF LOG(SIGMA NU SQRT(R)) TO LOG(SIGMA*R)
c     --------------------------------------------------------------
      function deriv_signusr_sigr(ir,r,tprsigma,w)
      implicit none
      integer ir,w
      doubleprecision deriv_signusr_sigr,r,tprsigma
c
      doubleprecision eps,sigma_nu_sqrtr,y1,y2,dum1,dum2,dum3,dum4
      parameter(eps=1.d-4)
c
#include "common_disk.h"
c
      y1 = sigma_nu_sqrtr(ir,r,tprsigma,dum1,dum2,w+1)
      y2 = sigma_nu_sqrtr(ir,r,tprsigma*(1.d0+eps),dum3,dum4,w+3)
      deriv_signusr_sigr = 2.d0*(y2-y1)/(eps*(y2+y1))
c
      return
      end


c     --------------------------------------------------------------
c             DERIVATIVE OF LOG(SIGMA NU SQRT(R)) TO LOG(R)
c     --------------------------------------------------------------
      function deriv_signusr_r(ir,r,tprsigma,w)
      implicit none
      integer ir,w
      doubleprecision deriv_signusr_r,r,tprsigma
c
      doubleprecision eps,sigma_nu_sqrtr,y1,y2,dum1,dum2,dum3,dum4
      parameter(eps=1.d-4)
c
      y1 = sigma_nu_sqrtr(ir,r,tprsigma,dum1,dum2,0)
      y2 = sigma_nu_sqrtr(ir,r*(1.d0+eps),tprsigma,dum3,dum4,0)
      deriv_signusr_r = 2.d0*(y2-y1)/(eps*(y2+y1))
c
      return
      end


c     --------------------------------------------------------------
c             COMPUTE RADIAL VELOCITY AND DIFFUSION COEFFICIENT
c
c     This routine computes the radial velocity vr according to the
c     physics of Shakura-Sunyaev disks. It also computes the 
c     effective radial velocity v0 and the diffusion coefficient
c     D that you get when you write the continuity equation out
c     into a diffusion-advection type equation. Because the pure
c     advection form with vr is numerically unstable, we use the
c     diffusion-advection type equation for the model.
c     --------------------------------------------------------------
      subroutine compute_vr_v0_difcoef()
      implicit none
c
#include "common_natconst.h"
#include "common_disk.h"
c
      doubleprecision sigma,mu,tillum,taues,temp_illum,r,tprsigma
      doubleprecision v0corr(FRSIZE_R),signusr,pd_r,pd_sigr,nu,temp
      doubleprecision deriv_signusr_r,deriv_signusr_sigr,sigma_nu_sqrtr
      doubleprecision pd_sigr1,eta_r,fourpi,dum
      doubleprecision x1,x2,y1,y2,dum1,dum2
      integer ir,number_invalid,nfix
      parameter(tillum=3.481278d9)
      parameter(eta_r=2d-3)
      parameter(fourpi=12.5663706143591729538505735331d0)
c
      mu    = global_mugas
      nfix  = 0
c
c     First check something
c
      if(global_sigma_min.eq.0.d0) then
          write(*,*) 'Error: Must define global_sigma_min'
          write(*,*) '       in order to keep log derivs intact'
          stop 13
      endif
c
c     For the numerical model we use v0 and D. These are found as
c     partial derivatives of Sigma * nu * sqrt(R) to (Sigma*R)_{const R}
c     and to (R)_{const Sigma*R}. This is simply using the chain rule.
c     
c                 / dlg(Sigma*nu*sqrt(R)) \           nu
c       v_0 = -3  | --------------------- |           --
c                 \        dlg(R)         /{Sigma*R}  R
c
c                 / dlg(Sigma*nu*sqrt(R)) \                        .
c       D   =  3  | --------------------- |     nu
c                 \      dlg(Sigma*R)     /{R}
c
      do ir=1,grid_nr
          r                 = grid_r(ir)
          tprsigma          = disk_tprsigma(ir) + 
     %                        2*pi*r*global_sigma_min
          signusr           = sigma_nu_sqrtr(ir,r,tprsigma,temp,nu,0)
          if(temp.le.0.d0) then 
              write(*,*) 'ERROR: Temperature <=0 found!'
              stop
          endif
          disk_temp(ir)     = temp
          disk_nu(ir)       = nu
          disk_teff(ir)     = st_tempeff
c
c         **** MAJOR BUG FIX: ORDER OF NEXT FEW LINES (23.05.05) ****
c
          global_alpnew(ir) = st_alpha
          if(time_time.gt.5.53d12.and.time_time.lt.5.55d12) then
              pd_r              = deriv_signusr_r(ir,r,tprsigma,1)
              pd_sigr           = deriv_signusr_sigr(ir,r,tprsigma,1)
          else
              pd_r              = deriv_signusr_r(ir,r,tprsigma,0)
              pd_sigr           = deriv_signusr_sigr(ir,r,tprsigma,0)
          endif
          disk_v0(ir)       = -3.d0*pd_r*nu/r
          disk_v0ori(ir)    = disk_v0(ir)
c
c     Note by Ruud: At some point, I had the feeling that pd_sigr at any
c     radius could never be less than pd_sigr at the inner grid point. I
c     do not remember why I thought so or what led me to that thought.
c     Also, I only have a vague understanding of what this routine does
c     and of the meaning of pd_sigr, so I can't justify this feeling
c     either. This if..endif loop was supposed to fix this, but tests
c     now don't show any problem when some pd_sigr(ir) are less than
c     pd_sigr(1). This may be because of the extra inward radial 
c     velocity term (see a few lines below), which was not yet there
c     when I first thought to see a problem with pd_sigr. Anyway, I'll
c     just disable this "fix" for now.
c
c          if(ir.eq.1) then
c              pd_sigr1 = pd_sigr
c          else
c              if(pd_sigr.lt.0.999d0*pd_sigr1) then
c                  write(*,20) 'pd_sigr',time_time/year,r/au,pd_sigr,
c     $                        pd_sigr1
c   20             format(A7,4(1PE13.5))
c                  pd_sigr = pd_sigr1
c              endif
c          endif
c
          if(time_time.lt.tcloudm) then
c
c            If ml_angmom_mode .ne. 0, then we modify the radial velocity
c            of the disk according to some recipe:
c
             if(ml_angmom_mode.eq.-1) then
c
c               If cloud_idistr is 3 and the system is still accreting, add the
c               extra inward radial term discussed in the paper.
c
                if(cloud_idistr.eq.3) then
                   disk_v0(ir) = disk_v0(ir) - 
     %                           eta_r*sqrt(GG*global_mstar/r)
                endif
             elseif(ml_angmom_mode.eq.1) then
c
c               Do the radial velocity modification using angular 
c               momentum conservation considerations
c
c                 if(ir.gt.iredge.and.r.lt.ml_rcentr) then
c                     nfix = nfix + 1
c                     v0corr(ir) = 0d0
c                 else
                     v0corr(ir) = fourpi * ml_psi(ir) * r*r *
     %                            (1.d0-ml_subkep(ir)) /
     %                            (disk_tprsigma(ir) +
     %                             2d0*pi*r*ml_psi(ir)*time_dt)
c                 endif
             endif
          endif
c
c         The disk diffusion coefficient
c
          disk_dcoef(ir)    = 3.d0*pd_sigr*nu
      enddo
c
c     For ml_angmom_mode=1, do the velocity correction.
c
      if(ml_angmom_mode.eq.1.and.ml_mdotdisk.gt.0d0.and.
     %                                        time_time.lt.tcloudm) then
c          if(nfix.gt.0) then
c              do ir=iredge,grid_nr
c                  if(grid_r(ir).lt.ml_rcentr) then
c                      v0corr(ir) = v0corr(iredge-1) *
c     %                             (1d0 - (1+ir-iredge)/(nfix+2d0))
c                  endif
c              enddo
c          endif
          do ir=1,grid_nr
              disk_v0(ir) = disk_v0(ir) - v0corr(ir)
          enddo
      endif
c
c     Compute the intervace values of of disk_v0
c
c     **** NOTE: **** 27-08-04   [20-02-05: I put this note here
c                                 later. Hope I interpreted it OK.]
c          I changed the indexing of vint, so that it is 
c          consistent with the indexing of 'flux' in the
c          aposteriori_flux_and_source() routine!!!
c
      do ir=2,grid_nr
          disk_iv0(ir) = 0.5d0 * ( disk_v0(ir) + disk_v0(ir-1) )
      enddo
      disk_iv0(1) = 0.d0
c
c     Done
c
      end




c     --------------------------------------------------------------
c          THE CRANK-NICHOLSON SCHEME FOR DIFFUSION-ADVECTION
c
c     Perform one time step for the following PDE:
c
c        du    d /     \    d /           d  / u  \ \
c        -- + -- | u v | - -- | h(x) D(x) -- |----| | = K + L u
c        dt   dx \     /   dx \           dx \h(x)/ /
c
c     with boundary conditions
c
c         du/h |            |
c       p ---- |      + q u |       = r
c          dx  |x=xbc       |x=xbc
c
c     ARGUMENTS:
c       n             Nr of grid points
c       x             The grid
c       u             The current values of u(x)
c       vint          The values for v(x) [NOTE: vint defined at 
c                                                cell interface!]
c       d             The values for D(x)
c       h             The values for h(x)
c       k             The values for K(x)
c       l             The values for L(x)
c       dt            The time step
c       pl,pr         The p value for the lhs/rhs of the BC equation
c       ql,qr         The q value for the lhs/rhs of the BC equation
c       rl,rr         The r value for the lhs/rhs of the BC equation
c       a,b,c,q       Dummy arrays used by the cranknich() routine
c
c     NOTE: This is the new version, in which the bugs are fixed,
c           and in which a K and L term are added on the right-
c           hand-side (01-11-03).
c
c     NOTE: An even newer version (changed name to cranknich_general)
c           which also allows `relative diffusion': the h(x) factor
c           in the equation (27-01-04).
c
c     *** IMPORTANT NOTE: ***
c           This is an adapted version of the cranknich_general
c           routine in which the velocity v is given at the cell
c           boundaries!!! Specially designed for the radial mixing
c           problem... 
c           **** NOTE: **** 27-08-04
c                I changed the indexing of vint, so that it is 
c                consistent with the indexing of 'flux' in the
c                aposteriori_flux_and_source() routine!!!
c
c     --------------------------------------------------------------
      subroutine cranknich_general(n,x,u,vint,d,h,k,l,dt,
     %                     pl,pr,ql,qr,rl,rr,a,b,c,q,store)
      implicit none
      integer n,store
      doubleprecision pl,pr,ql,qr,rl,rr,dt
      doubleprecision x(n),u(n),vint(n+1),d(n),a(n),b(n),c(n),q(n)
      doubleprecision k(n),l(n),h(n)
c
#include "common_disk.h"
c
      integer i
      doubleprecision eps,eps1,dum,dpl,dmn,vpl,vmn
      parameter(eps=THEEPS) ! eps=0 -> exp, eps=1 -> imp, eps=0.5 -> CN
c
      if(vint(1).ne.0.d0) then
          write(*,*) 'ERROR: This is a special version of '
          write(*,*) '   cranknich_general: v defined at interface!'
          write(*,*) '   Found vint(1).ne.0... Should not be!'
          stop
      endif
c
c     The main body of the equation
c
      eps1=1.d0-eps
      do i=2,n-1
          dum  = 2.0*dt/(x(i+1)-x(i-1))
          dpl  = dum*0.25*(h(i+1)+h(i))*(d(i+1)+d(i))/(x(i+1)-x(i))
          dmn  = dum*0.25*(h(i)+h(i-1))*(d(i)+d(i-1))/(x(i)-x(i-1))
          vpl  = dt*vint(i+1)/(x(i+1)-x(i-1))
          vmn  = dt*vint(i)/(x(i+1)-x(i-1))
          b(i) = 1.d0+eps*((dpl+dmn)/h(i)-vmn+vpl-dt*l(i))
          a(i) = -eps*dmn/h(i-1)-eps*vmn
          c(i) = -eps*dpl/h(i+1)+eps*vpl
          q(i) = u(i)+k(i)*dt
     %            -eps1*((dpl+dmn)/h(i)-vmn+vpl-dt*l(i))*u(i)
     %            +eps1*dmn*u(i-1)/h(i-1)+eps1*dpl*u(i+1)/h(i+1)
     %            +eps1*vmn*u(i-1)-eps1*vpl*u(i+1)
      enddo
c
c     The left (i=1) boundary condition
c
      b(1) = ql-pl/(h(1)*(x(2)-x(1)))
      c(1) = pl/(h(2)*(x(2)-x(1)))
      q(1) = rl
c
c     The right (i=n) boundary condition
c
      b(n) = qr+pr/(h(n)*(x(n)-x(n-1)))
      a(n) = -pr/(h(n-1)*(x(n)-x(n-1)))
      q(n) = rr
c
c     Solve the implicit differencing system
c
      call tridag(a,b,c,q,u,n)
c
      end


c     --------------------------------------------------------------
c            A-POSTERIORI COMPUTATION OF FLUXES AND SOURCES
c
c     This routine computes to machine precision what the fluxes
c     at the cell interfaces were over the last time step and what
c     the source terms were. For the cell interfaces between two
c     cell centers this is straight-forward. For the cell interfaces
c     at the edges this is not straightforward, but the fluxes can
c     be computed indirectly. For this to work one has to specify
c     where the boundaries of the outermost cells are. This can be
c     directly at the location of the grid point (possibility 1) or
c     a linear extrapolation (possibility 2):
c
c          *---*-----*-----*
c          | |    |     |      (possibility 1)
c
c          *---*-----*-----*
c        |   |    |     |      (possibility 2)
c
c     The resulting flux and source are such that one can verify:
c
c                            2*dt     /             \
c       unew_i - uold_i = ----------- | F_i+1 - F_i | + dt * Src_i
c                         x_i+1-x_i-1 \             /
c
c     *** IMPORTANT NOTE: ***
c           This is an adapted version of the cranknich_general
c           routine in which the velocity v is given at the cell
c           boundaries!!! Specially designed for the radial mixing
c           problem... 
c           **** NOTE: **** 27-08-04
c                I changed the indexing of vint, so that it is 
c                consistent with the indexing of 'flux' in the
c                aposteriori_flux_and_source() routine!!!
c
c     --------------------------------------------------------------
      subroutine aposteriori_flux_and_source(n,x,uold,unew,vint,
     %                d,h,k,l,dt,ibnd,flux,src)
      implicit none
c
#include "common_disk.h"
c
      integer n,ibnd
      doubleprecision dt,flux(n+1),src(n+1)
      doubleprecision x(n),uold(n),unew(n)
      doubleprecision vint(n+1),d(n),k(n),l(n),h(n)
c
      integer i
      doubleprecision eps,eps1,flold,flnew,v_av,d_av,uold_av,unew_av,dx
      doubleprecision dxl,dxr,h_av
      parameter(eps=THEEPS) ! eps=0 -> exp, eps=1 -> imp, eps=0.5 -> CN
c
      if(vint(1).ne.0.d0) then
          write(*,*) 'ERROR: This is a special version of '
          write(*,*) '   cranknich_general: v defined at interface!'
          write(*,*) '   Found vint(1).ne.0... Should not be!'
          stop
      endif
c
c     Compute the fluxes for i=2,n
c
      eps1=1.d0-eps
      do i=2,n
          v_av    = vint(i)
          d_av    = 0.5*(d(i-1)+d(i))
          h_av    = 0.5*(h(i-1)+h(i))
          uold_av = 0.5*(uold(i-1)+uold(i))
          unew_av = 0.5*(unew(i-1)+unew(i))
          dx      = x(i)-x(i-1)
          flold   = v_av*uold_av - h_av*d_av*
     %                   (uold(i)/h(i)-uold(i-1)/h(i-1))/dx 
          flnew   = v_av*unew_av - h_av*d_av*
     %                   (unew(i)/h(i)-unew(i-1)/h(i-1))/dx 
          flux(i) = eps1*flold + eps*flnew
          src(i)  = k(i) + l(i) * ( eps1*uold(i) + eps*unew(i) )
      enddo
c
c     For i=1 and i=n+1 (the boundaries) we know what the source
c     should be, and we know what the boundary conditions have 
c     produced (unew(1) and unew(n+1)). So we can compute implicitly
c     what the flux must have been.
c
c     BUGFIX 14-Dec-2007: The arrays k and l only have dimension n, not
c     n+1. Both have all zero values, so src(1) and src(n+1) always
c     become zero.
c
c     First compute the source terms:
c
      src(1)   = k(1) + l(1) *  ( eps1*uold(1) + eps*unew(1) )
c      src(n+1) = k(n+1) + l(n+1) *  ( eps1*uold(n+1) + eps*unew(n+1) )
      src(n+1) = k(n) + l(n) *  ( eps1*uold(n) + eps*unew(n) )
c
c     Now specify the dx of the boundary cells according to possibility
c     1 or 2
c
      if(ibnd.eq.1) then
          dxl = 0.5*(x(2)-x(1))
          dxr = 0.5*(x(n)-x(n-1))
      elseif(ibnd.eq.2) then
          dxl = x(2)-x(1)
          dxr = x(n)-x(n-1)
      else
          write(*,*) 'ERROR in aposteriori_flux_and_source():'
          write(*,*) '   Dont know ibnd=',ibnd
          stop
      endif
c
c     Now find the flux at the left and right boundaries
c
      flux(1) = (dxl/dt)*(unew(1)-uold(1))-dxl*src(1)+flux(2)
      flux(n+1) = -(dxr/dt)*(unew(n)-uold(n))+dxr*src(n)+flux(n)
c
c     Done...
c
      end


c     --------------------------------------------------------------
c                   NUMERICAL RECIPES ROUTINE: TRIDAG
c     --------------------------------------------------------------
      SUBROUTINE tridag(a,b,c,r,u,n)
      INTEGER n
      DOUBLEPRECISION a(n),b(n),c(n),r(n),u(n)
      INTEGER j
      DOUBLEPRECISION bet,gam(FRSIZE_R+2)
      if(b(1).eq.0.) then
          write(*,*) 'tridag: rewrite equations'
          stop
      endif
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) then
            write(*,*) 'tridag failed'
            stop
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END
     
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

c     --------------------------------------------------------------
c                   NUMERICAL RECIPES ROUTINE: HUNT
c     --------------------------------------------------------------
      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER jlo,n
      DOUBLEPRECISION x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..


c     --------------------------------------------------------------
c                        INTEGRATION ROUTINE 
c     --------------------------------------------------------------
      function integrate(n,x,f)
      implicit none
      integer n,i
      doubleprecision x(n),f(n),integrate,int
      int=0.d0
      do i=2,n
          int=int+0.5d0*(f(i)+f(i-1))*(x(i)-x(i-1))
      enddo
      if(x(n).gt.x(1)) then
          integrate = int
      else
          integrate = -int
      endif
      return
      end


c     --------------------------------------------------------------
c                     UPDATE THE SIGMA FROM DIFFUSION
c
c     This is the main routine for the viscous evolution of the 
c     disk. 
c     --------------------------------------------------------------
      subroutine update_sigma_diff()
      implicit none
c
#include "common_natconst.h"
#include "common_disk.h"
c
      doubleprecision a(FRSIZE_R),b(FRSIZE_R)
      doubleprecision c(FRSIZE_R),q(FRSIZE_R)
      doubleprecision null(FRSIZE_R),onee(FRSIZE_R)
      doubleprecision infall(FRSIZE_R)
      doubleprecision src(FRSIZE_R+1)
      doubleprecision backup(FRSIZE_R)
      doubleprecision r,tprsigma,signusr,temp,nu
      doubleprecision sigma_nu_sqrtr,deriv_signusr_sigr
      doubleprecision deriv_signusr_r
      doubleprecision pl,pr,ql,qr,rl,rr,dt
      doubleprecision x1,x2,y1,y2,dum1,dum2
      integer ir,iflag,number_invalid
c
c     Set null to null
c
      do ir=1,grid_nr
          null(ir) = 0.d0
          onee(ir) = 1.d0
      enddo
c
c     Backup current tprsigma
c
      do ir=1,grid_nr
          backup(ir) = disk_tprsigma(ir)
      enddo
c
c     Set boundary conditions at inner edge
c
      pl      = 0.d0
      ql      = 2.d0
      rl      = 2.d0*global_sigma_min*2*pi*grid_r(1)
c
c     Set boundary conditions at outer edge
c
      r       = grid_r(grid_nr)
      tprsigma= disk_tprsigma(grid_nr)
      signusr = sigma_nu_sqrtr(grid_nr,r,tprsigma,temp,nu,0)
      pr      = deriv_signusr_sigr(grid_nr,r,tprsigma,0) * nu
      qr      = deriv_signusr_r(grid_nr,r,tprsigma,0) * nu / r
      rr      = 0.d0
c
c     Include mass loading elements 
c
      if((cloud_mass.gt.0.d0).and.(ml_implicit.eq.1)) then
          do ir=1,grid_nr
              infall(ir) = 2*pi*grid_r(ir)*ml_psi(ir)
          enddo
      else
          do ir=1,grid_nr
              infall(ir) = 0.d0
          enddo
      endif
c
c     Perform Crank-Nicholson step and the sources
c
      call cranknich_general(grid_nr,grid_r,disk_tprsigma,
     %           disk_iv0,disk_dcoef,onee,infall,null,
     %           time_dt,pl,pr,ql,qr,rl,rr,a,b,c,q,1,1)
c
c     First check if the results are reasonable
c
      do ir=1,grid_nr
          if(number_invalid(disk_tprsigma(ir)).ne.0) then
              write(*,*) 'ERROR: Invalid disk_tprsigma()'
              write(*,*) '     ir    = ',ir
              write(*,*) '     nu    = ',nu
              write(*,*) '     alpha = ',global_alpha(ir)
              write(*,*) '     T     = ',disk_temp(ir)
              write(*,*) '     tprs  = ',disk_tprsigma(ir)
              write(*,*) pl,ql,rl,pr,qr,rr
              stop
          endif
      enddo
c
c     Check for positiveness
c
      iflag=0
      do ir=1,grid_nr
          if(disk_tprsigma(ir).lt.2.d0*pi*grid_r(ir)*
     %             global_sigma_min) then
              disk_tprsigma(ir) = 2.d0*pi*grid_r(ir)*
     %             global_sigma_min
              iflag=1
          endif
      enddo
      if(iflag.ne.0) then
c          write(*,*) 'Negative Sigma detected'
      endif
c
c     Now a-posteriori exact computation of mass flux at interfaces
c
      call aposteriori_flux_and_source(grid_nr,grid_r,
     %           backup,disk_tprsigma,disk_iv0,disk_dcoef,
     %           onee,null,null,time_dt,2,disk_mdot,src)
c
c     Compute the exact v_R from flux/(2*pi*r*sigma)
c
      do ir=2,grid_nr
          disk_vr(ir) = disk_mdot(ir) / ( 0.5*(disk_tprsigma(ir)+
     %                                    disk_tprsigma(ir-1)) )
      enddo
c      disk_vr(1)         = disk_mdot(1) / disk_tprsigma(1)
c      disk_vr(grid_nr+1) = disk_mdot(grid_nr+1) / 
c     %                     disk_tprsigma(grid_nr)
      disk_vr(1) = 0.d0
      disk_vr(grid_nr+1) = 0.d0
c
c     Done
c
      end


c     --------------------------------------------------------------
c         ROUTINE FOR RADIAL MIXING OF TRACE SPECIES IN A DISK
c
c     This routine solves one time step of the following equation
c     (solving for Sgi, given Sg)
c
c      d(Sgi*R)   d(Sgi*R*vR)    d /               d / Sgi \ \
c      -------- + ----------- - -- | Sg(R)*R*D(R) -- |-----| | = Src
c         dt           dR       dR \              dR \Sg(R)/ /
c
c     which is the conservation equation of a trace element in an
c     active accretion disk according to the equations of Morfill
c     and Voelk. As boundary conditions we take a zero derivative
c     to the Sgi/Sg ratio.
c
c     ARGUMENTS:
c       nr           Nr of radial grid points of the model
c       r            Radial grid
c       sigmatot     The total surface density Sg(R)
c       sigmatrace   The surface density of the trace element
c       src          Source function for dust (infall from a cloud)
c       vr           The radial velocity of the total material
c                    [NOTE: vr should be given at the cell
c                           interfaces!!! vr(i) is velo between
c                           cell i and i+1...]
c       d            The diffusion constant
c       dt           The time step
c       bci          Boundary condition at inner edge. If 1d99, then
c                    take the normal boundary condition
c
c     RESULT:
c       sigmatrace   The sigmatrace will be updated to the new
c                    values.
c
c     IMPORTANT NOTE:
c       It is of *crucial* importance to have an *exact* expression
c       for vr at the cell boundaries. If the vr you give is only
c       formally OK, but not to machine precision equal to the value
c       Mdot(ir+1/2)/(pi*(Sigma(ir)+Sigma(ir+1))*R) at the cell 
c       interface, then the ratio of sigmatrace/sigmatot would start
c       to vary, even when one would start with constant ratio. This
c       is obviously inconsistent.
c     --------------------------------------------------------------
      subroutine dustradmix(nr,r,sigmatot,sigmatrace,src,vr,d,dt,bci)
      implicit none
c
      integer nr,ir
      doubleprecision r(FRSIZE_R),sigmatot(FRSIZE_R),src(FRSIZE_R)
      doubleprecision sigmatrace(FRSIZE_R),bci
      doubleprecision vr(FRSIZE_R+1),d(FRSIZE_R),dt
      doubleprecision u(FRSIZE_R),ubar(FRSIZE_R)
      doubleprecision k(FRSIZE_R),l(FRSIZE_R)
      doubleprecision a(FRSIZE_R),b(FRSIZE_R)
      doubleprecision c(FRSIZE_R),q(FRSIZE_R)
      doubleprecision pl,pr,ql,qr,rl,rr
c
      if(nr.gt.FRSIZE_R) stop 802
      if(vr(1).ne.0.d0) then
          write(*,*) 'ERROR: vr must be cell interface velocity'
          stop
      endif
c
c     Conversion to u and ubar, and setting of dummy arrays
c
      do ir=1,nr
          u(ir)    = sigmatrace(ir)*r(ir)
          ubar(ir) = sigmatot(ir)*r(ir)
          k(ir)    = src(ir)*r(ir)
          l(ir)    = 0.d0
      enddo
c
c     Boundary conditions
c
      if(bci.ge.1d98) then
          pl = 2.d0
          ql = 0.d0
          rl = 0.d0
      else
          pl = 0.d0
          ql = 2.d0
          rl = 2.d0*bci*r(1)
      endif
      pr = 2.d0
      qr = 0.d0
      rr = 0.d0
c     
c     Now call the generalized cranknich routine
c
      call cranknich_general(nr,r,u,vr,d,ubar,k,l,dt,
     %                       pl,pr,ql,qr,rl,rr,a,b,c,q,0)      
c
c     Put result back into sigmatrace
c
      do ir=1,nr
          sigmatrace(ir) = u(ir)/r(ir)
      enddo
c
c     Done...
c
      end



c     --------------------------------------------------------------
c         ROUTINE FOR RADIAL MIXING OF CRYSTALLINE DUST IN A DISK
c              INCLUDING TIME-DEPENDENT CRYSTALLIZATION
c
c     This routine solves one time step of the following equation
c     (solving for Sgi, given Sg)
c
c      d(Sgi*R)   d(Sgi*R*vR)    d /               d / Sgi \ \
c      -------- + ----------- - -- | Sg(R)*R*D(R) -- |-----| | = 
c         dt           dR       dR \              dR \Sg(R)/ /
c
c        (gtd*Sg-Sgi)*R*rate 
c
c     where rate is the crystallization rate and Src the source
c     function.
c
c     ARGUMENTS:
c       nr           Nr of radial grid points of the model
c       r            Radial grid
c       sigmatot     The total surface density Sg(R)
c       sigmatrace   The surface density of the trace element
c       vr           The radial velocity of the total material
c                    [NOTE: vr should be given at the cell
c                           interfaces!!! vr(i) is velo between
c                           cell i and i+1...]
c       d            The diffusion constant
c       dt           The time step
c       bci          Boundary condition at inner edge. If 1d99, then
c                    take the normal boundary condition
c
c     RESULT:
c       sigmatrace   The sigmatrace will be updated to the new
c                    values.
c
c     IMPORTANT NOTE:
c       It is of *crucial* importance to have an *exact* expression
c       for vr at the cell boundaries. If the vr you give is only
c       formally OK, but not to machine precision equal to the value
c       Mdot(ir+1/2)/(pi*(Sigma(ir)+Sigma(ir+1))*R) at the cell 
c       interface, then the ratio of sigmatrace/sigmatot would start
c       to vary, even when one would start with constant ratio. This
c       is obviously inconsistent.
c     --------------------------------------------------------------
      subroutine crystradmix(nr,r,sigmatot,sigmatrace,vr,d,
     %                       crystrate,gtd,dt,bci)
      implicit none
c
      integer nr,ir
      doubleprecision r(FRSIZE_R),sigmatot(FRSIZE_R)
      doubleprecision sigmatrace(FRSIZE_R),bci,gtd
      doubleprecision vr(FRSIZE_R+1),d(FRSIZE_R),dt
      doubleprecision u(FRSIZE_R),ubar(FRSIZE_R)
      doubleprecision crystrate(FRSIZE_R)
      doubleprecision k(FRSIZE_R),l(FRSIZE_R)
      doubleprecision a(FRSIZE_R),b(FRSIZE_R)
      doubleprecision c(FRSIZE_R),q(FRSIZE_R)
      doubleprecision pl,pr,ql,qr,rl,rr
c
      if(nr.gt.FRSIZE_R) stop 803
      if(vr(1).ne.0.d0) then
          write(*,*) 'ERROR: vr must be cell interface velocity'
          stop
      endif
c
c     Conversion to u and ubar, and setting of dummy arrays
c
      do ir=1,nr
          u(ir)    = sigmatrace(ir)*r(ir)
          ubar(ir) = sigmatot(ir)*r(ir)
          k(ir)    = crystrate(ir)*sigmatot(ir)*r(ir)/gtd
          l(ir)    = -crystrate(ir)
      enddo
c
c     Boundary conditions
c
      if(bci.ge.1d98) then
          pl = 2.d0
          ql = 0.d0
          rl = 0.d0
      else
          pl = 0.d0
          ql = 2.d0
          rl = 2.d0*bci*r(1)
      endif
      pr = 2.d0
      qr = 0.d0
      rr = 0.d0
c     
c     Now call the generalized cranknich routine
c
      call cranknich_general(nr,r,u,vr,d,ubar,k,l,dt,
     %                       pl,pr,ql,qr,rl,rr,a,b,c,q,0)      
c
c     Put result back into sigmatrace
c
      do ir=1,nr
          sigmatrace(ir) = u(ir)/r(ir)
      enddo
c
c     Done...
c
      end



c     --------------------------------------------------------------
c                     DO RADIAL MIXING THE OLD WAY 
c     --------------------------------------------------------------
      subroutine dustradmix_old(nr,r,sigmatrace,src,iv0,d,dt,bci)
      implicit none
c
      integer nr,ir
      doubleprecision r(FRSIZE_R),sigmatot(FRSIZE_R),src(FRSIZE_R)
      doubleprecision sigmatrace(FRSIZE_R),bci
      doubleprecision iv0(FRSIZE_R+1),d(FRSIZE_R),dt
      doubleprecision u(FRSIZE_R),ubar(FRSIZE_R)
      doubleprecision k(FRSIZE_R),l(FRSIZE_R)
      doubleprecision a(FRSIZE_R),b(FRSIZE_R)
      doubleprecision c(FRSIZE_R),q(FRSIZE_R)
      doubleprecision pl,pr,ql,qr,rl,rr
c
      if(nr.gt.FRSIZE_R) stop 804
      if(iv0(1).ne.0.d0) then
          write(*,*) 'ERROR: iv0 must be cell interface velocity'
          stop
      endif
c
c     Conversion to u and ubar, and setting of dummy arrays
c
      do ir=1,nr
          u(ir)    = sigmatrace(ir)*r(ir)
          ubar(ir) = 1.d0
          k(ir)    = src(ir)*r(ir)
          l(ir)    = 0.d0
      enddo
c
c     Boundary conditions
c
      if(bci.ge.1d98) then
          pl = 2.d0
          ql = 0.d0
          rl = 0.d0
      else
          pl = 0.d0
          ql = 2.d0
          rl = 2.d0*bci*r(1)
      endif
      pr = 2.d0
      qr = 0.d0
      rr = 0.d0
c     
c     Now call the generalized cranknich routine
c
      call cranknich_general(nr,r,u,iv0,d,ubar,k,l,dt,
     %                       pl,pr,ql,qr,rl,rr,a,b,c,q,0)      
c
c     Put result back into sigmatrace
c
      do ir=1,nr
          sigmatrace(ir) = u(ir)/r(ir)
      enddo
c
c     Done...
c
      end


c     --------------------------------------------------------------
c               CALCULATE THE INFALL ONTO THE DISK SURFACE
c
c     Use Shu model coupled to Ulrich model to calculate the 
c     accretion rate onto the disk surface.
c     --------------------------------------------------------------
      subroutine shuulrich_infall(nr,r,mcloud,cs,omega,time,
     %                            psi,subkep,mdotcap,mdotdisk,mdotflow,
     %                            rcentr,mcentr,ismooth,smoothparam1,
     %                            smoothparam2,idistr)
      implicit none
      integer nr,nrcap,icen,ir,number_invalid,ismooth,idistr,iz,jr,lz
      integer docosth0,icen2,i,nroot,irzero,irmax
      doubleprecision r(FRSIZE_R),psi(FRSIZE_R),rcap(FRSIZE_R)
      doubleprecision subkep(FRSIZE_R)
      doubleprecision dum(FRSIZE_R),fact,integrate,dm
      doubleprecision smoothparam1,smoothparam2,check
      doubleprecision mcloud,cs,omega,time,mcentr,psicap
      doubleprecision tcloud,mfree,mcoll,dnedz,dnddz,dz,dt,scale
      doubleprecision rcoll,rcoll0,mdotcap,mdotdisk,rcentr,mmu0
      doubleprecision rhocmueq,vz,zeta,legendre,vx,rsph,theta,outangle
      doubleprecision sinth,costh,rrc,costh0,aout,surfmax
      doubleprecision distort,accden(FRSIZE_R),flim,ztemp,vzacc,ttemp
      doubleprecision accdenf(FRSIZE_R),psif(FRSIZE_R),mdotflow,acub(4)
      doubleprecision xcub(3)
      parameter(nrcap=30)
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
c
c     Set source function to zero.
c
      do ir=1,nr
          accden(ir) = 0d0
          accdenf(ir) = 0d0
          psi(ir) = 0d0
          psif(ir) = 0d0
      enddo
c
c     Calculate the global numbers for the Shu infall model with
c     rotation.
c
      mdot   = m0*cs**3/GG
      rcloud = GG * mcloud / ( 2*cs**2 )
      tcloud = rcloud/cs
      tcloudm= tcloud*(2./m0)
      if(time.lt.tcloudm) then
          mcentr = (m0/2.)*mcloud*(time/tcloud)
          if(time.lt.tcloud) then 
              mcoll  = mcloud*(time/tcloud)
          else
              mcoll  = mcloud
          endif
      else
          mcentr = mcloud
          mcoll  = mcloud
      endif
      mfree  = mcoll-mcentr
      rcoll  = cs*time 
      rcoll0 = rcoll * (m0/2.d0)
      if(time.gt.tcloud) then
          rcoll = rcloud
      endif
      if(time.gt.tcloudm) then
          rcoll0 = rcloud
      endif
      rcentr = omega**2*rcoll0**4/(GG*(mcentr+1d0))
c
c     Now include, if asked for, the smoothing-off of the infall
c
      if(ismooth.eq.0) then 
c
c         Original Shu type model: an abrupt end to the infall
c
          if(time.gt.tcloudm) then
              mdot   = 0.d0
          endif
          if(idistr.eq.3.and.global_mstar+global_mdisk.ge.mcloud) then
              mdot   = 0.d0
          endif
      elseif(ismooth.eq.1) then
c
c         A linearly smoothed-off version
c
          if(smoothparam1.gt.1.d0) stop 7290
          if(smoothparam1.lt.0.d0) stop 7291
          if(time/tcloudm.gt.1.d0-smoothparam1) then
              dm = (time/tcloudm-1.d0+smoothparam1)/(2*smoothparam1)
              if(dm.lt.0.d0) stop 7292
              mdot = (1.d0-dm)*mdot
              if(mdot.lt.0.d0) mdot=0.d0
          endif
      elseif(ismooth.eq.11) then
c
c         A constant tail of infall
c
          if(time.gt.tcloudm) then
              mdot   = smoothparam1*mdot
              rcentr = smoothparam2
          endif
      elseif(ismooth.eq.21) then
c
c         A constant tail of infall
c
          if(time.gt.tcloudm) then
              mdot   = smoothparam1*mdot*(tcloudm/time)
              rcentr = smoothparam2
          endif
      endif
c
c     Compute the accretion rate onto the disk
c
      do ir=1,nr
          subkep(ir) = 1d0
      enddo
      if(mdot.gt.0.d0) then
c
c         Yes, the mass-loading is active at present...
c
          if(time.le.2.01d4*year) then
              mdotcap = mdot
              goto 300
          endif
          icen = 0
          do ir=1,nr
              if(r(ir).lt.rcentr) icen=ir
          enddo
          icen2 = 0
          do ir=1,nr
              if(r(ir).lt.2d0*rcentr) icen2=ir
          enddo
          psi(1) = 0.d0
          mdotcap = 0d0
          mdotdisk = 0d0
          mdotflow = 0d0
          if(idistr.eq.0) then
c
c             Distribute the matter over the disk in the Ulrich way
c
              do ir=1,icen-1
                  mmu0      = sqrt(1-r(ir)/rcentr)
                  rhocmueq  = mdot*(r(ir)/rcentr)/
     %                 ((4*pi*sqrt(GG*mcentr*r(ir)**3))*(2*mmu0**2))
                  vz        = sqrt(GG*mcentr/r(ir))*mmu0
                  psi(ir)   = 2*rhocmueq*vz
                  subkep(ir) = sqrt(1d0-mmu0**2)
              enddo
          else if(idistr.eq.1) then
c
c             Distribute the matter over the disk in such a way that
c             there is no friction of the infalling matter with the disk
c             (see Hueso & Guillot 2005)
c
              do ir=1,icen-1
                  psi(ir)   = (mdot/(8*pi*r(ir)*rcentr)) *
     %                        ((r(ir)/rcentr)*
     %                        (1.d0-sqrt(r(ir)/rcentr)))**(-0.5)
                  subkep(ir) = 1.d0
              enddo
          else if(idistr.eq.2) then
c
c     Distribute the matter over the disk in the Cassen & Moosman (1981)
c     way as described in Terebey, Shu & Cassen (1984).
c
c              do ir=1,icen-1
              do ir=1,icen
                  zeta = rcentr / r(ir)
                  mmu0 = sqrt(1d0-1d0/zeta)
                  legendre = 0.5d0 * (3d0*mmu0**2-1d0)
                  vx = -sqrt(GG*mcentr/r(ir))
                  vz = -mmu0 * vx
                  rhocmueq = -mdot / (4d0 * pi * r(ir)**2 * vx *
     $                                (1d0 + 2d0*zeta*legendre))
                  psi(ir) = 2d0 * rhocmueq * vz
                  subkep(ir) = sqrt(1d0-mmu0**2)
              enddo
          else if(idistr.eq.3) then
c
c     Take the vertical profile of the disk into account for computing
c     the mass loading. First, compute the "surface" of the disk: the
c     point where disk density and envelope density are equal. If no
c     such point exists for a given R, place the surface at 1e90 cm.
c
             iredge = nr+1
             do ir=nr,1,-1
                 disk_mu0(ir) = -1d0
                 if(disk_midprs(ir).le.env_midprs(ir)) then
c                 if(disk_midden(ir).le.env_midden(ir)) then
c
c     If the disk density at the midplane is lower than the envelope
c     density, the disk may actually still exist at this radius. So,
c     scan upwards to see if the disk density exceeds the envelope
c     density at any other z.
c
                     do iz=1,grid_nz
                         if(disk_vertprs(ir,iz).gt.1d-50.and.
     $                       disk_vertprs(ir,iz).gt.env_vertprs(ir,iz))
c                         if(disk_vertden(ir,iz).gt.1d-50.and.
c     $                       disk_vertden(ir,iz).gt.env_vertden(ir,iz))
     $                                                              then
                             jz = iz
                             goto 150
                         endif
                     enddo
c
c     If we get here, it means there is no disk at this radius. Compute
c     the density and vertical velocity at the midplane according to
c     Cassen & Moosman (1981). If we're outside the centrifugal radius,
c     set these to zero.
c
                     disk_surface(ir) = 0d0
                     acc_surface(ir) = 0d0
                     if(r(ir).lt.rcentr) then
                         zeta = rcentr / r(ir)
                         disk_mu0(ir) = sqrt(1d0-1d0/zeta)
                         call get_vertden(r(ir),0d0,r(ir),0.5d0*pi,
     $                                    time,disk_mu0(ir),0)
                         accden(ir) = mden * mpeffinv
                     else
                         accden(ir) = 0d0
                         goto 160
                     endif
                     goto 160
                 endif
                 jz = 1
  150            continue
c
c     Scan upwards from the midplane to find the point where the disk
c     density equals the envelope density. Also, make sure that the disk
c     actually extends to this point.
c
                 do iz=jz,grid_nz
                    if(iz.eq.1) then
                       if(disk_midprs(ir).gt.1d-50.and.
     $                    disk_vertprs(ir,1).le.env_vertprs(ir,1)) then
c                       if(disk_midden(ir).gt.1d-50.and.
c     $                    disk_vertden(ir,1).le.env_vertden(ir,1)) then
c
c     At this point, the envelope density just above the midplane is
c     higher than the disk density, but we still need to check if the
c     disk density exceeds the envelope density at a higher altitude. If
c     so, the disk still exists at this radius.
c
                          do lz=iz,grid_nz
                             if(disk_vertprs(ir,lz).gt.1d-50.and.
     $                                          disk_vertprs(ir,lz).gt.
     $                                          env_vertprs(ir,lz)) then
c                             if(disk_vertden(ir,lz).gt.1d-50.and.
c     $                                          disk_vertden(ir,lz).gt.
c     $                                          env_vertden(ir,lz)) then
                                jz = lz
                                goto 150
                             endif
                          enddo
                          dnedz = (log(env_vertprs(ir,1)) -
     $                             log(env_midprs(ir))) / grid_z(1)
c                          dnedz = (log(env_vertden(ir,1)) -
c     $                             log(env_midden(ir))) / grid_z(1)
                          scale = sqrt(kk*disk_temp(ir)*(grid_r(ir)**3)*
     $                                 mpeffinv*GGinv/disk_mgrav(ir))
                          disk_surface(ir) = scale**2 * (-dnedz +
     $                       sqrt(dnedz**2-2d0*log(env_midprs(ir)/
     $                            disk_midprs(ir))/scale**2))
c                          disk_surface(ir) = scale**2 * (-dnedz +
c     $                       sqrt(dnedz**2-2d0*log(env_midden(ir)/
c     $                            disk_midden(ir))/scale**2))
                          acc_surface(ir) = disk_surface(ir)
c     if pressure:
                          rsph = sqrt(r(ir)**2+disk_surface(ir)**2)
                          theta = atan2(r(ir),disk_surface(ir))
                          call get_vertden(r(ir),disk_surface(ir),rsph,
     $                                     theta,time,-1d0,0)
                          accden(ir) = mden * mpeffinv
c     if density:
c                          accden(ir) = disk_midden(ir) *
c     $                       exp(-0.5d0*(disk_surface(ir)/scale)**2)
c     ------------
                          if(iredge.gt.nr) iredge = ir
                          goto 160
                       endif
                    else
                       if(disk_vertprs(ir,iz-1).gt.1d-50.and.
     $                   disk_vertprs(ir,iz).le.env_vertprs(ir,iz)) then
c                       if(disk_vertden(ir,iz-1).gt.1d-50.and.
c     $                   disk_vertden(ir,iz).le.env_vertden(ir,iz)) then
                          do lz=iz,grid_nz
                             if(disk_vertprs(ir,lz).gt.1d-50.and.
     $                                          disk_vertprs(ir,lz).gt.
     $                                          env_vertprs(ir,lz)) then
c                             if(disk_vertden(ir,lz).gt.1d-50.and.
c     $                                          disk_vertden(ir,lz).gt.
c     $                                          env_vertden(ir,lz)) then
                                jz = lz
                                goto 150
                             endif
                          enddo
                          dnedz = (log(env_vertprs(ir,iz)) -
     $                             log(env_vertprs(ir,iz-1))) /
     $                            (grid_z(iz) - grid_z(iz-1))
c                          dnedz = (log(env_vertden(ir,iz)) -
c     $                             log(env_vertden(ir,iz-1))) /
c     $                            (grid_z(iz) - grid_z(iz-1))
                          scale = sqrt(kk*disk_temp(ir)*(grid_r(ir)**3)*
     $                                 mpeffinv*GGinv/disk_mgrav(ir))
                          disk_surface(ir) = scale**2 * (-dnedz +
     $sqrt(dnedz**2-2d0*(log(env_vertprs(ir,iz-1)/disk_midprs(ir))-
     $dnedz*grid_z(iz-1))/scale**2))
c                          disk_surface(ir) = scale**2 * (-dnedz +
c     $sqrt(dnedz**2-2d0*(log(env_vertden(ir,iz-1)/disk_midden(ir))-
c     $dnedz*grid_z(iz-1))/scale**2))
                          acc_surface(ir) = disk_surface(ir)
c     if pressure:
                          rsph = sqrt(r(ir)**2+disk_surface(ir)**2)
                          theta = atan2(r(ir),disk_surface(ir))
                          call get_vertden(r(ir),disk_surface(ir),rsph,
     $                                     theta,time,-1d0,0)
                          accden(ir) = mden * mpeffinv
c     if density:
c                          accden(ir) = disk_midden(ir) *
c     $                       exp(-0.5d0*(disk_surface(ir)/scale)**2)
c     ------------
                          if(iredge.gt.nr) iredge = ir
                          goto 160
                       endif
                    endif
                 enddo
  160            continue
             enddo
c
c     Find the radius where the surface reaches its maximum.
c
             irmax = 0
             surfmax = 0d0
             irzero = nr+1
             jr = 1
  191        continue
             do ir=jr,nr
                 if(disk_surface(ir).gt.surfmax) then
                     irmax = ir
                     surfmax = disk_surface(ir)
                 endif
             enddo
             if(irmax.eq.0) goto 190
             if(irmax.eq.nr) stop 19011
c
c     Check if there is a point outside R(irmax) where the surface goes
c     up again.
c
             do ir=irmax+1,nr
                 if(disk_surface(ir).gt.disk_surface(ir-1)) then
                     if(disk_surface(ir).lt.
     $                                   0.5d0*disk_surface(irmax)) then
                         irzero = ir
                         goto 190
                     else
                         surfmax = 0d0
                         irzero = nr+1
                         jr = ir
                         goto 191
                     endif
                 endif
             enddo
  190        continue
c
c     Check if material can indeed reach each point on the surface.
c
             do ir=nr,1,-1
                 if(accden(ir).eq.0d0) goto 100
                 if(ir.ge.irzero) then
                     acc_surface(ir) = 0d0
                     if(r(ir).lt.rcentr) then
                         zeta = rcentr / r(ir)
                         disk_mu0(ir) = sqrt(1d0-1d0/zeta)
                         call get_vertden(r(ir),0d0,r(ir),0.5d0*pi,
     $                                    time,disk_mu0(ir),0)
                         accden(ir) = mden * mpeffinv
                     else
                         accden(ir) = 0d0
                     endif
                     goto 100
                 endif
  232            continue
                 docosth0 = 1
                 costh0 = -1d0
                 jr = ir
                 ztemp = acc_surface(ir)
                 ttemp = time
  230            continue
                 if(jr.ge.irzero) goto 100
                 rsph = sqrt(r(jr)**2+ztemp**2)
                 theta = atan2(r(jr),ztemp)
                 costh = ztemp / rsph
                 sinth = r(jr) / rsph
                 call get_vertden(r(jr),ztemp,rsph,theta,time,costh0,0)
                 if(docosth0.eq.1) then
                     costh0 = global_mu0
                     docosth0 = 0
                 endif
                 vx = ux
                 vz = uz
c
c     If the velocities are zero, there is no accretion, so skip ahead.
c
                 if(jr.eq.ir.and.vx.eq.0d0.and.vz.eq.0d0) goto 231
c
c     Now take small steps back along the infall trajectory to see if it
c     intersects the disk anywhere. If so, this trajectory cannot reach
c     the point where we started, so we have to raise the disk surface
c     there. If the trajectory originated from inside the outflow
c     region, set the accretion density to zero and skip ahead.
c
                 if(jr.lt.nr) then
                     dt = abs((r(jr+1)-r(jr)) / vx)
                     ztemp = ztemp - vz * dt
                     ttemp = ttemp - dt
                     if(ttemp.le.0d0) goto 100
                     if(number_invalid(vx).ne.0.or.
     $                                     number_invalid(vz).ne.0) then
                         write(*,*) time/year,r(jr)/au,ztemp/au,
     $                              vx,vz,dt
                         write(*,*) ir,disk_surface(ir)/au,
     $rsph/(cloud_asound*ttemp),(cloud_omega*ttemp)**2
                         stop 229
                     endif
                     outangle = 0.5519d0*(ttemp/tcloudm)**0.26d0
                     if(outangle.gt.pi/2) outangle = pi/2
                     aout = 0.191d0 / (ttemp/tcloudm)**3
                     if(aout.lt.(0.191d0)) aout = 0.191d0
                     if(ztemp.gt.aout*au*(r(jr+1)/au)**1.5d0.and.
     $                                       time.gt.0.2d0*tcloudm) then
                         accdenf(ir) = accden(ir)
                         accden(ir) = 0d0
                         goto 100
                     endif
                     if(ztemp.lt.disk_surface(jr+1)) then
                         accden(ir) = 0d0
                         if(acc_surface(ir).gt.0d0) then
                             acc_surface(ir) = acc_surface(ir) +
     $                                        (disk_surface(jr+1)-ztemp)
                             goto 232
                         endif
                         goto 100
                     endif
c
c     Take another step if we haven't reached the edge of the grid.
c
                     if(jr.lt.nr-1) then
                         jr = jr + 1
                         goto 230
                     endif
                 endif
                 goto 100
  231            continue
                 accden(ir) = 0d0
  100            continue
             enddo
c
c     Now that we know the disk surface at each radius, we can compute
c     the vertical velocity and the mass loading (psi) at those points.
c
             do ir=1,nr
                if(acc_surface(ir).gt.1d80) then
                   psi(ir) = 0d0
                   subkep(ir) = 1d0
                   goto 110
                endif
c
c     If the disk surface falls below the first vertical grid point, the
c     disk height is effectively zero and we use the velocities at the
c     midplane.
c
                if(acc_surface(ir).eq.0d0) then
                   if(r(ir).ge.rcentr) then
                      psi(ir) = 0d0
                      subkep(ir) = 1d0
                      goto 110
                   endif
                   zeta = rcentr / r(ir)
                   disk_mu0(ir) = sqrt(1d0-1d0/zeta)
                   call get_vertden(r(ir),0d0,r(ir),0.5d0*pi,
     $                              time,disk_mu0(ir),0)
                   vzacc = uz
c                   if(time.gt.1.44d12.and.time.lt.1.45d12.and.
c     $                                     ir.gt.770.and.ir.lt.800) then
c                       write(*,*) 2,ir,vzacc,0d0
c                   endif
                   subkep(ir) = uphi / sqrt(GG*global_mstar/r(ir))
c                   subkep(ir) = sqrt(1d0-disk_mu0(ir)**2)
                   goto 111
                endif
c
c     Now do normal accretion.
c
                ztemp = acc_surface(ir)
                ttemp = time
                rsph = sqrt(r(ir)**2+ztemp**2)
                theta = atan2(r(ir),ztemp)
                costh = ztemp / rsph
                sinth = r(ir) / rsph
                call get_vertden(r(ir),ztemp,rsph,theta,time,-1d0,0)
                disk_mu0(ir) = global_mu0
                vzacc = uz
c                if(time.gt.1.44d12.and.time.lt.1.45d12.and.
c     $                                     ir.gt.770.and.ir.lt.800) then
c                    call get_vertden(r(ir),ztemp,rsph,theta,time,-1d0,1)
c                    disk_mu0(ir) = global_mu0
c                    vzacc = uz
c                    write(*,*) 3,ir,vzacc,ztemp/au
c                else
c                    call get_vertden(r(ir),ztemp,rsph,theta,time,-1d0,0)
c                    disk_mu0(ir) = global_mu0
c                    vzacc = uz
c                endif
                subkep(ir) = uphi / sqrt(GG*global_mstar/r(ir))
  111           continue
                psi(ir) = abs(2d0 * accden(ir) * mp * global_mugas *
     $                        vzacc)
                psif(ir) = abs(2d0 * accdenf(ir) * mp * global_mugas *
     $                        vzacc)
c                if(time.gt.1.44d12.and.time.lt.1.45d12.and.
c     $                                     ir.gt.770.and.ir.lt.800) then
c                    write(*,*) 1,ir,psi(ir),accden(ir),vzacc
c                endif
  110           continue
                if(psi(ir).eq.0d0) subkep(ir) = 1d0
             enddo
             psi(1) = 0d0
             psif(1) = 0d0
             subkep(1) = 1d0
c
c     Adjust the accretion rate for the area between the outer edge of
c     the disk and the centrifugal radius.
c
c             if(time.gt.3d4*year) then
c                 if(psi(iredge).gt.0d0) then
c                     do ir=iredge+1,nr
c                         if(r(ir).gt.rcentr) goto 120
c                         psi(ir) = psi(iredge)
c                     enddo
c                 endif
c             endif
c  120        continue
          else
              write(*,*) 'Distribution mode ',idistr,'not supported.'
              stop
          endif
          do ir=1,nr
              if(psi(ir).eq.0d0.and.subkep(ir).ne.1d0) then
                  write(*,*) time/year,ir,psi(ir),subkep(ir)
                  stop 9922
              endif
          enddo
c
c     For idistr.ne.3, set accretion outside the centrifugal radius to
c     zero. For idistr.eq.3, accretion can take place on the part of the
c     disk extending beyond the centrifugal radius.
c
          if(idistr.ne.3) then
             do ir=max(icen,1),nr
                psi(ir) = 0.d0
             enddo
          endif
c
c         Compute rate of accretion inward of inner R coordinate
c         This mass is assumed to be directly loaded onto the star
c
          if(rcentr.le.r(1)) then
c
c             All matter falls within R_in, which is by definition the
c             capture radius of the central star.
c
              if(icen.ne.0) stop 19009
              mdotcap = mdot
          else
              if(nrcap.gt.FRSIZE_R) stop 71200
              if(idistr.eq.0) then
c
c                 Ulrich way
c
                  do ir=1,nrcap
                      rcap(ir)   = r(1)*ir/(1.d0*nrcap)
                      mmu0       = sqrt(1-rcap(ir)/rcentr)
                      rhocmueq   = mdot*(rcap(ir)/rcentr)/
     %                     ((4*pi*sqrt(GG*mcentr*rcap(ir)**3))*
     %                     (2*mmu0**2))
                      vz         = sqrt(GG*mcentr/rcap(ir))*mmu0
                      psicap     = 2*rhocmueq*vz
                      dum(ir)    = 2*pi*rcap(ir)*psicap
                  enddo
              else if(idistr.eq.1) then
c
c                 Hueso & Guillot way
c
                  do ir=1,nrcap
                      rcap(ir)   = r(1)*ir/(1.d0*nrcap)
                      psicap     = (mdot/(8*pi*rcap(ir)*rcentr)) *
     %                        ((rcap(ir)/rcentr)*
     %                        (1.d0-sqrt(rcap(ir)/rcentr)))**(-0.5)
                      dum(ir)    = 2*pi*rcap(ir)*psicap
                  enddo
              else if(idistr.eq.2) then
c
c                 Cassen & Moosman (1981) way.
c
                  do ir=1,nrcap
                      rcap(ir) = r(1)*ir/(1.d0*nrcap)
                      zeta = rcentr / rcap(ir)
                      mmu0 = sqrt(1d0-1d0/zeta)
                      legendre = 0.5d0 * (3d0*mmu0**2-1d0)
                      vx = -sqrt(GG*mcentr/rcap(ir))
                      vz = -mmu0 * vx
                      rhocmueq = -mdot / (4d0 * pi * rcap(ir)**2 * vx *
     $                                    (1d0 + 2d0*zeta*legendre))
                      psicap = 2d0 * rhocmueq * vz
                      dum(ir) = 2d0 * pi * rcap(ir) * psicap
                  enddo
              else if(idistr.eq.3) then
c
c                 Visser et al. (2009) way.
c
                  do ir=1,nrcap
                      rcap(ir) = r(1)*ir/(1.d0*nrcap)
                      zeta = rcentr / rcap(ir)
                      mmu0 = sqrt(1d0-1d0/zeta)
                      call get_vertden(rcap(ir),0d0,rcap(ir),0.5d0*pi,
     $                                 time,mmu0,0)
                      psicap = abs(2d0 * mden * uz)
                      dum(ir) = 2d0 * pi * rcap(ir) * psicap
                  enddo
              else
                  write(*,*) 'Distribution mode ',idistr,
     $                       'not supported.'
                  stop
              endif
              mdotcap = integrate(nrcap,rcap,dum)
              if(mdotcap.gt.2.*mdot) stop 90909
              if(mdotcap.lt.0.d0) then 
                  write(*,*) (rcap(ir),ir=1,nrcap)
                  write(*,*) (dum(ir),ir=1,nrcap)
                  stop 811
              endif
          endif
c
c         Normalize this to make sure the mass loading is perfect
c
          do ir=1,nr
             dum(ir) = 2*pi*r(ir)*psi(ir)
          enddo
          mdotdisk = integrate(nr,r,dum)
          do ir=1,nr
             dum(ir) = 2*pi*r(ir)*psif(ir)
          enddo
          mdotflow = integrate(nr,r,dum)
          if(mdotflow.lt.0d0) mdotflow = 0d0
          fact = mdot / ( mdotcap + mdotdisk + mdotflow)
          if(number_invalid(fact).ne.0) then
              write(*,*) 'fact = ',fact
              stop
          endif
          flim = 1.5d0
          if(idistr.ne.3.and.((fact.gt.flim).or.(fact.lt.1d0/flim)))
     $                                                              then
              write(*,*) 'ERROR: Something wrong with normalization'
              write(*,*) '       of the mass loading...'
              write(*,*) '       Factor wrong is: ',fact
              write(*,*) mdot,mdotcap,mdotdisk,rcentr/au
              stop
          endif
          do ir=1,nr
              psi(ir) = psi(ir) * fact
          enddo
          mdotcap = mdotcap * fact
          mdotdisk = mdotdisk * fact
          mdotflow = mdotflow * fact
  300     continue
      else
c
c         Apparently there is no mass loading (anymore)
c
          mdotcap = 0.d0
          do ir=1,nr
              psi(ir) = 0.d0
          enddo
          mdotdisk = 0.d0
          mdotflow = 0.d0
      endif
c
c     Done...
c
      end


c     --------------------------------------------------------------
c                 SAVE THE DATA BY APPENDING TO FILE
c     --------------------------------------------------------------
      subroutine save_data
      implicit none
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
c
      integer ir,iz,ith,idust
      doubleprecision rcav1,rcav2,rcav
c
      open(unit=1,file='sigma.dat',status='unknown',access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,10) disk_tprsigma(ir)/
     %                (2*pi*grid_r(ir))+1.d-98
 10       format(1PE13.6)
 11       format(2(1PE13.6,1X))
 12       format(3(1PE13.6,1X))
      enddo
      if(global_ndust.gt.0) then
          do idust=1,global_ndust
              do ir=1,grid_nr
                  write(1,10) disk_tprsigmadust(ir,idust)/
     %                        (2*pi*grid_r(ir))+1.d-98
              enddo
          enddo
      endif
      close(1)
c
      if(dovert.eq.1) then
          open(unit=1,file='midden.dat',status='unknown',
     %         access='append')
          write(1,*) 
          do ir=1,grid_nr
              write(1,11) disk_midden(ir),env_midden(ir)
          enddo
          close(1)
c
          open(unit=1,file='midprs.dat',status='unknown',
     %         access='append')
          write(1,*) 
          do ir=1,grid_nr
              write(1,11) disk_midprs(ir),env_midprs(ir)
          enddo
          close(1)
c
          if(write_vertden.eq.1) then
              open(unit=1,file='vertden.dat',status='unknown',
     $             access='append')
              write(1,*) 
              do ir=1,grid_nr
                  do iz=1,grid_nz
                      write(1,11) disk_vertden(ir,iz)+1d-98,
     $                            env_vertdenf(ir,iz)+1d-98
                  enddo
              enddo
              close(1)
c
c   Added by Daniel
c   Write the velocities
c

              open(unit=1,file='env_vr.dat',status='unknown',
     $            access='append')
              write(1,*)
              do ir=1,grid_nr
                  do iz=1,grid_nz
c                      write(*,*) env_ur(ir,iz)
                      write(1,10) env_ur(ir,iz)
                  enddo
              enddo
              close(1)

              open(unit=1,file='env_vth.dat', status='unknown',
     $             access='append')
              write(1,*)
              do ir=1,grid_nr
                  do iz=1,grid_nz
c                     write(*,*) env_uth(ir,iz)
                      write(1,10) env_uth(ir,iz)
                  enddo
              enddo
              close(1)

              open(unit=1,file='env_vphi.dat', status='unknown',
     $             access='append')
              write(1,*)
              do ir=1,grid_nr
                  do iz=1,grid_nz
c                        write(*,*) env_uphi(ir,iz)
                        write(1,10) env_uphi(ir,iz)
                  enddo
              enddo
              close(1)

c
c              open(unit=1,file='vertprs.dat',status='unknown',
c     $             access='append')
c              write(1,*) 
c              do ir=1,grid_nr
c                  do iz=1,grid_nz
c                      write(1,11) disk_vertprs(ir,iz)+1d-98,
c     $                            env_vertprs(ir,iz)+1d-98
c                  enddo
c              enddo
c              close(1)
          endif
      endif
c
      open(unit=1,file='velo.dat',status='unknown',access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,10) disk_vr(ir)
      enddo
!       if(global_ndust.gt.0) then
!           do idust=1,global_ndust
!               do ir=1,grid_nr
!                   write(1,10) disk_vrdust(ir,idust)
!               enddo
!           enddo
!       endif
      close(1)
c
      open(unit=1,file='temperature.dat',status='unknown',
     %           access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,10) disk_temp(ir)
      enddo
      close(1)
c
      open(unit=1,file='dcoef.dat',status='unknown',
     %           access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,10) disk_dcoef(ir)
      enddo
      close(1)
c
      open(unit=1,file='v0.dat',status='unknown',
     %           access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,11) disk_v0(ir),disk_v0ori(ir)
      enddo
      close(1)
c
      open(unit=1,file='mgrav.dat',status='unknown',
     %           access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,10) disk_mgrav(ir)
      enddo
      close(1)
c
      open(unit=1,file='visc.dat',status='unknown',
     %           access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,10) disk_nu(ir)
      enddo
      close(1)
c
      open(unit=1,file='alphagrav.dat',status='unknown',
     %           access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,10) global_alpnew(ir)
      enddo
      close(1)
c
      open(unit=1,file='mdot.dat',status='unknown',
     %           access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,10) disk_mdot(ir)
      enddo
      close(1)
c
      open(unit=1,file='mstar.dat',status='unknown',
     %           access='append')
      write(1,10) global_mstar
      close(1)
c
      open(unit=1,file='mdisk.dat',status='unknown',
     %           access='append')
      write(1,10) global_mdisk
      close(1)
c
      open(unit=1,file='msyst.dat',status='unknown',
     %           access='append')
      write(1,10) global_msyst
      close(1)
c
      open(unit=1,file='lstar.dat',status='unknown',
     %           access='append')
      write(1,10) global_lstar
      close(1)
c
      open(unit=1,file='rstar.dat',status='unknown',
     %           access='append')
      write(1,10) global_rstar
      close(1)
c
      open(unit=1,file='tstar.dat',status='unknown',
     %           access='append')
      write(1,10) global_tstar
      close(1)
c
      open(unit=1,file='surface.dat',status='unknown',
     %           access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,11) disk_surface(ir),acc_surface(ir)
      enddo
      close(1)
c
      open(unit=1,file='psi.dat',status='unknown',
     %           access='append')
      do ir=1,grid_nr
          write(1,10) ml_psi(ir)
      enddo
      close(1)
c
      open(unit=1,file='subkep.dat',status='unknown',
     %           access='append')
      do ir=1,grid_nr
          write(1,10) ml_subkep(ir)
      enddo
      close(1)
c
      open(unit=1,file='mu0.dat',status='unknown',
     %           access='append')
      do ir=1,grid_nr
          write(1,10) disk_mu0(ir)
      enddo
      close(1)
c
      open(unit=1,file='j_outflow.dat',status='unknown',
     %           access='append')
      write(1,10) global_j_in,global_j_out
      close(1)
c
      open(unit=1,file='infall.dat',status='unknown',access='append')
      write(1,12) ml_mdotcap,ml_rcentr,ml_mdotdisk
      close(1)
c
      open(unit=1,file='infall2.dat',status='unknown',access='append')
      write(1,11) ml_mdot
      close(1)
c
      open(unit=1,file='time.dat',status='unknown',access='append')
      write(1,10) time_time
      close(1)
c
      time_nr_save = time_nr_save + 1
      open(unit=1,file='time.info',status='unknown')
      write(1,*) time_nr_save
      close(1)
c
c     Done
c
      end



c     --------------------------------------------------------------
c                             DO THE RUN
c     --------------------------------------------------------------
      subroutine do_run
      implicit none
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
c
      integer it,itsave,ir,idust,irs,iz,ith,j
      doubleprecision timesave,dummy,xp,get_lstar,scale,rhom
      doubleprecision sigma(FRSIZE_R),diffcoef(FRSIZE_R)
      doubleprecision sigmadust(FRSIZE_R),vr(FRSIZE_R)
      doubleprecision src(FRSIZE_R),crystrate(FRSIZE_R)
      doubleprecision xpos,zpos,eps,epsn,temptemp,tempmgrav,temptprsigma
      doubleprecision rrc,costh,costh0,xx,rsph,theta
      doubleprecision legendre,check,outangle,sinth,aout,mden81,mden84
      doubleprecision prs81,prs84,tau,tau2,ratio,maxratio
      logical simstop,dump,ttauri
c
c     Reset time stuff
c
      itsave   = 1
      timesave = 0.d0
      global_msyst = cloud_mass
c
c     Set the initial 2-D density to a singular isothermal sphere.
c
      if(dovert.eq.1) then
          do ir=1,grid_nr
              do iz=1,grid_nz
                  call get_tscvelden(sqrt(grid_r(ir)**2+grid_z(iz)**2),
     $                               atan2(grid_r(ir),grid_z(iz)),
     $                               0d0,1d0,0)
                  env_vertden(ir,iz) = max(mden * mpeffinv,1d-50)
                  disk_vertden(ir,iz) = 1d-50
                  env_vertprs(ir,iz) = env_vertden(ir,iz) * mp *
     $                                 global_mugas * cloud_asound**2
                  disk_vertprs(ir,iz) = 1d-50
                  disk_ur(ir,iz) = 0d0
                  disk_uth(ir,iz) = 0d0
                  disk_uphi(ir,iz) = 0d0
c                 write(*,*) ur,uth,uphi
                  env_ur(ir,iz) = ur
                  env_uth(ir,iz) = uth
                  env_uphi(ir,iz) = uphi
              enddo
          enddo
      endif
c
c     First save data
c
      call save_data()
c
c     Then do the time loop
c
      simstop = .false.
      dump    = .false.
      ttauri  = .false.
      do it=1,time_maxstep
c
c         Reset flags
c
          global_flag_lowtemp = 0
c
c         Set the time step
c
          time_dt = max(time_dtmin,time_time*time_dtrel)
c  Testing for parcels: 11/10/2011
c          time_dt = 0.25*time_dt
c          write(*,*) time_dt
c
c         Check if we are about to cross the save time
c
          if(timesave+time_dt.ge.time_tsave) then
              time_dt = time_tsave - timesave
              dump = .true.
          endif
c
c         Check if we are about to cross the time-limit
c
          if(time_time+time_dt.gt.time_end) then
              time_dt = time_end - time_time
              simstop = .true.
          endif
c
c         Message
c
c          write(*,*) 'Time step ',it,', Time = ',time_time/year,' year'
c
c         Prepare mass loading sources
c
          if(cloud_mass.gt.0.d0) then
c
c             Check if cloud parameters are consistent 
c
              if((cloud_asound.eq.0.d0).or.(cloud_omega.eq.0.d0)) then
                  write(*,*) 'ERROR: Input parameters of cloud contain'
                  write(*,*) '       values of 0.'
                  stop
              endif
c
c             Compute the rotating Shu/Ulrich infall model, and compute
c             the mass accretion rate onto the surface of the disk
c             (psi(R)) and the mass accretion rate inward of the capture
c             radius (the `effective radius' of the central star).
c
              call shuulrich_infall(grid_nr,grid_r,cloud_mass,
     %                      cloud_asound,cloud_omega,time_time,
     %                      ml_psi,ml_subkep,ml_mdotcap,ml_mdotdisk,
     %                      ml_mdotflow,ml_rcentr,ml_mcentr,
     %                      cloud_ismooth,cloud_smoothparam1,
     %                      cloud_smoothparam2,cloud_idistr)
              ml_mdot = mdot
          endif
c
c     Check if we're not dumping too much mass in too short a time.
c
c          if(iredge.ge.grid_nr) goto 710
c          if(grid_r(iredge).ge.ml_rcentr) goto 710
c          maxratio = 1d-90
c          do ir=iredge+1,grid_nr
c              if(grid_r(ir).gt.ml_rcentr) goto 700
c              ratio = 2d0*pi*grid_r(ir)*ml_psi(ir)*time_dt/
c     $                disk_tprsigma(ir)
c              if(ratio.gt.maxratio) maxratio = ratio
c          enddo
c  700     continue
c          if(maxratio.gt.1d-1) then
c              time_dt = time_dt / (1d1*maxratio)
c          endif
c  710     continue
c
c     If we just crossed into the T Tauri phase, store the current
c     star and disk mass and inform the user.
c
          if((time_time.gt.tcloudm.or.
     $            global_mstar+global_mdisk.ge.cloud_mass.or.
     $            global_mstar+global_mdisk.ge.global_msyst).and.
     $            .not.ttauri) then
              mstar_ttauri = global_mstar
              mdisk_ttauri = global_mdisk
              write(*,*) 'Masses at start of T Tauri phase:'
              write(*,*) '    star = ',mstar_ttauri/MS,'MS'
              write(*,*) '    disk = ',mdisk_ttauri/MS,'MS'
              write(*,*) 'envelope = ',(cloud_mass-mstar_ttauri-
     $                                  mdisk_ttauri)/MS,'MS'
              ttauri = .true.
          endif

c
c         First (re)compute the gravitational mass
c
          if(isw_incmdisk.ne.0) then
c
c             Include the mass of the disk inward of R to the total 
c             mass inducing the gravity.
c
              disk_mgrav(1)  = global_mstar
              do ir=2,grid_nr
                  disk_mgrav(ir) = disk_mgrav(ir-1) + 0.5 * 
     %                 ( disk_tprsigma(ir) + disk_tprsigma(ir-1) ) *
     %                 ( grid_r(ir) - grid_r(ir-1) )
              enddo
          else
c
c             Only use the star as the source of gravity
c     
              do ir=1,grid_nr
                  disk_mgrav(ir) = global_mstar
              enddo
          endif
c
c         For idistr=3, update the stellar radius, luminosity and
c         temperature.
c
          if(cloud_idistr.eq.3) then
              call get_rstar()
              global_lstar = get_lstar()
              global_tstar = sqrt(sqrt(global_lstar /
     $                            (4d0*pi*ss*global_rstar**2)))
          endif
c
c         First compute radial velocity
c
          call compute_vr_v0_difcoef()
c
c         Check if we are about to cross the save time
c
c          if(timesave+time_dt.ge.time_tsave) then
c              time_dt = time_tsave - timesave
c              dump = .true.
c          endif
c
c         Check if we are about to cross the time-limit
c
c          if(time_time+time_dt.gt.time_end) then
c              time_dt = time_end - time_time
c              simstop = .true.
c          endif
c
c         Now do the diffusion step for the gas+dust
c
          call update_sigma_diff()
c
c         Warn, when necessary
c
          if(global_flag_lowtemp.eq.1) then
              write(*,*) "WARNING: TOO LOW TEMP"
          endif
c
c         If matter infalling onto disk:
c         Update the mass of the central star by matter falling directly
c         from the envelope onto the star (i.e. into the capture radius
c         of the star).
c
c         DANIEL CLEANUP
          write(*,*)
          write(*,*) 'MDOTS'
          write(*,*) time_time/year, time_dt/year
          write(*,*) ml_mdotcap, ml_mdotdisk, ml_mdotflow
          write(*,*) global_mstar/ms, global_mdisk/ms, global_msyst/ms
          write(*,*) cloud_mass/ms
          write(*,*)
          if (time_time/year > 50000d0) then 
            stop
          endif
          if(cloud_mass.gt.0.d0) then
              global_mstar = global_mstar + time_dt * ml_mdotcap
              global_mdisk = global_mdisk + time_dt * ml_mdotdisk
              global_msyst = global_msyst - time_dt * ml_mdotflow
          endif
c
c         Update the surface density of the disk by the matter falling
c         onto the disk from the envelope.
c         (OLD EXPLICIT WAY BEFORE 22.05.05)
c
          if((cloud_mass.gt.0.d0).and.(ml_implicit.eq.0)) then
              do ir=1,grid_nr
                  disk_tprsigma(ir) = disk_tprsigma(ir) + 
     %                 time_dt * 2*pi*grid_r(ir)*ml_psi(ir)
                  if(global_ndust.ge.1) then
                      disk_tprsigmadust(ir,1) = 
     %                     disk_tprsigmadust(ir,1) +
     %                     time_dt * 2*pi*grid_r(ir)*
     %                     ml_psi(ir) / global_gtd
                  endif
              enddo
          endif
c
c         Find out until which radius we should do crystallization
c
          irs = 1
          if((dpr_mode.ge.1).and.(dpr_mode.lt.100)) then
              do ir=1,grid_nr-4
                  if(disk_temp(ir).gt.dpr_tempcryst) then
                      irs = ir + 1
                  endif
              enddo
          endif
c
c         Now do the radial mixing of the dust species
c         I use for this a subroutine I programmed some time ago
c         for the stationary background disk, but that routine is
c         also valid for evolving disks.
c
          if(global_ndust.gt.0) then 
              do ir=1,grid_nr
                  sigma(ir)    = disk_tprsigma(ir)/(2*pi*grid_r(ir))
                  diffcoef(ir) = global_prandtlinv*disk_nu(ir)
                  vr(ir)       = disk_vr(ir)
              enddo
              if((dpr_mode.ge.1).and.(dpr_mode.lt.100)) then 
c
c                 Crystalline silicate with discrete temperature
c                 activation.
c
                  do idust=1,global_ndust
c
c                     Special radial mixing for crystallization
c
c                     idust=1  :  Amorphous dust
c                     idust=2  :  Crystalline dust
c
                      do ir=1,grid_nr
                          sigmadust(ir) = disk_tprsigmadust(ir,idust)/
     %                                    (2*pi*grid_r(ir))
                      enddo
                      if(idust.eq.1) then 
c
c                         Amorphous dust
c
                          do ir=1,irs-1
                              sigmadust(ir) = 0.d0
                          enddo
                          if((cloud_mass.gt.0.d0).and.
     %                       (ml_implicit.eq.1)) then
                              do ir=irs,grid_nr
                                  src(ir) = ml_psi(ir) / global_gtd
                              enddo
                              do ir=1,irs-1
                                  src(ir) = 0.d0
                              enddo
                          else
                              do ir=1,grid_nr
                                  src(ir) = 0.d0
                              enddo
                          endif
                      elseif(idust.eq.2) then
c
c                         Crystalline dust
c
                          do ir=1,irs-1
                              sigmadust(ir) = sigma(ir) / global_gtd
                          enddo
                          if((cloud_mass.gt.0.d0).and.
     %                       (ml_implicit.eq.1)) then
                              do ir=irs,grid_nr
                                  src(ir) = 0.d0
                              enddo
                              do ir=1,irs-1
                                  src(ir) = ml_psi(ir) / global_gtd
                              enddo
                          else
                              do ir=1,grid_nr
                                  src(ir) = 0.d0
                              enddo
                          endif
                      else
                          stop 89236
                      endif
                      if(irs.eq.1) then
                          dummy = 1d99
                      else
                          dummy = sigmadust(irs-1)
                      endif
                      vr(irs) = 0.d0
c
c                     Do the radial mixing
c
                      if(dpr_mode.eq.1) then
                          call dustradmix(grid_nr-irs+1,
     %                       grid_r(irs),sigma(irs),sigmadust(irs),
     %                       src(irs),vr(irs),diffcoef(irs),
     %                       time_dt,dummy)
                      elseif(dpr_mode.eq.2) then
                          disk_iv0(irs) = 0.d0
                          call dustradmix_old(grid_nr-irs+1,
     %                       grid_r(irs),sigmadust(irs),src(irs),
     %                       disk_iv0(irs),disk_dcoef(irs),
     %                       time_dt,dummy)
                      else
                          stop 88299
                      endif
c
c                     Put things back into the main array
c
                      do ir=1,grid_nr
                          disk_tprsigmadust(ir,idust) = sigmadust(ir) *
     %                         (2*pi*grid_r(ir))
                      enddo
                  enddo
c
              elseif(dpr_mode.ge.100) then
c
c                 Crystalline silicate with time-dependent 
c                 crystallization
c
c                     idust=1  :  Amorphous dust
c                     idust=2  :  Crystalline dust
c
                  if(global_ndust.ne.2) stop 9892
c
c                 Compute the crystallization rate
c
                  do ir=1,grid_nr
                      xp  = dpr_tempcryst/disk_temp(ir)
                      if(xp.gt.1.d2) then 
                          crystrate(ir) = 0.d0
                      else
                          crystrate(ir) = dpr_nu * exp(-xp)
                      endif
                  enddo
c                 
c                 Crystalline dust
c
                  do ir=1,grid_nr
                      sigmadust(ir) = disk_tprsigmadust(ir,2)/
     %                                    (2*pi*grid_r(ir))
                  enddo
                  dummy = sigma(ir) / global_gtd
c
c                 Do the radial mixing
c
                  call crystradmix(grid_nr,grid_r,sigma,sigmadust,
     %                       vr,diffcoef,crystrate,global_gtd,
     %                       time_dt,dummy)
c
c                 Put things back into the main array
c                     
                  do ir=1,grid_nr
                      disk_tprsigmadust(ir,2) = sigmadust(ir) *
     %                         (2*pi*grid_r(ir))
                      disk_tprsigmadust(ir,1) = 
     %                         ((sigma(ir)/global_gtd)-sigmadust(ir)) *
     %                         (2*pi*grid_r(ir))
                  enddo
c
              elseif(dpr_mode.le.0) then
c
c                 Old mode
c
                  do idust=1,global_ndust
c
c                     Normal routine for radial mixing
c     
c                     If infall, then also add to the dust of course
c
                      if((cloud_mass.gt.0.d0).and.
     %                   (ml_implicit.eq.1)) then
                          do ir=1,grid_nr
                              src(ir) = ml_psi(ir) / global_gtd
                          enddo
                      else
                          do ir=1,grid_nr
                              src(ir) = 0.d0
                          enddo
                      endif
c
c                     Do mixing
c     
                      call dustradmix(grid_nr,grid_r,sigma,sigmadust,
     %                     src,disk_vr,diffcoef,time_dt,1d99)
c
c                     Put things back into the main array
c
                      do ir=1,grid_nr
                          disk_tprsigmadust(ir,idust) = sigmadust(ir) *
     %                         (2*pi*grid_r(ir))
                      enddo
                  enddo
              else
                  stop 771
              endif
          endif
c
c         Update mass of the star (if active)
c
          if(isw_growstar.ne.0) then
              global_mstar = global_mstar + time_dt * abs(disk_mdot(1))
              global_mdisk = global_mdisk - time_dt * abs(disk_mdot(1))
          endif
c
c         Update the total amount of outflowing angular momentum through
c         the inner and outer edge
c         NOTE 1: For the outer edge this is inaccurate, since there is
c                 no explicit no-friction BC imposed. But the purpose of
c                 this variable is anyway to make sure it is near 0. 
c         NOTE 2: I take the absolute value of this flux of angular 
c                 momentum, since both these variables are only used for
c                 testing.
c
          global_j_in  = global_j_in + abs(disk_mdot(1)) * 
     %                   sqrt(GG*disk_mgrav(1)*grid_r(1)) * 
     %                   time_dt
          global_j_out = global_j_out + abs(disk_mdot(grid_nr+1)) * 
     %                   sqrt(GG*disk_mgrav(grid_nr)*grid_r(grid_nr)) * 
     %                   time_dt
c
c         Get the midplane density and the vertical density profile.
c         Take: a Gaussian for the disk and the Terebey, Shu & Cassen
c         (1984) or Cassen & Moosman (1981) solution for the envelope.
c         Also take the outflow into account.
c
          if(dovert.eq.1) then
              do ir=1,grid_nr
                  rsph = grid_r(ir)
                  theta = 0.5d0*pi
                  scale = sqrt(kk*disk_temp(ir)*(grid_r(ir)**3)*
     $                mpeffinv*GGinv/disk_mgrav(ir))
                  rhom = disk_tprsigma(ir) /
     $                (2d0 * pi * grid_r(ir) * scale * sqrt(2d0*pi))
                  disk_midden(ir) = max(1d-50,rhom*mpeffinv)
                  disk_midprs(ir) = disk_midden(ir) * kk*disk_temp(ir)
                  call get_vertden(grid_r(ir),0d0,rsph,theta,
     $                             time_time,-1d0,0)
                  env_midden(ir) = max(1d-50,mden*mpeffinv)
                  env_midprs(ir) = prs
c
c     Now loop over the vertical grid to get the vertical density
c     profile.
c
                  do iz=1,grid_nz
                      rsph = sqrt(grid_r(ir)**2+grid_z(iz)**2)
                      theta = atan2(grid_r(ir),grid_z(iz))
                      disk_vertden(ir,iz) = max(1d-50,rhom * mpeffinv *
     $                                exp(-0.5d0*(grid_z(iz)/scale)**2))
                      disk_vertprs(ir,iz) = disk_vertden(ir,iz) *
     $                                      kk*disk_temp(ir)
                      call get_vertden(grid_r(ir),grid_z(iz),rsph,theta,
     $                                 time_time,-1d0,0)
                      env_vertden(ir,iz) = max(1d-50,mden*mpeffinv)
                      env_vertdenf(ir,iz) = max(1d-50,mdenf*mpeffinv)
                      env_vertprs(ir,iz) = prs
                      disk_ur(ir,iz) = disk_vr(ir)
                      disk_uth(ir,iz) = 0.0
                      disk_uphi(ir,iz) = sqrt(GG*global_mstar/(rsph))
c                      write(*,*) ur,uth,uphi,max(1d-50,mdenf*mpeffinv)
                      env_ur(ir,iz) = ur
                      env_uth(ir,iz) = uth
                      env_uphi(ir,iz) = uphi
c                      write(*,*) env_ur(ir,iz),env_uth(ir,iz)
c     $                           env_uphi(ir,iz)
                  enddo
              enddo
          endif
c
c         Then save when necessary
c
          timesave = timesave + time_dt
          itsave   = itsave+1
          if(dump) then
c              if(time_time.gt.1.27d12.and.time_time.lt.1.30d12) then
                  write(*,*) 'Saving data at time step ',it,
     $                       time_time/year,' year'
                  call save_data()
c              endif
              timesave=0.d0
              itsave=0
    
c              write(*,*) 'Masses :'
c              write(*,*) '    star = ',global_mstar/MS,'MS'
c              write(*,*) '    disk = ',global_mdisk/MS,'MS'
c              write(*,*) 'envelope = ',(cloud_mass-global_mstar-
c     $                                global_mdisk)/MS,'MS'
c              write(*,*)

          endif
c          if(.not.dump.and.time_time.gt.4.17d12.and.
c     $                                        time_time.lt.4.50d12) then
c              write(*,*) 'Saving data at time step ',it,
c     $                    time_time/year,' year (t>3.38e12)'
c              call save_data()
c          endif
c
c         Update time
c
          time_time = time_time + time_dt
c
c         If variable save time, then recompute the time_tsave
c
          if(varsdt.and.dump) then
              timesave   = timesave + time_tsave
              time_tsave = time_tsave * time_inctsave
c Added by DH Aug 18, 2011:
c Need a constant increase of tsave by a function of accretion time:
c better and faster for radiative transfer...
              if ( (time_time).lt.(1.2 * (GG*cloud_mass) / 
     $          (m0*cloud_asound*cloud_asound*cloud_asound) ) ) then
                  time_tsave = min(time_tsave, ( (time_time - time_dt) + 
     $            0.01 * (GG*cloud_mass) / 
     $            (m0*cloud_asound*cloud_asound*cloud_asound) ) )
              endif
c              write(*,*) '***********',time_tsave/year
          endif
c
c         Check if we're done
c
          if(simstop) goto 20
c
c         Reset things
c
          dump=.false.
c
      enddo
c
 20   continue
      if(.not.varsdt) then
          write(*,*) 'Saving last data'
          call save_data()
      endif
      write(*,*) 'Finished,     Time = ',time_time/year,
     %           ' year'
c
c     Done
c
      end


c     --------------------------------------------------------------
c                   FUNCTION: TEST VALIDITY OF NUMBER
c
c     0 = Number is okay
c     1 = Number is INF
c     2 = Number is NAN
c     --------------------------------------------------------------
      function number_invalid(a)
      implicit none
      doubleprecision a,b,c
      integer number_invalid
      logical div,sub
c
      b=a*2.d0
      b=b/2.d0
      c=a-1.d100
c
      div = (b.eq.a)
      sub = (c.lt.a)
c
      if(div.and.sub) then
          number_invalid = 0
      elseif(div) then
          number_invalid = 1
      else
          number_invalid = 2
      endif
c
      return
      end


c     --------------------------------------------------------------
c                      NUMERICAL RECIPES: ZBRENT
c     --------------------------------------------------------------
      FUNCTION zbrent(func,x1,x2,tol,w,rr)
      INTEGER ITMAX,w
      DOUBLEPRECISION zbrent,tol,x1,x2,func,EPS,rr
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.d-8)
      INTEGER iter
      DOUBLEPRECISION a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a,0d0,0)
      fb=func(b,0d0,0)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then
          write(*,*) 'root must be bracketed for zbrent'
          stop
      endif
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b,0d0,0)
11    continue
      write(*,*) 'zbrent exceeding maximum iterations'
      stop
      zbrent=b
      return
      END

C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..
c
c     --------------------------------------------------------------
c                 READ THE DISK EVOLUTION RESULTS
c     --------------------------------------------------------------
      subroutine read_data()
      implicit none
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
c
      integer ir,iz,it,idust,i1,i2
      doubleprecision dum,dum2
c
   10 format(E13.6)
   11 format(2(E13.6,1X))
   12 format(3(E13.6,1X))
c
      de_ngrid = grid_nr
      open(unit=1,file='time.info',status='old')
      read(1,*) de_ntime
      if(de_ntime.gt.FRSIZE_TIME) stop 42699
      if(de_ntime.ne.time_nr_save.and.time_nr_save.gt.0) stop 42700
      close(1)
c
      open(unit=1,file='sigma.dat',status='old')
      read(1,*) i1,i2
      if(i1.ne.grid_nr.or.i2.ne.global_ndust) stop 42701
      do it=1,de_ntime
          read(1,*)
          do ir=1,de_ngrid
              read(1,10) dum
              if(dum.le.global_sigma_min) then
                  de_lnsigma(it,ir) = log(global_sigma_min)
              else
                  de_lnsigma(it,ir) = log(dum)
              endif
          enddo
          if(global_ndust.gt.0) then
              do idust=1,global_ndust
                  do ir=1,de_ngrid
                      read(1,10) dum
                  enddo
              enddo
          endif
      enddo
      close(1)
c
      open(unit=1,file='velo.dat',status='old')
      read(1,*) i1,i2
      if(i1.ne.grid_nr.or.i2.ne.global_ndust) stop 42702
      do it=1,de_ntime
          read(1,*)
          do ir=1,de_ngrid
              read(1,10) de_vr(it,ir)
          enddo
!           if(global_ndust.gt.0) then
!               do idust=1,global_ndust
!                   do ir=1,de_ngrid
!                       read(1,10) dum
!                   enddo
!               enddo
!           endif
      enddo
      close(1)
c
      open(unit=1,file='temperature.dat',status='old')
      read(1,*) i1
      if(i1.ne.grid_nr) stop 42703
      do it=1,de_ntime
          read(1,*) 
          do ir=1,de_ngrid
              read(1,10) de_temp(it,ir)
          enddo
      enddo
      close(1)
c
      open(unit=1,file='mdot.dat',status='old')
      read(1,*) i1
      if(i1.ne.grid_nr) stop 42704
      do it=1,de_ntime
          read(1,*) 
          read(1,10) de_mdot_d2s(it)
          do ir=2,de_ngrid
              read(1,10) dum
          enddo
      enddo
      close(1)
c
      open(unit=1,file='mgrav.dat',status='old')
      read(1,*) i1
      if(i1.ne.grid_nr) stop 42705
      do it=1,de_ntime
          read(1,*) 
          do ir=1,de_ngrid
              read(1,10) de_mgrav(it,ir)
          enddo
      enddo
      close(1)
c
      open(unit=1,file='surface.dat',status='old')
      read(1,*) i1
      if(i1.ne.grid_nr) stop 42706
      do it=1,de_ntime
          read(1,*) 
          do ir=1,de_ngrid
              read(1,11) de_surface(it,ir),dum
          enddo
      enddo
      close(1)
c
      open(unit=1,file='infall.dat',status='old')
      read(1,*)
      do it=1,de_ntime
          read(1,12) de_mdot_c2s(it),dum,dum2
      enddo
      close(1)
c
      open(unit=1,file='mstar.dat',status='old')
      read(1,*)
      do it=1,de_ntime
          read(1,10) de_mstar(it)
      enddo
      close(1)
c
      open(unit=1,file='mdisk.dat',status='old')
      read(1,*)
      do it=1,de_ntime
          read(1,10) de_mdisk(it)
      enddo
      close(1)
c
      open(unit=1,file='time.dat',status='old')
      read(1,*)
      do it=1,de_ntime
          read(1,10) de_time(it)
          if(it.eq.1) then
              de_lntime(it) = 0d0
          else
              de_lntime(it) = log(de_time(it))
          endif
      enddo
      close(1)
c
c     Done
c
      end



c     ------------------------------------------------------------------
c     This subroutine computes the stellar radius as a function of time
c     and stellar mass (Palla & Stahler 1991).
c     ------------------------------------------------------------------
      subroutine get_rstar()
c
      implicit none
      integer j
      doubleprecision rstarps,eps
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
c
      j = 1
      call hunt(mps,NPS,global_mstar,j)
      if((j.le.0).or.(j.ge.NPS)) then
          write(*,*) j,global_mstar
          stop 8370
      endif
      eps = (global_mstar-mps(j)) / (mps(j+1)-mps(j))
      if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 8371
      rstarps = (1.d0-eps) * rps(j) + eps * rps(j+1)
      if(time_time.le.2d4*year) then
          global_rstar = 5d0 * au
      else if(time_time.lt.2.01d4*year) then
          global_rstar = 5d0 * au * (1d0-sqrt((time_time*yearinv-2d4)/
     $                   1d2)) + rstarps
      else
          global_rstar = rstarps
      endif
c
      return
      end
c
c
c
c     ------------------------------------------------------------------
c     This function computes the stellar luminosity as a function of
c     time and stellar mass (D'Antona & Mazzitelli 1994). A negative
c     value for mmmd means we're in the T Tauri phase, where we set the
c     accretion luminosity to zero.
c     ------------------------------------------------------------------
      function get_lstar()
c
      implicit none
      integer j,k
      doubleprecision get_lstar,mmmd,laccs1,laccs2
      doubleprecision ustar,dum,lum0,laccs,lphot,time5,eps,epst
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
c
c     Compute the accretion luminosity.
c
      mmmd = -1d0
      ustar = global_rstar / ml_rcentr
      lum0 = GG * ml_mcentr * mdot / global_rstar
      if(ustar.ge.1d0) then
          laccs1 = lum0
          laccs2 = 0d0
          laccs = lum0
      else
          dum = sqrt(1d0-ustar)
          mmmd = (1d0 - global_mstar/ml_mcentr) / (1d0 - etadisk)
          laccs1 = ((6d0*ustar - 2d0 + dum*(2d0-5d0*ustar))/(6d0*ustar))
     $             * lum0
          laccs2 = ((1d0-(1d0-etadisk)*mmmd)*(1d0-(1d0-etadisk)*dum)*
     $             0.5d0*etastar) * lum0
          laccs = laccs1 + laccs2
      endif
c
c     Compute the photospheric luminosity.
c
      lphot = 0d0
      if(global_mstar.gt.0.02d0*MS) then
          time5 = time_time - 0.4d0 * tcloudm
          if((time5.gt.0d0).and.(time5.le.year7d4)) then
              j = 1
              call hunt(dmmass,NDMM,global_mstar,j)
              if((j.le.0).or.(j.ge.NDMM)) then
                  j = NDMM-1
                  eps = 1d0
                  write(*,*) j,NDMM,global_mstar/MS
                  stop 8360
              else
                  eps = (global_mstar-dmmass(j)) /
     $                  (dmmass(j+1)-dmmass(j))
                  if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 8361
              endif
              lphot = (1d0-eps) * dmlum(j,1) + eps * dmlum(j+1,1)
              lphot = lphot * (time5/year7d4)**5
          else if(time5.gt.year7d4) then
              lphot = 0d0
              j = 1
              call hunt(dmmass,NDMM,global_mstar,j)
              if(j.le.0) stop 8362
              eps = (global_mstar-dmmass(j)) / (dmmass(j+1)-dmmass(j))
              if(j.ge.NDMM) then
                  j = NDMM-1
                  eps = 1d0
                  if(dmpassmax.eq.0) then
                      write(*,*) 'Stellar mass in get_lstar is higher ',
     $                    "than the highest stellar mass in D'Antona &"
                      write(*,*) 'Mazzitelli (1994). The code will ',
     $                    'continue with the luminosity table for a'
                      write(*,*) 'stellar mass of 2.5 M_sun.'
                      dmpassmax = 1
                  endif
              endif
              if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 8363
              k = 1
              call hunt(dmtime,NDMT,time5,k)
              if((k.le.0).or.(k.ge.NDMT)) stop 8364
              epst = (time5-dmtime(k)) / (dmtime(k+1)-dmtime(k))
              if((epst.lt.0.d0).or.(epst.gt.1.d0)) stop 8365
              lphot = (1d0-epst) * ((1d0-eps) * dmlum(j,k) + eps *
     $                dmlum(j+1,k)) + epst * ((1d0-eps) * dmlum(j,k+1) +
     $                eps * dmlum(j+1,k+1))
          else
              lphot = 0d0
          endif
          get_lstar = laccs + lphot
      else
          get_lstar = laccs
      endif
c
      return
      end
c
c
c
c     ------------------------------------------------------------------
c     This subroutine computes the infall trajectory for a single
c     parcel. If the dowritechem flag is active, it writes the output to
c     a file that can be read by Steve Doty's chemical code.
c     ------------------------------------------------------------------
      subroutine make_trajectory(ip,dowritechem,do_tau,do_cdh)
c
      implicit none
      integer ip,dowritechem,do_tau,j,itime,jtime,zmethod,i,nroot
      integer do_cdh
      doubleprecision tnew,rnew,thnew,mdennew,rnewinv,tnewinv,dt,xnew
      doubleprecision mu0,ustar,surf,rcen,nden,ndenenv,ndendisk,nden0
      doubleprecision znew,told,rold,thold,rsave,thsave,eps,deeps,deepsn
      doubleprecision erfi,zold,surferf,dr,dum,scaleold,rrc,xx,sedum
      doubleprecision costh,check,munew,tflow
      doubleprecision prsdisk, prsenv
      doubleprecision get_dust_temp,temp,lntnew,xcub(3),acub(4)
      character*4 loc
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
c
c     Initialize several parameters. The time (tnew) is initialized to
c     the point where the expansion wave reaches the initial position
c     (r = c_s t).
c
      rnew = r0(ip)
      thnew = th0(ip)
      xnew = rnew * sin(thnew)
      znew = rnew * cos(thnew)
      tnew = 1.00000001d0 * rnew * cloud_asound_inv
      lntnew = log(tnew)
      nden0 = mpeffinv * 0.5d0 * piinv * GGinv *
     $        (cloud_asound/rnew)**2
      rold = rnew
      thold = thnew
      told = tnew
      zold = znew
      scaleold = -1d0
      outflow = 0
      time_time = 0d0
c
c     If the parcel's initial angle is less than the maximum outflow
c     angle, compute the time at which the outflow angle equals the
c     parcel's initial angle.
c
      xnew = rnew*sin(thnew)
      znew = rnew*cos(thnew)
      if(znew.gt.0.191d0*au*(xnew/au)**1.5d0) then
          tflow = (znew/(0.191d0*au*(xnew/au)**1.5d0))**(-f1o3) *
     $            tcloudm
      else
          tflow = 10d0 * tmax
      endif
c
c     If tnew exceeds the maximum output time step or if the outflow
c     will reach the parcel before it begins to fall in, write the
c     initial position for each time step and skip straight to the end.
c     Terminate if and when the outflow reaches the parcel.
c
      if(tnew.gt.tmax.or.tnew.gt.tflow) then
          step = 2
          tstore(1) = 0d0
          rstore(1) = rnew
          thstore(1) = thnew
          nstore(1) = nden0
          tstore(2) = min(tmax,tflow)
          rstore(2) = rnew
          thstore(2) = thnew
          nstore(2) = nden0
          do j=1,de_ntime
              if(nparc.lt.100000.or.j.eq.102.or.j.eq.104) then
                  if(de_time(j).lt.tflow) then
                      write(10,100) j,de_time(j),rnew,thnew,nden0
                  else
                      write(10,100) j,de_time(j),1d0,0.25d0*pi,1d0
                  endif
              endif
          enddo
          itime = de_ntime + 1
          goto 1020
      endif
c
c     Find the initial time step on the disk evolution grid. Set the
c     stellar mass accordingly.
c
      jtime = 1
      call hunt(de_lntime,de_ntime,lntnew,jtime)
      if((jtime.le.1).or.(jtime.ge.de_ntime)) then
          write(*,*) jtime,de_ntime
          stop 8310
      endif
      deeps = (lntnew-de_lntime(jtime)) /
     $        (de_lntime(jtime+1)-de_lntime(jtime))
      if((deeps.lt.0.d0).or.(deeps.gt.1.d0)) stop 8311
      mstar = (1.d0-deeps) * de_mstar(jtime) + deeps * de_mstar(jtime+1)
c
c     Locate tnew on the diskevol time grid.
c
      itime = 1
      call hunt(de_time,de_ntime,tnew,itime)
      if((itime.le.0).or.(itime.ge.de_ntime)) stop 8350
c
c     Write the initial values at all time steps up to itime.
c
      do j=1,itime
          if(nparc.lt.100000.or.j.eq.102.or.j.eq.104) then
              temp = get_dust_temp(de_time(j),rnew,thnew)
              write(10,100) j,de_time(j),rnew,thnew,nden0,temp
          endif
      enddo
  100 format(I5,2(1PE16.8),0PF12.8,1PE16.8,1PE13.5)
  101 format(I5,2(1PE16.8),0PF12.8,6(1PE16.8),1X,A5,I2)
      itime = itime + 1
c
c     If the infall starts after the first hydrostatic core phase
c     [assumed to be 2.01e4 year], add two steps to the
c     storage arrays to get the correct density and temperature profile.
c
      if(tnew.gt.2.01d4*year) then
          tstore(1) = 0.999d0*2.01d4*year
          rstore(1) = rnew
          thstore(1) = thnew
          nstore(1) = nden0
          if(write_traj.eq.1) then
              write(20,101) 0,0.999d0*2.01d4*year,rnew,thnew,nden0,
     $                      0d0,0d0,0d0,0d0,zerotemp,'init',0
          endif
          tstore(2) = 2.01d4*year
          rstore(2) = rnew
          thstore(2) = thnew
          nstore(2) = nden0
          if(write_traj.eq.1) then
              temp = get_dust_temp(2.01d4*year,rnew,thnew)
              write(20,101) 1,2.01d4*year,rnew,thnew,nden0,0d0,0d0,
     $                      0d0,0d0,temp,'init',0
          endif
          step = 3
      else
          step = 1
      endif
      if(write_traj.eq.1) then
          temp = get_dust_temp(tnew,rnew,thnew)
          write(20,101) step-1,tnew,rnew,thnew,nden0,0d0,0d0,0d0,0d0,
     $                  temp,'init',0
      endif
c
c     For any reasonable set of initial conditions and position, the
c     parcel will start in the TSC84 or CM81 regime. Check whether that
c     is indeed the case.
c
      rcen = 0.0625d0*cloud_asound*((m0*tnew)**3)*(cloud_omega**2)
      if(rcen.ge.rnew) then
          write(*,*) 'Centrifugal radius exceeds expansion wave at ',
     $               'the beginning of the infall of parcel ip=',ip
          stop
      endif
      loc = 'tsc_'
      if(rnew/(cloud_asound*tnew).lt.(cloud_omega*tnew)**2) then
          costh = cos(thnew)
          rrc = rnew / rcen
          if(rrc.eq.1d0) then
              mu0 = costh**f1o3
          else if(costh.eq.0d0) then
              if(rrc.lt.1d0) then
                  mu0 = sqrt(1d0-rrc)
              else
                  mu0 = 0d0
              endif
          else
              acub(4) = 1d0
              acub(3) = 0d0
              acub(2) = rrc-1d0
              acub(1) = -rrc*costh
              call cbcrt(acub,xcub)
              mu0 = xcub(1)
          endif
          if(costh.gt.0d0) then
              check = mu0 * (1d0-(1d0-mu0**2)/rrc)
              if(abs(1d0-check/costh).gt.1d-6) stop 3362
          endif
          loc = 'cm__'
      endif
c
c     Start the infall.
c
      j = 1
      ux = 0d0
      uz = 0d0
      mden = 0d0
      outflow = 0
 1000 continue
c
c     Compute 1/rnew and 1/tnew.
c
      rnewinv = 1d0/rnew
      tnewinv = 1d0/tnew
      lntnew = log(tnew)
c
c     Find the current time step on the disk evolution grid. Set the
c     mass of the star and the disk accordingly.
c
      call hunt(de_lntime,de_ntime,lntnew,jtime)
      if((jtime.le.1).or.(jtime.ge.de_ntime)) stop 8312
      deeps = (lntnew-de_lntime(jtime)) /
     $        (de_lntime(jtime+1)-de_lntime(jtime))
      deepsn = 1d0 - deeps
      if((deeps.lt.0.d0).or.(deeps.gt.1.d0)) stop 8313
      mstar = deepsn * de_mstar(jtime) + deeps * de_mstar(jtime+1)
      mdisk = deepsn * de_mdisk(jtime) + deeps * de_mdisk(jtime+1)
c
c     The method of computing the velocity depends on the type of regime
c     the parcel is currently in. If we're in the TSC84 or CM81 regime,
c     we also need to get the disk density to determine if the
c     parcel has entered the disk. The 'velden' routines overwrite the
c     value of the global parameter mden, so we have to store the
c     densities in other variables.
c
      if(loc.eq.'tsc_') then
          call get_tscvelden(rnew,thnew,tnew,tnewinv,0)
          ndenenv = mden * mpeffinv
          prsenv = prs
          call get_diskvelden(xnew,znew,tnew,jtime,deeps,deepsn,0)
          ndendisk = mden * mpeffinv
          prsdisk = ndendisk * kk * dktemp
      endif
      if(loc.eq.'cm__') then
          call get_cmvelden(rnew,thnew,tnew,mu0,rnewinv,0)
          ndenenv = mden * mpeffinv
          prsenv = prs
          call get_diskvelden(xnew,znew,tnew,jtime,deeps,deepsn,0)
          ndendisk = mden * mpeffinv
          prsdisk = ndendisk * kk * dktemp
      endif
c
c     Check if the parcel entered the outflow. If so, delete the last
c     step and stop the calculation for this parcel.
c
      if((loc.eq.'tsc_'.or.loc.eq.'cm__').and.outflow.eq.1) then
          step = step - 1
          goto 1020
      endif
c
c     If we were in the TSC84 or CM81 regime, did we end up below the
c     surface of the disk? If so, the parcel is in the disk from now on.
c     Compute surferf, which we need for the vertical motion later.
c
      if((loc.eq.'tsc_').or.(loc.eq.'cm__')) then
c          if(znew.lt.surface) then
c          if(ndendisk.gt.ndenenv) then
          if (prsdisk.gt.prsenv) then
              loc = 'disk'
              write(12,*) ip,tnew,rnew,thnew
              dr = xepsn * grid_dr(jgrid) + xeps * grid_dr(jgrid+1)
              surferf = surf2 * erf(znew/(scale2*sqrt(2d0)))
          else
              nden = ndenenv
          endif
      endif
c
c     Compute the velocity and density for the disk. Do this also if we
c     crossed into the disk regime in this iteration. In that case, we
c     don't need to recompute the density, but since this only happens
c     once for every parcel, we're not wasting a lot of time.
c
      if(loc.eq.'disk') then
          call get_diskvelden(xnew,znew,tnew,jtime,deeps,deepsn,1)
          nden = mden * mpeffinv
      endif
c
c     Check if the parcel entered the outflow. If so, delete the last
c     step and stop the calculation for this parcel.
c
      if((loc.eq.'tsc_'.or.loc.eq.'cm__').and.outflow.eq.1) then
          step = step - 1
          goto 1020
      endif
c
c     Now that we have the velocities, we can compute the new radial and
c     polar position. We first try a fairly large time step and reduce
c     it until the position changes by at most 1%, to obtain an accurate
c     trajectory. We check for changes in the x and z direction in the
c     TSC84 and CM81 regimes, and for changes in only the x direction in
c     the disk regime.
c
      dt = 1d-1 * tnew
      if(loc.eq.'disk') then
 1010     dt = 1d-1 * dt
          if(abs(ux*dt).gt.abs(1d-2*xnew)) goto 1010
      else if(loc.eq.'cm__') then
 1012     dt = 1d-1 * dt
          if(abs(ur*dt).gt.abs(1d-2*rnew).or.
     $       (abs(dt*uth*rnewinv).gt.abs(1d-2*uth))) goto 1012
      else
 1011     dt = 1d-1 * dt
          if((abs(ux*dt).gt.abs(1d-2*xnew)).or.
     $        (abs(uz*dt).gt.abs(1d-2*znew))) goto 1011
      endif
      if(dt*tnewinv.lt.1d-14) then
          write(*,500) 'dt/tnew < 1d-14 at tnew =',tnew
          write(*,*) 'ip,loc = ',ip,loc
          write(*,500) 'x,z,ux,uz = ',xnew,znew,ux,uz
  500     format(A,4(1PE13.5))
          stop
      endif
c
c     In the TSC84 and CM81 regimes, use ux and uz to get the new
c     position. In the disk regime, use ux for the radial position and
c     get the vertical position from the erf formula.
c
      zmethod = 0
      tnew = tnew + dt
      if(loc.eq.'cm__') then
          rnew = rnew + ur*dt
          munew = mu0*(1d0-(1d0-mu0**2)*rdisk/rnew)
          if(munew.le.0d0) then
              thnew = 0.4999d0*pi
          else
              thnew = acos(munew)
          endif
          xnew = rnew * sin(thnew)
          znew = rnew * cos(thnew)
      else
          xnew = xnew + ux*dt
          if(loc.eq.'disk') then
              lntnew = log(tnew)
              call hunt(de_lntime,de_ntime,lntnew,jtime)
              if((jtime.le.0).or.(jtime.ge.de_ntime)) then
                  write(*,*) jtime,tnew/year,tmax/year
                  stop 83121
              endif
              deeps = (lntnew-de_lntime(jtime)) /
     $                (de_lntime(jtime+1)-de_lntime(jtime))
              deepsn = 1d0 - deeps
              if((deeps.lt.0.d0).or.(deeps.gt.1.d0)) stop 83131
              call get_diskvelden(xnew,znew,tnew,jtime,deeps,deepsn,0)
              dr = xepsn * grid_dr(jgrid) + xeps * grid_dr(jgrid+1)
              dum = surferf / surf2
              if(dum.lt.1d0) then
                  znew = scale2 * sqrt(2d0) * erfi(dum)
                  zmethod = 1
              else
                  if(scaleold.lt.0d0) then
                      write(*,*) ip,scaleold
                      stop 442
                  endif
                  znew = scale2 * zold / scaleold
                  zmethod = 2
              endif
              uz = (znew-zold) / (tnew-told)
          else
              znew = znew + uz*dt
          endif
          if(xnew.lt.0d0.or.znew.lt.0d0) stop 565
          rnew = sqrt(xnew**2+znew**2)
          thnew = atan2(xnew,znew)
      endif
c
c     Make sure we don't cross the pole or the equator.
c
      if(thnew.gt.0.4999d0*pi) then
          thnew = 0.4999d0*pi
          xnew = rnew * sin(thnew)
          znew = rnew * cos(thnew)
      endif
      if(thnew.lt.0.0001d0) then
          thnew = 0.0001d0
          xnew = rnew * sin(thnew)
          znew = rnew * cos(thnew)
      endif
c
c     If we crossed a diskevol time step, compute and write the position
c     corresponding to this step. Also check if we perhaps crossed two
c     or more diskevol steps at once.
c
 1030 if(tnew.gt.de_time(itime)) then
          eps = (de_time(itime)-told) / (tnew-told)
          if((eps.lt.0.d0).or.(eps.gt.1.d0)) then
              write(*,*) ip,eps,de_time(itime),told,tnew,itime
              write(*,*) de_time(de_ntime)/year,tcloudm/year
              stop 8341
          endif
          rsave = (1.d0-eps) * rold + eps * rnew
          thsave = (1.d0-eps) * thold + eps * thnew
          if(nparc.lt.100000.or.itime.eq.102.or.itime.eq.104) then
              write(10,100) itime,de_time(itime),rsave,thsave,nden
          endif
          itime = itime + 1
          if(itime.le.de_ntime) goto 1030
      endif
c
c     If we were in the TSC84 regime, did we cross into the CM81 regime?
c
      if(loc.eq.'tsc_') then
          xx = rnew / (cloud_asound*tnew)
          if(xx.lt.(cloud_omega*tnew)**2) then
              loc = 'cm__'
              rrc = rnew / rdisk
              costh = cos(thnew)
              if(rrc.eq.1d0) then
                  mu0 = costh**f1o3
              else if(costh.eq.0d0) then
                  if(rrc.lt.1d0) then
                      mu0 = sqrt(1d0-rrc)
                  else
                      mu0 = 0d0
                  endif
              else
                  acub(4) = 1d0
                  acub(3) = 0d0
                  acub(2) = rrc-1d0
                  acub(1) = -rrc*costh
                  call cbcrt(acub,xcub)
                  mu0 = xcub(1)
              endif
              if(costh.gt.0d0) then
                  check = mu0 * (1d0-(1d0-mu0**2)/rrc)
                  if(abs(1d0-check/costh).gt.1d-6) stop 3332
              endif
          endif
      endif
c
c     Get the correct density.
c
      if(loc.eq.'tsc_') then
          call get_tscvelden(rnew,thnew,tnew,tnewinv,0)
      else if(loc.eq.'cm__') then
          call get_cmvelden(rnew,thnew,tnew,mu0,rnewinv,0)
      else
          call get_diskvelden(xnew,znew,tnew,jtime,deeps,deepsn,0)
      endif
      nden = mden * mpeffinv
c
c     Store the results from this step.
c
      tstore(step) = tnew
      rstore(step) = rnew
      thstore(step) = thnew
      nstore(step) = nden
      if(write_traj.eq.1) then
          temp = get_dust_temp(tnew,rnew,thnew)
          write(20,101) step,tnew,rnew,thnew,nden,scale2,surf2,ux,uz,
     $                  temp,loc,zmethod
      endif
c
c     Did we reach the star or cross the maximum time? If not, save the
c     current values and do another loop. If we did, stop the iteration
c     for this parcel.
c
      if((rnew.lt.rin).or.(tnew.gt.tmax)) goto 1020
      told = tnew
      rold = rnew
      thold = thnew
      zold = znew
      scaleold = scale2
c
c     Did we reach the max number of steps? If so, report the problem
c     and stop the iteration. If not, go to the next step.
c
      if(step.eq.NSTORE) then
          write(*,*) 'Max. number of steps reached. You may want to ',
     $               'recompile with a larger NSTORE.'
          write(*,*) ip,loc,ux,uz
          write(*,*) step, nstore, rold/au, zold/au
!           stop 900
c        DH Change Dec 2013
          goto 1020 
      endif
      step = step + 1
      goto 1000
 1020 continue
c
c     Write values for the remaining steps (radial distance of 1 cm
c     and angle of pi/4).
c
      if(itime.le.de_ntime) then
          do j=itime,de_ntime
              if(nparc.lt.100000.or.j.eq.102.or.j.eq.104) then
                  write(10,100) j,de_time(j),5d0,0.78539816d0,1d0
              endif
          enddo
      endif
c
c     Create an input file for the chemistry with the temperature,
c     density and time derivative of the density at a number of time
c     steps.
c
      if(dowritechem.eq.1) then
          call write_chem_input(ip,do_tau,do_cdh)
      endif
c
      return
      end
c
c
c
c     ------------------------------------------------------------------
c     This subroutine computes the radial and angular velocity and the
c     density of an infalling parcel in the TSC84 regime.
c     ------------------------------------------------------------------
      subroutine get_tscvelden(r,th,time,timeinv,w)
c
      implicit none
      integer w,j
      doubleprecision r,th,time,timeinv,x,xlog,tau,tau2,costh,legendre
      doubleprecision eps,epsn,f,fn,xinv,xinvsq,xinv1p5,xinv2p5,sinth
      doubleprecision v0,vq,vm,wq,a0,aq,am,outangle,gamma,aout
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
c
c     Compute the dimensionless collapse variables (S77,TSC84).
c
      if(time.gt.0d0) then
          x = r * cloud_asound_inv * timeinv
          xlog = log(x)
      else
          x = 10d0
      endif
      if(x.ge.1d0) then
          ux = 0d0
          uz = 0d0
          mden = 0.5d0 * piinv * GGinv * (cloud_asound/r)**2
          mdenf = mden
          if(time.eq.0d0) goto 100
          outangle = 0.5519d0*(time/tcloudm)**0.26d0
          if(outangle.gt.pi/2) outangle = pi/2
          bgden = 3d3 * mp * global_mugas
          aout = 0.191d0 / (time/tcloudm)**3
          if(aout.lt.(0.191d0)) aout = 0.191d0
          if(r*cos(th).gt.aout*au*(r*sin(th)/au)**1.5d0) then
              mdenf = bgden
              outflow = 1
          endif
  100     continue
          return
      endif
      tau = cloud_omega * time
      tau2 = tau * tau
      costh = cos(th)
      sinth = sin(th)
      legendre = 0.5d0 * (3d0*costh**2-1d0)
c
c     If x > 2e-3, interpolate the TSC84 solution. If x < 1e-3, use the
c     asymptotic values. For intermediate values of x, make a transition
c     from the interpolated to the asymptotic values. The interpolation
c     is done linearly or logarithmically, depending on the variable and
c     the value of x.
c
      if(x.gt.1d-3) then
          call hunt(tscxlog,NTSCLOG,xlog,j)
          if((j.le.0).or.(j.ge.NTSCLOG)) then
              write(*,*) j,NTSCLOG,xlog
              stop 8380
          endif
          eps = (xlog-tscxlog(j)) / (tscxlog(j+1)-tscxlog(j))
          epsn = 1.d0 - eps
          if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 8381
          v0 = epsn * tscv0(j) + eps * tscv0(j+1)
          wq = -exp(epsn * tscwqlog(j) + eps * tscwqlog(j+1))
          a0 = exp(epsn * tsca0log(j) + eps * tsca0log(j+1))
          am = exp(epsn * tscamlog(j) + eps * tscamlog(j+1))
          if(x.lt.3d-1) then
              vm = exp(epsn * tscvmlog(j) + eps * tscvmlog(j+1))
              vq = -exp(epsn * tscvqlog(j) + eps * tscvqlog(j+1))
              aq = -exp(epsn * tscaqlog(j) + eps * tscaqlog(j+1))
          else
              call hunt(tscxlin,NTSCLIN,x,j)
              if((j.le.0).or.(j.ge.NTSCLIN)) stop 8382
              eps = (x-tscxlin(j)) / (tscxlin(j+1)-tscxlin(j))
              epsn = 1.d0 - eps
              if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 8383
              vm = epsn * tscvmlin(j) + eps * tscvmlin(j+1)
              vq = epsn * tscvqlin(j) + eps * tscvqlin(j+1)
              aq = epsn * tscaqlin(j) + eps * tscaqlin(j+1)
          endif
          if(x.lt.2d-3) then
              f = 1d3 * (x-1d-3)
              fn = 1d0 - f
              xinv = 1d0 / x
              xinvsq = sqrt(xinv)
              xinv1p5 = xinv * xinvsq
              xinv2p5 = xinv * xinv1p5
              v0 = f * v0 - fn * m0twosq * xinvsq
              vm = f * vm + fn * 0.5d0 * m0half3p5o3 * xinv1p5
              vq = f * vq - fn * 0.5d0 * m0half3p5o3 * xinv1p5
              wq = f * wq - fn * m0half3p5o3 * xinv1p5
              a0 = f * a0 + fn * m0halfsq * xinv1p5
              am = f * am - fn * 0.25d0 * m0half3p5o3 * xinv2p5
              aq = f * aq - fn * 13d0*0.25d0 * m0half3p5o3 * xinv2p5
          endif
      else
          xinv = 1d0 / x
          xinvsq = sqrt(xinv)
          xinv1p5 = xinv * xinvsq
          xinv2p5 = xinv * xinv1p5
          v0 = -m0twosq * xinvsq
          vm = 0.5d0 * m0half3p5o3 * xinv1p5
          vq = -0.5d0 * m0half3p5o3 * xinv1p5
          wq = -m0half3p5o3 * xinv1p5
          a0 = m0halfsq * xinv1p5
          am = -0.25d0 * m0half3p5o3 * xinv2p5
          aq = -13d0*0.25d0 * m0half3p5o3 * xinv2p5
      endif
c
c     Convert the dimensionless velocities and density to dimensional
c     values. For small x, the TSC84 solution doesn't work for the
c     radial velocity, and we use the S77 solution instead.
c
      if(x.gt.1d-3) then
          ur = cloud_asound * (v0 + tau2 * (vm + vq*legendre))
          uth = -1.5d0 * cloud_asound * (tau2 * wq * sin(2d0*th))
      else
          ur = cloud_asound * v0
          uth = -1.5d0 * cloud_asound * (tau2 * wq * sin(2d0*th))
      endif
      gamma = tau2 * (0.5d0*(x**2*a0*(x-v0))*sinth)**2
      uphi = (cloud_asound/(x*tau*sinth)) * gamma
      ux = uth * costh + ur * sinth
      uz = -uth * sinth + ur * costh
      
      mden = (a0 + tau2 * (am + aq*legendre)) *
     $       0.25d0 * piinv * GGinv * timeinv**2
      mdenf = mden
      if(w.eq.1) then
          write(*,*) mden,a0,tau2,am,aq,legendre
      endif
      if(time.eq.0d0) goto 101
      outangle = 0.5519d0*(time/tcloudm)**0.26d0
      if(outangle.gt.pi/2) outangle = pi/2
      bgden = 3d3 * mp * global_mugas
      aout = 0.191d0 / (time/tcloudm)**3
      if(aout.lt.(0.191d0)) aout = 0.191d0
c     sept13
      if(r*costh.gt.aout*au*(r*sinth/au)**1.5d0) then
c          write(*,*) outflow, r*costh/au, 'tsc solution'
c          write(*,*) aout*(r*sinth/au)**(1.5d0)
          mdenf = bgden
          outflow = 1
      endif
  101 continue
c
      return
      end
c
c
c
c     ------------------------------------------------------------------
c     This subroutine computes the radial and angular velocity and the
c     density of an infalling parcel in the CM81 regime.
c     ------------------------------------------------------------------
      subroutine get_cmvelden(r,th,time,costh0,rinv,w)
c
      implicit none
      integer w,number_invalid
      doubleprecision r,th,time,costh0,rinv,mtot,costh,legendre,zeta
      doubleprecision sinth,outangle,aout,rcentr,mdotcm
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
c
      if(time.gt.tcloudm) then
          rcentr = 0.0625d0 * cloud_asound * (m0*tcloudm)**3 *
     $             cloud_omega**2
      else
          rcentr = 0.0625d0 * cloud_asound * (m0*time)**3 *
     $             cloud_omega**2
      endif
      mdotcm = m0 * cloud_asound**3 * GGinv
      mtot = mdotcm * time
      costh = cos(th)
      sinth = sin(th)
      legendre = 0.5d0 * (3d0*costh0**2-1d0)
      zeta = rcentr * rinv
      if(costh0.eq.0d0) then
          ur = -sqrt(GG*global_mstar*rinv)
          uth = 0d0
          uphi = -ur
          ux = ur
          uz = 0d0
          mden = -mdotcm / (4d0*pi*r**2*ur*(2d0*zeta))
      else
          ur = -sqrt((GG*global_mstar*rinv) * (1d0+costh/costh0))
          uth = -ur * (costh0-costh) / sqrt(1d0-costh**2)
          uphi = sqrt((GG*global_mstar*rinv) * (1d0-costh/costh0) *
     $                (1d0-costh0**2)/(1d0-costh**2))
          ux = uth * costh + ur * sinth
          uz = -uth * sinth + ur * costh
          mden = -mdotcm / (4d0 * pi * r**2 * ur * (2d0*zeta))
      endif
      mdenf = mden
      if(time.eq.0d0) goto 100
      outangle = 0.5519d0*(time/tcloudm)**0.26d0
      if(outangle.gt.pi/2) outangle = pi/2
      bgden = 3d3 * mp * global_mugas
      aout = 0.191d0 / (time/tcloudm)**3
      if(aout.lt.(0.191d0)) aout = 0.191d0
c     sept13
      if(r*costh.gt.aout*au*(r*sinth/au)**1.5d0) then
c          write(*,*) outflow, r*costh/au, 'cm solution'
c          write(*,*) r/au, aout*(r*sinth/au)**(1.5d0)
c          write(*,*) time/year, tcloudm/year, time/tcloudm
c          write(*,*) r*costh.gt.aout*au*(r*sinth/au)**(1.5d0)
          mdenf = bgden
          outflow = 1
      endif
  100 continue
c
      return
      end
c
c
c
c     ------------------------------------------------------------------
c     This subroutine computes the radial and angular velocity of an
c     infalling parcel in the disk regime.
c     ------------------------------------------------------------------
      subroutine get_diskvelden(x,z,time,jtime,eps,epsn,dovel)
c
      implicit none
      integer dovel,jtime,nlayers,iz,niz
      parameter(nlayers=1000)
      doubleprecision r,th,time,eps,epsn,scale,costh,sinth,x,xinv,dtemp
      doubleprecision rhom,surf,temp,dsigma,mdotstar,tempnew,mgrav
      doubleprecision mlow,mlayer,zz1(2*nlayers),zz1sq(2*nlayers),z,znew
      doubleprecision dzdt(2*nlayers),sum1,sum2,scaleinv,scalesqinv
      doubleprecision hlayer,zohsq,dt,surfnew,scalenew,mstarnew,erfi,dum
      doubleprecision lnx
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
c
c     Compute the disk radius.
c
      rdisk = 0.0625d0 * cloud_asound * (m0*time)**3 * cloud_omega**2
      if(time.gt.tcloudm) then
          rdisk = 0.0625d0 * cloud_asound * (m0*tcloudm)**3 *
     $            cloud_omega**2
      endif
      rin = 0.1d0 * au
      if(rin.lt.global_rstar) rin = global_rstar
      if(rdisk.lt.rin) then
          rdisk = rin
      endif
c
c     If the horizontal position falls outside of the grid, return zero
c     velocities.
c
      if(x.lt.grid_r(1).or.x.gt.grid_r(de_ngrid)) then
          mden = 0d0
          scale2 = 0d0
          surface = 0d0
          ux = 0d0
          uz = 0d0
          return
      endif
c
c     Compute the inverse of the horizontal position.
c
      xinv = 1d0 / x
      lnx = log(x)
c
c     Locate the current position on the disk evolution grid.
c
      jgrid = 1
      call hunt(grid_lnr,de_ngrid,lnx,jgrid)
      if(jgrid.ge.de_ngrid) then
          write(*,500) jgrid,de_ngrid,x,grid_r(1),grid_r(de_ngrid),dovel
  500     format(I5,I5,3(1PE13.5),I5)
          stop 8300
      else if(jgrid.le.0) then
          mden = -1d0
          scale2 = 0d0
          surface = 0d0
          ux = 0d0
          uz = 0d0
          return
      endif
      xeps = (lnx-grid_lnr(jgrid)) / (grid_lnr(jgrid+1)-grid_lnr(jgrid))
      xepsn = 1d0 - xeps
      if((xeps.lt.0.d0).or.(xeps.gt.1.d0)) stop 8301
c      write(*,*) jtime,eps,epsn,jgrid,xeps,xepsn
c
c     Compute the disk density.
c
      temp  = epsn * (xepsn * de_temp(jtime,jgrid) +
     $                xeps  * de_temp(jtime,jgrid+1)) +
     $        eps  * (xepsn * de_temp(jtime+1,jgrid) +
     $                xeps  * de_temp(jtime+1,jgrid+1))
      dktemp = temp
      mgrav = epsn * (xepsn * de_mgrav(jtime,jgrid) +
     $                xeps  * de_mgrav(jtime,jgrid+1)) +
     $        eps  * (xepsn * de_mgrav(jtime+1,jgrid) +
     $                xeps  * de_mgrav(jtime+1,jgrid+1))
      scale = sqrt(kk*temp*(x**3)*mpeffinv*GGinv/mgrav)
      scale2 = scale
      if(mdisk.gt.0d0) then
          surf = exp(epsn * (xepsn * de_lnsigma(jtime,jgrid) +
     $                       xeps  * de_lnsigma(jtime,jgrid+1)) +
     $               eps  * (xepsn * de_lnsigma(jtime+1,jgrid) +
     $                       xeps  * de_lnsigma(jtime+1,jgrid+1)))
          surf2 = surf
          rhom = surf / (scale * sqrt(2d0*pi))
          mden = rhom * exp(-0.5d0*(z/scale)**2)
      else
          mden = 0d0
      endif
c      write(*,47) x,z,time,temp,mgrav,scale,surf,mden
c   47 format(8(1PE13.5))
c
c     Compute velocity.
c
      if(dovel.eq.1.and.mdisk.gt.0d0) then
          ux = epsn * (xepsn * de_vr(jtime,jgrid) +
     $                 xeps  * de_vr(jtime,jgrid+1)) +
     $         eps  * (xepsn * de_vr(jtime+1,jgrid) +
     $                 xeps  * de_vr(jtime+1,jgrid+1))
          if(1.eq.1) then
          dt = de_time(jtime+1)-de_time(jtime)
          dsigma = exp(xepsn * (de_lnsigma(jtime+1,jgrid)-
     $                          de_lnsigma(jtime,jgrid)) +
     $                 xeps  * (de_lnsigma(jtime+1,jgrid+1)-
     $                          de_lnsigma(jtime,jgrid+1)))
          dtemp  = xepsn * (de_temp(jtime+1,jgrid)-
     $                      de_temp(jtime,jgrid)) +
     $             xeps  * (de_temp(jtime+1,jgrid+1)-
     $                      de_temp(jtime,jgrid+1))
          mdotstar = epsn * (de_mdot_c2s(jtime)+de_mdot_d2s(jtime)) +
     $               eps  * (de_mdot_c2s(jtime+1)+de_mdot_d2s(jtime+1))
          surfnew = surf + dsigma
          tempnew = temp + dtemp
          mstarnew = mstar + mdotstar * dt
          scalenew = sqrt(kk*tempnew*(x**3)*mpeffinv*GGinv/mstarnew)
          znew = scalenew * sqrt(2d0) *
     $           erfi(erf(z/(scale*sqrt(2d0)))*surf/surfnew)
          uz = (znew-z) / dt
          endif
      endif
c
c     Compute the location of the envelope/disk boundary.
c
      surface = epsn * (xepsn * de_surface(jtime,jgrid) +
     $                  xeps  * de_surface(jtime,jgrid+1)) +
     $          eps  * (xepsn * de_surface(jtime+1,jgrid) +
     $                  xeps  * de_surface(jtime+1,jgrid+1))
      global_loc = 'disk'
c
      return
      end
c
c
c
c     ------------------------------------------------------------------
c     This subroutine computes the radial and angular velocity and the
c     density of an infalling parcel in the CM81 regime.
c     ------------------------------------------------------------------
      subroutine get_vertden(x,z,r,th,t,costh0,w1)
c
      implicit none
      integer w1
      doubleprecision x,z,r,th,t,costh0
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
c
      integer number_invalid,i,nroot,w
      doubleprecision costh,rrc,check,mden81,mden84,ux81,uz81,uphi81
      doubleprecision mdenf81,mdenf84,prs81,prs84,mu0,xx,tau,tau2,eps
      doubleprecision xcub(3),acub(4),ux84,uz84,uphi84
      doubleprecision ur81,uth81,ur84,uth84
c     Added by Daniel, Jun 8, 2011
c     Need to scale the boundary conditions of different solution
c     such that the it scales with rotation
      doubleprecision scalexx
c     Added by Daniel, Jun 11, 2011
c     Need variables to define the transition boundary between
c     TSC and CM solutions
      doubleprecision upperbnd, lowerbnd,norm1,norm2
      doubleprecision offset
      upperbnd = 1d1
      lowerbnd = 1d-2
      offset = 2.0d0
      norm1 = 0.3333333d0
c      norm2 = (1.0d0)/(-1.0d0*log10(lowerbnd))
      norm2 = 0.5d0
c
c     If t = 0, return the TSC84 densities.
c
      if(t.eq.0d0) then
          call get_tscvelden(r,th,t,1d0,0)
          prs = mden * (ux**2+uz**2)
          return
      endif
c
c     Initialize some variables.
c
      xx = r / (cloud_asound * t)
      tau = cloud_omega * t
      tau2 = tau * tau
c     Added by Daniel: Jun 8, 2011
c     need to scale the xx for the transition between solution
c     this depends on the rc -> smaller rotation lead to smaller rc
c     hence this means that tsc84 solution is valid into smaller scale
      if (cloud_omega.gt.8d-14) then
        scalexx = 1d0
      elseif (cloud_omega.lt.1.2d-14) then
        scalexx = 1d-1
      else
        scalexx = (cloud_omega/1d-13)*(cloud_omega/1d-13)
c         scalexx = 1d0
      endif
      
c     
c     Nov 18, 2013
c     stupid code: cylinder symmetry assured
c      
!       if (z.lt.0) then
!         z = abs(z)
!       else if ((z/au).lt.1d-2) then
!         write(*,*) z/au, z
!         z = abs(0.1*au)
!       endif
      
      w = w1
      if(w.eq.1) then 
        write(*,*) 'cloud_asound,omega,t',cloud_asound,cloud_omega,t
        write(*,*) 'xx,tau,tau2',xx,tau,tau2
      endif
        
c
c     If costh0 < 0, it is a signal to compute costh0.
c
      if(xx.gt.upperbnd*tau2*scalexx) goto 10
      if(costh0.lt.0d0) then
          costh = z / r
          rrc = r / ml_rcentr
          if(z.eq.0d0) then
              if(rrc.lt.1d0) then
                  mu0 = sqrt(1d0-rrc)
              else
                  mu0 = 0d0
              endif
              goto 210
          endif
          if(rrc.eq.1d0) then
              mu0 = costh**f1o3
          else if(costh.eq.0d0) then
              if(rrc.lt.1d0) then
                  mu0 = sqrt(1d0-rrc)
              else
                  mu0 = 0d0
              endif
          else
              acub(4) = 1d0
              acub(3) = 0d0
              acub(2) = rrc-1d0
              acub(1) = -rrc*costh
              call cbcrt(acub,xcub)
              mu0 = xcub(1)
          endif
          if(costh.gt.0d0) then
              check = mu0 * (1d0-(1d0-mu0**2)/rrc)
              if(abs(1d0-check/costh).gt.1d-5) then
                  write(*,*) x/au, z/au, r/au, t/year
                  write(*,3342) r,ml_rcentr,rrc,costh,mu0,check,
     $                          check/costh,1d0-check/costh
 3342             format(8(1PE13.5))
                  stop 33423
              endif
          endif
      else
          mu0 = costh0
      endif
      global_mu0 = mu0
c
c     Get the CM81 density.
c
  210 continue
      call get_cmvelden(r,th,t,mu0,1d0/r,0)
      mden81 = mden
      mdenf81 = mdenf
      prs81 = mden * (ux**2+uz**2)
      ux81 = ux
      uz81 = uz
      ur81 = ur
      uth81 = uth
      uphi81 = uphi
      if(w.eq.1) write(*,*) 'nden81,ux81,uz81',mden81*mpeffinv,ux81,uz81
c
c     Get the TSC84 density.
c
      if(xx.lt.lowerbnd*tau2) goto 11
   10 continue
      if(t.gt.0.97d12.and.t.lt.0.98d12.and.x.gt.6.06d12.and.
     $                x.lt.6.07d12.and.z.gt.4.5d12.and.z.lt.7.1d12) then
          call get_tscvelden(r,th,t,1d0/t,0)
      else
          call get_tscvelden(r,th,t,1d0/t,0)
      endif
      mden84 = mden
      mdenf84 = mdenf
      prs84 = mden * (ux**2+uz**2)
      ux84 = ux
      uz84 = uz
      ur84 = ur
      uth84 = uth
      uphi84 = uphi
      if(w.eq.1) write(*,*) 'nden84,ux84,uz84',mden84*mpeffinv,ux84,uz84
c
c     Take the CM81 or the TSC84 density.
c
   11 continue
      eps = -1d0
c      write(*,101) t/year,x/au,z/au,xx,tau2,mden81,mden84
c  101 format(7(1PE13.5))
c If the xx << tau2 -> cm81
c NOTE: (jun 8, 2011)
c Older version has smaller range of smoothing instead of 
c xx < 0.01 tau2 and xx > 100 tau2
c
      if(xx.lt.lowerbnd*tau2) then
          mden = mden81
          prs = prs81
          ux = ux81
          uz = uz81
          ur = ur81
          uth = uth81
          uphi = uphi81
          global_loc = 'cm__'
c If the xx >> tau2 -> tsc84
      elseif(xx.gt.upperbnd*tau2*scalexx) then
          mden = mden84
          prs = prs84
          ux = ux84
          uz = uz84
          ur = ur84
          uth = uth84
          uphi = uphi84
          global_loc = 'tsc_'
      else
          if(mden81.lt.1d-50) then
              if(mden84.lt.1d-50) then
                  mden = 1d-50
                  prs = 1d-50
                 if (xx.lt.tau2) then
                      ux = ux81
                      uz = uz81
                      ur = ur81
                      uth = uth81
                      uphi = uphi81
                      global_loc = 'cm__'
                  else
                      ux = ux84
                      uz = uz84
                      ur = ur84
                      uth = uth84
                      uphi = uphi84
                      global_loc = 'tsc_'
                  endif
              else
                  mden = mden84
                  prs = prs84
                  ux = ux84
                  uz = uz84
                  ur = ur84
                  uth = uth84
                  uphi = uphi84
                  global_loc = 'tsc_'
              endif
          elseif(mden84.lt.1d-50) then
              mden = mden81
              prs = prs81
              ux = ux81
              uz = uz81
              ur = ur81
              uth = uth81
              uphi = uphi81
              global_loc = 'cm__'
          else
c              eps = (0.25d0) * (log10(xx/tau2)+2d0)
              if (cloud_omega.gt.8d-14) then
                 eps = (norm1) * (log10(xx/tau2) + offset)
              elseif (cloud_omega.lt.1.2d-14) then
                  eps = (norm2) * (log10(xx/tau2) + offset)
c              else
c                  eps = (1d0/(2d0-log10(1d2 * scalexx))) * 
c     $                  (log10(xx/tau2)+1d0)
              endif

              if(eps.lt.0d0.or.eps.gt.1d0) then
                  write(*,100) t/year,x/au,z/au,xx,tau2,mden81,mden84,
     $                         eps,scalexx
  100             format(8(1PE13.5))
                  stop 5000
              endif
              mden = eps * mden84 + (1d0-eps) * mden81
              prs = eps * prs84 + (1d0-eps) * prs81
              ux = eps * ux84 + (1d0-eps) * ux81
              uz = eps * uz84 + (1d0-eps) * uz81
              ur = eps * ur84 + (1d0-eps) * ur81
              uth = eps * uth84 + (1d0-eps) * uth81
              uphi = eps * uphi84 + (1d0-eps) * uphi81
              if(xx.lt.tau2) then
                  global_loc = 'cm__'
              else
                  global_loc = 'tsc_'
              endif
          endif
      endif
c     Need to check if mden is greater than 1d25

      if (mden.gt.1d25) then
        
        write(*,*) x/au,z/au,r/au,xx,tau2
        write(*,*) global_loc, mden,prs,ux,uz
        stop 50002

      endif

c
c     Add the thermal pressure (10 K) to the ram pressure.
c
      prs = max(prs,mden*mpeffinv*kk*10d0)
c      prs = prs + mden * mpeffinv * kk * 10d0
c
c     Get the density with correction for the outflow.
c
      if(xx.lt.lowerbnd*tau2) then
          mdenf = mdenf81
      elseif(xx.gt.upperbnd*tau2*scalexx) then
          mdenf = mdenf84
      else
          if(mdenf81.lt.1d-50) then
              if(mdenf84.lt.1d-50) then
                  mdenf = 1d-50
              else
                  mdenf = mdenf84
              endif
          elseif(mdenf84.lt.1d-50) then
              mdenf = mdenf81
          else
c              eps = 0.25d0 * (log10(xx/tau2)+2d0)
              if (cloud_omega.gt.8d-14) then
                  eps = (norm1) * (log10(xx/tau2) + offset)
              elseif (cloud_omega.lt.1.2d-14) then
                  eps = (norm2) * (log10(xx/tau2) + offset)
c              else
c                  eps = (1d0/(2d0-log10(1d2 * scalexx))) * 
c     $                  (log10(xx/tau2)+1d0)
              endif

              if(eps.lt.0d0.or.eps.gt.1d0) then
                  write(*,100) t/year,x/au,z/au,xx,tau2,mdenf81,mdenf84,
     $                         eps
                  stop 5001
              endif
              mdenf = eps * mdenf84 + (1d0-eps) * mdenf81
          endif
      endif

c     Added Jun 15, 2011

      if (t.gt.tcloudm) then
        
        if(mdenf81.lt.1d-50) then
          if(mdenf84.lt.1d-50) then
            mdenf = 1d-50
          else
            mdenf = mdenf84
          endif
        elseif(mdenf84.lt.1d-50) then
          mdenf = mdenf81
          
          if(mdenf.gt.8d3) then
c            write(*,*) 'I am here...'
            mdenf = (1.0d6/mpeffinv)*(r/(1d3*au))**(-1.5d0)
            write(*,*) mdenf81, mdenf84, mdenf
          endif
        
        endif
        
      endif

      if (number_invalid(mden).ne.0 .and. 
     $                number_invalid(mden81).ne.0) then
        mden = mden84
        mdenf = mden84
        prs = prs84
      endif
 
c
      if(mden.le.0d0.or.prs.le.0d0.or.number_invalid(mden).ne.0.or.
     $   mden.gt.1d25) then
          write(*,110) t/year,x/au,z/au,xx,tau2,mden81,mdenf81,mden84,
     $                 mdenf84,mden,mdenf,eps,prs81,prs84,prs
  110     format(15(1PE13.5))
          stop 1614
      endif
c
      return
      end
c
c
c
c     ------------------------------------------------------------------
c     This function takes the temperature at time t and position r,theta
c     from a grid of RADMC results with a three-dimensional
c     interpolation.
c     ------------------------------------------------------------------
      function get_dust_temp(t,r,th)
c
      implicit none
      doubleprecision t,r,th,x,z,epst,epsx,epsz,get_dust_temp
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
c
c     If time_tdust(1) is negative, this means no temperatures were read
c     in from RADMC, and we simply return the pre-collapse temperature
c     of 10 K.
c
      if(time_tdust(1).lt.0d0) then
          get_dust_temp = zerotemp
          return
      endif
c
c     If we're still before the end of the first hydrostatic core
c     [assumed to be time_tdust(1)], return the base temperature.
c
      if(t.lt.time_tdust(1)) then
          get_dust_temp = zerotemp
          return
      endif
c
      x = r * sin(th)
      z = r * cos(th)
c
c     Find the right point on the time grid.
c
      call hunt(time_tdust,nr_time_tdust,t,jt)
      if(jt.le.0) then
          jt = 1
          epst = 0d0
      else
          epst = (t-time_tdust(jt)) * dtime_tdust_inv(jt)
      endif
c
c     If our time exceeds the maximum time from RADMC, we will use the
c     results from that maximum time. Warn the user once about this.
c
      if(jt.ge.nr_time_tdust) then
         if(warning.eq.0) then
            write(*,*) t/year,time_tdust(1)/year
            write(*,*) jt,nr_time_tdust
            write(*,*) 'Passing the maximum time for which ',
     $                 'temperatures are defined.'
            warning = 1
         endif
         jt = nr_time_tdust-1
         epst = 1d0
      endif
      if((epst.lt.0.d0).or.(epst.gt.1.d0)) stop 8321
c
c     Find the right point on the spatial grids and do the
c     interpolation.
c
      call hunt(grid_r,grid_nr,x,jx)
      call hunt(grid_z,grid_nz,z,jz)
      if((jx.ge.grid_nr).or.(jz.ge.grid_nz)) then
          get_dust_temp = (1d0-epst) * dust_temp(grid_nr,grid_nz,jt) +
     $                  epst * dust_temp(grid_nr,grid_nz,jt+1)
          get_dust_temp = 10d0**get_dust_temp
      else
          if(jx.le.0) then
              jx = 1
              epsx = 0d0
          else
              epsx = (x-grid_r(jx)) / (grid_r(jx+1)-grid_r(jx))
          endif
          if(jz.le.0) then
              jz = 1
              epsz = 0d0
          else
              epsz = (z-grid_z(jz)) / (grid_z(jz+1)-grid_z(jz))
          endif
          if((epsx.lt.0.d0).or.(epsx.gt.1.d0).or.
     $        (epsz.lt.0.d0).or.(epsz.gt.1.d0)) stop 8323
          get_dust_temp = (1d0-epsz) * ((1d0-epsx) * 
     $               ((1d0-epst) * dust_temp(jx,jz,jt) +
     $                     epst  * dust_temp(jx,jz,jt+1)) +
     $        epsx * ((1d0-epst) * dust_temp(jx+1,jz,jt) +
     $                     epst  * dust_temp(jx+1,jz,jt+1))) +
     $        epsz * ((1d0-epsx) * 
     $               ((1d0-epst) * dust_temp(jx,jz+1,jt) +
     $                     epst  * dust_temp(jx,jz+1,jt+1)) +
     $        epsx * ((1d0-epst) * dust_temp(jx+1,jz+1,jt) + 
     $                     epst  * dust_temp(jx+1,jz+1,jt+1)))
          get_dust_temp = 10d0**get_dust_temp
      endif
c
c     Take a minimum value of 8 K to implicitly account for cosmic-ray
c     heating.
c
      if(get_dust_temp.lt.8d0) get_dust_temp = 8d0
c
      return
      end
c
c
c
c     ------------------------------------------------------------------
c     This function takes the UV field at time t and position r,theta
c     from a grid of RADMC results with a three-dimensional
c     interpolation.
c     ------------------------------------------------------------------
      subroutine get_rad2d(t,r,th)
c
      implicit none
      doubleprecision t,r,th,x,z,epst,epsx,epsz
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
#include "common_chem.h"
c
c     If time_tdust(1) is negative, we should not have gotten here.
c
      if(time_tdust(1).lt.0d0) then
          write(*,*) 'time_tdust(1) cannot be negative in get_rad2d.'
          stop
      endif
c
c     If we're still before the end of the first hydrostatic core
c     [assumed to be time_tdust(1)], return a zero radiation field.
c
      if(t.lt.time_tdust(1)) then
          global_g0 = 1d-90
          global_tau = 1d9
          return
      endif
c
      x = r * sin(th)
      z = r * cos(th)
c
c     Find the right point on the time grid.
c
      call hunt(time_tdust,nr_time_tdust,t,jt)
      if(jt.le.0) then
          jt = 1
          epst = 0d0
      else
          epst = (t-time_tdust(jt)) * dtime_tdust_inv(jt)
      endif
c
c     If our time exceeds the maximum time from RADMC, we will use the
c     results from that maximum time. Warn the user once about this.
c
      if(jt.ge.nr_time_tdust) then
         if(warning2.eq.0) then
            write(*,*) 'Passing the maximum time for which ',
     $                 '2D radiation fields are defined.'
            warning2 = 1
         endif
         jt = nr_time_tdust-1
         epst = 1d0
      endif
      if((epst.lt.0.d0).or.(epst.gt.1.d0)) stop 8621
c
c     Find the right point on the spatial grids and do the
c     interpolation.
c
      call hunt(grid_r,grid_nr,x,jx)
      call hunt(grid_z,grid_nz,z,jz)
      if((jx.ge.grid_nr).or.(jz.ge.grid_nz)) then
          global_g0 = (1d0-epst) * rad2d(grid_nr,grid_nz,jt) +
     $                epst * rad2d(grid_nr,grid_nz,jt+1)
          global_g0 = 10d0**global_g0
      else
          if(jx.le.0) then
              jx = 1
              epsx = 0d0
          else
              epsx = (x-grid_r(jx)) / (grid_r(jx+1)-grid_r(jx))
          endif
          if(jz.le.0) then
              jz = 1
              epsz = 0d0
          else
              epsz = (z-grid_z(jz)) / (grid_z(jz+1)-grid_z(jz))
          endif
          if((epsx.lt.0.d0).or.(epsx.gt.1.d0).or.
     $        (epsz.lt.0.d0).or.(epsz.gt.1.d0)) stop 8623
          global_g0 = (1d0-epsz) * ((1d0-epsx) * 
     $               ((1d0-epst) * rad2d(jx,jz,jt) +
     $                     epst  * rad2d(jx,jz,jt+1)) +
     $        epsx * ((1d0-epst) * rad2d(jx+1,jz,jt) +
     $                     epst  * rad2d(jx+1,jz,jt+1))) +
     $        epsz * ((1d0-epsx) * 
     $               ((1d0-epst) * rad2d(jx,jz+1,jt) +
     $                     epst  * rad2d(jx,jz+1,jt+1)) +
     $        epsx * ((1d0-epst) * rad2d(jx+1,jz+1,jt) + 
     $                     epst  * rad2d(jx+1,jz+1,jt+1)))
          global_g0 = 10d0**global_g0
      endif
      global_tau = 1d-90
c
      return
      end
c
c
c
c     ------------------------------------------------------------------
c     This function takes the H2 column density at time t and position
c     r,theta from a grid of RADMC results with a three-dimensional
c     interpolation.
c     ------------------------------------------------------------------
      subroutine get_cdh(t,r,th)
c
      implicit none
      doubleprecision t,r,th,x,z,epst,epsx,epsz
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
#include "common_chem.h"
c
c     If time_tdust(1) is negative, we should not have gotten here.
c
      if(time_tdust(1).lt.0d0) then
          write(*,*) 'time_tdust(1) cannot be negative in get_cdh.'
          stop
      endif
c
c     If we're still before the end of the first hydrostatic core
c     [assumed to be time_tdust(1)], return a column density of 1.
c
      if(t.lt.time_tdust(1)) then
          global_cdh = 1.59d23
          return
      endif
c
      x = r * sin(th)
      z = r * cos(th)
c
c     Find the right point on the time grid.
c
      call hunt(time_tdust,nr_time_tdust,t,jt)
      if(jt.le.0) then
          jt = 1
          epst = 0d0
      else
          epst = (t-time_tdust(jt)) * dtime_tdust_inv(jt)
      endif
c
c     If our time exceeds the maximum time from RADMC, we will use the
c     results from that maximum time. Warn the user once about this.
c
      if(jt.ge.nr_time_tdust) then
         if(warning3.eq.0) then
            write(*,*) 'Passing the maximum time for which ',
     $                 '2D H2 column densities are defined.'
            warning3 = 1
         endif
         jt = nr_time_tdust-1
         epst = 1d0
      endif
      if((epst.lt.0.d0).or.(epst.gt.1.d0)) stop 8621
c
c     Find the right point on the spatial grids and do the
c     interpolation.
c
      call hunt(grid_r,grid_nr,x,jx)
      call hunt(grid_z,grid_nz,z,jz)
      if((jx.ge.grid_nr).or.(jz.ge.grid_nz)) then
          global_cdh = (1d0-epst) * cdh(grid_nr,grid_nz,jt) +
     $                  epst * cdh(grid_nr,grid_nz,jt+1)
          global_cdh = 10d0**global_cdh
      else
          if(jx.le.0) then
              jx = 1
              epsx = 0d0
          else
              epsx = (x-grid_r(jx)) / (grid_r(jx+1)-grid_r(jx))
          endif
          if(jz.le.0) then
              jz = 1
              epsz = 0d0
          else
              epsz = (z-grid_z(jz)) / (grid_z(jz+1)-grid_z(jz))
          endif
          if((epsx.lt.0.d0).or.(epsx.gt.1.d0).or.
     $        (epsz.lt.0.d0).or.(epsz.gt.1.d0)) stop 8623
          global_cdh = (1d0-epsz) * ((1d0-epsx) * 
     $                 ((1d0-epst) * cdh(jx,jz,jt) +
     $                       epst  * cdh(jx,jz,jt+1)) +
     $          epsx * ((1d0-epst) * cdh(jx+1,jz,jt) +
     $                       epst  * cdh(jx+1,jz,jt+1))) +
     $          epsz * ((1d0-epsx) * 
     $                 ((1d0-epst) * cdh(jx,jz+1,jt) +
     $                       epst  * cdh(jx,jz+1,jt+1)) +
     $          epsx * ((1d0-epst) * cdh(jx+1,jz+1,jt) + 
     $                       epst  * cdh(jx+1,jz+1,jt+1)))
          global_cdh = 10d0**global_cdh
      endif
      if(global_cdh.lt.1d0) global_cdh = 1d0
c
      return
      end
c
c
c
c     ------------------------------------------------------------------
c     This subroutine creates an input file for the chemistry with the
c     temperature, density and time derivative of the density at a
c     number of time steps.
c     ------------------------------------------------------------------
      subroutine write_chem_input(ip,do_tau,do_cdh)
c
      implicit none
      integer ip,do_tau,do_cdh,i,nfco,i1,ir,iz,it,number_invalid,nextra
      doubleprecision tchem(NCHEM),tempchem(NCHEM)
      doubleprecision ndenchem(NCHEM),dndtgchem(NCHEM)
      doubleprecision rchem(NCHEM),thchem(NCHEM),get_dust_temp
      doubleprecision g0chem(NCHEM),tauchem(NCHEM),cdhchem(NCHEM)
      doubleprecision dn,dux,duz,dx,dz,dndx,dndz,duxdx,duzdz,dndtg(step)
      doubleprecision lnadndtg(step),dum,dt,dspline(NCHEM),d1(NCHEM)
      doubleprecision d2(NCHEM),dlnn,dlnndtg(step),dlnndtgchem(NCHEM)
      doubleprecision temptemp(step),disttol,gcode1(step),gcode2(step)
      doubleprecision gcode3(step),g0temp(step),tautemp(step),lng0(step)
      doubleprecision cdhtemp(step),extime,extemp,excdh,exg0,extau
      character*15 filename
      parameter(disttol=3d-1)
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
#include "common_chem.h"
c
c     Get temperature, UV intensity and H2 column density at each time
c     step.
c
      do i=1,step
          temptemp(i) = get_dust_temp(tstore(i),rstore(i),thstore(i))
          if(do_tau.eq.1) then
              write(*,*) 'Error: Invalid value for do_tau:',do_tau
              stop
          elseif(do_tau.eq.2) then
              call get_rad2d(tstore(i),rstore(i),thstore(i))
              g0temp(i) = global_g0
              tautemp(i) = global_tau
          else
              g0temp(i) = 1d-90
              tautemp(i) = 1d9
          endif
          lng0(i) = log(g0temp(i))
          if(do_cdh.eq.1) then
              call get_cdh(tstore(i),rstore(i),thstore(i))
              cdhtemp(i) = global_cdh
          else
              cdhtemp(i) = 1d0
          endif
      enddo
c
c     Compute dn/dt.
c
      do i=1,step-1
          dn = nstore(i+1)-nstore(i)
          dt = tstore(i+1)-tstore(i)
          dndtg(i) = dn/dt
          if(number_invalid(dndtg(i)).ne.0) then
              write(*,*) 'dndtg = ',dndtg(i)
              write(*,*) step,dn,dt
              write(*,*) nstore(i+1),nstore(i),nstore(i-1)
              write(*,*) tstore(i+1),tstore(i),tstore(i-1)
              write(*,*) 'at ip,i = ',ip,i
              stop 15231
          endif
          lnadndtg(i) = log(abs(dndtg(i)))
      enddo
      dndtg(step) = dndtg(step-1)
      lnadndtg(step) = lnadndtg(step-1)
c
c     Reduce the number of points.
c
      call generalize(tstore,temptemp,step,gcode1)
      call generalize(tstore,lnadndtg,step,gcode2)
      call generalize(tstore,lng0,step,gcode3)
      nfco = 0
      do i=1,step
c          if(gcode1(i).gt.disttol.or.gcode2(i).gt.disttol.or.
c     $                                        gcode3(i).gt.disttol) then
              nfco=nfco+1
              if(nfco.gt.NCHEM) then
                  write(*,*) 'nfco > NCHEM for ip =',ip,nfco,NCHEM
                  stop 88665
              endif
              tchem(nfco) = tstore(i)
              rchem(nfco) = rstore(i)
              thchem(nfco) = thstore(i)
              if (nstore(i).lt.1d0) then
                nstore(i) = 1d0
              endif
              ndenchem(nfco) = nstore(i)
              tempchem(nfco) = temptemp(i)
              dndtgchem(nfco) = dndtg(i)
              dlnndtgchem(nfco) = dlnndtg(i)
              g0chem(nfco) = g0temp(i)
              tauchem(nfco) = tautemp(i)
              cdhchem(nfco) = cdhtemp(i)
c          endif
      enddo
      dndtgchem(1) = 0d0
      dlnndtgchem(1) = 0d0
c
c     Make the filename and open the file for writing.
c
c      write(*,*) 'Writing doty file...'
      write(filename,119) 'doty.',ip-1,'.dat'
  119 format(A5,I6.6,A4)
      open(unit=1,file=filename)
c
c     Write the steps to the file.
c
      nextra = floor(0.999d0*tchem(1)/(1d4*year))
      write(1,120) nfco+nextra+1
  120 format(I4)
      write(1,122) 0d0,rstore(1),thstore(1),zerotemp,
     $             nstore(1)+1d-90,0d0,cdhchem(1),1d-90,1d9
c      write(1,122) precoll,rstore(1),thstore(1),zerotemp,
c     $             nstore(1)+1d-90,0d0,cdhchem(1),1d-90,1d9
      do i=1,nextra
          extime = i*1d4*year
          extemp = get_dust_temp(extime,rstore(1),thstore(1))
          if(do_cdh.eq.1) then
              call get_cdh(tstore(i),rstore(i),thstore(i))
              excdh = global_cdh
          else
              excdh = 1d0
          endif
          if(do_tau.eq.2) then
              call get_rad2d(tstore(i),rstore(i),thstore(i))
              exg0 = global_g0
              extau = global_tau
          else
              exg0 = 1d-90
              extau = 1d9
          endif
          
          write(1,122) precoll+extime,rstore(1),thstore(1),extemp,
     $                 nstore(1)+1d-90,0d0,excdh,exg0,extau
      enddo
      do i=1,nfco
          write(1,122) precoll+tchem(i),rchem(i),
     $                 thchem(i),tempchem(i),ndenchem(i)+1d-90,
     $                 dndtgchem(i),cdhchem(i),g0chem(i),tauchem(i)
      enddo
  122 format(1PE24.16,1PE13.5,0PF9.5,1PE13.5,2(1PE24.16),3(1PE10.3))
      close(1)
c
c     Use the nncoll variable to track across all parcels what the
c     largest number of data points in the chemistry file is.
c
      if(nfco+2.gt.nncoll) nncoll = nfco+nextra+2
c
c      stop 4477
      return
      end
c
c
c
c     ------------------------------------------------------------------
c     Added: DH, Sep 9 2011
c     This subroutine computes the evolution of the parcels at fixed
c     positions.
c     ------------------------------------------------------------------
      subroutine make_evolvepts(ip,dowritechem,do_tau,do_cdh)
c
      implicit none
      integer ip,dowritechem,do_tau,do_cdh
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
#include "common_traj.h"
c
      integer j,leniwd,lenwd,ipar,jtime
      parameter(leniwd=32)
      parameter(lenwd=62) !48
      integer itime,zmethod,i,nroot
      doubleprecision xz(2),t,tinv,t00,told,rpar,lnt,deeps,deepsn
      doubleprecision dr,dum,erfi,scaleold,tflow,rold,thold
      doubleprecision check,costh,eps,nden,nden0,r,ndenold,ndensave
      doubleprecision rcen,rrc,rsave,temp,theta,thsave,get_dust_temp
      doubleprecision ndenenv,ndendisk,prsenv,prsdisk
      doubleprecision ndenf,dummu0,xx,tau,tau2,nden81,nden84
      doubleprecision xcub(3),acub(4),rinv
c Added by DH Sep 15, 2011 -> scale the CM81 boundary
      doubleprecision scalexx
      character*4 loc
      external fparcel,jacobn
c     Added by Daniel, Sep 15, 2011
c     Need variables to define the transition boundary between
c     TSC and CM solutions
      doubleprecision upperbnd, lowerbnd,norm1,norm2
      doubleprecision offset
      upperbnd = 1d1
      lowerbnd = 1d-2
      offset = 2.0d0
      norm1 = 0.3333333d0
c      norm2 = (1.0d0)/(-1.0d0*log10(lowerbnd))
      norm2 = 0.5d0

c     Added by Daniel: Sep 15, 2011
c     need to scale the xx for the transition between solution
c     this depends on the rc -> smaller rotation lead to smaller rc
c     hence this means that tsc84 solution is valid into smaller scale
      
      if (cloud_omega.gt.8d-14) then
        scalexx = 1d0
      elseif (cloud_omega.lt.1.2d-14) then
        scalexx = 1d-1
      else
        scalexx = (cloud_omega/1d-13)*(cloud_omega/1d-13)
      endif
c
c     Initialize several parameters. The time is initialized to the
c     point where the expansion wave reaches the initial position
c     (r = c_s t).
c
      xz(1) = r0(ip)
      xz(2) = th0(ip)
      r0(ip) = sqrt(xz(1)**2+xz(2)**2)
      th0(ip) = atan2(xz(1),xz(2))
      t = t0(ip)
      if(t.lt.r0(ip)*cloud_asound_inv) then
          t = 1.00000001d0 * r0(ip) * cloud_asound_inv
      endif
      t00 = t
      rold = r0(ip)
      theta = th0(ip)
      thold = theta
      told = t
      nden0 = mpeffinv * 0.5d0 * piinv * GGinv *
     $        (cloud_asound/r0(ip))**2
      ndenold = nden0
      scaleold = -1d0
      tcm81 = 1d90
      tdisk = 1d90
c
c     Locate t00 on the diskevol time grid.
c
      itime = 1
c      write(*,*) t00, t00/year
      call hunt(de_time,de_ntime,t00,itime)
      if((itime.le.0).or.(itime.ge.de_ntime)) stop 8350
      temp = get_dust_temp(t00,r0(ip),th0(ip))
c
c     Write the initial values at all time steps up to itime.
c
      do j=1,itime
          if(nparc.lt.25000.or.j.eq.102.or.j.eq.104) then
              write(10,100) j,de_time(j),r0(ip),th0(ip),nden0,temp
          endif
      enddo
  100 format(I5,2(1PE16.8),0PF12.8,1PE16.8,1PE13.5)
  101 format(I5,4(1PE16.8),0PF12.8,6(1PE16.8),1X,A5,I2)
      itime = itime + 1
c
c     If the infall starts after the first hydrostatic core phase
c     [assumed to correspond to 2.01d4*year], add two steps to the
c     storage arrays to get the correct density and temperature profile.
c
      if(t.gt.2.01d4*year) then
          tstore(1) = 0.999d0*2.01d4*year
          rstore(1) = r0(ip)
          thstore(1) = th0(ip)
          nstore(1) = nden0
          if(write_traj.eq.1) then
              write(20,101) 0,0.999d0*2.01d4*year,xz(1),xz(2),r0(ip),
     $                      th0(ip),nden0,0d0,0d0,0d0,0d0,zerotemp,
     $                      'init',0
          endif
          tstore(2) = 2.01d4*year
          rstore(2) = r0(ip)
          thstore(2) = th0(ip)
          nstore(2) = nden0
          if(write_traj.eq.1) then
              temp = get_dust_temp(2.01d4*year,r0(ip),th0(ip))
              write(20,101) 1,2.01d4*year,xz(1),xz(2),r0(ip),th0(ip),
     $                      nden0,0d0,0d0,0d0,0d0,temp,'init',0
          endif
          step = 3
      else
          step = 1
      endif
      if(write_traj.eq.1) then
          temp = get_dust_temp(t00,r0(ip),th0(ip))
          write(20,101) step-1,t00,xz(1),xz(2),r0(ip),th0(ip),nden0,0d0,
     $                  0d0,0d0,0d0,temp,'init',0
      endif
c
c     If the parcel's initial angle is less than the maximum outflow
c     angle, compute the time at which the outflow angle equals the
c     parcel's initial angle.
c
      if(xz(2).gt.0.191d0*au*(xz(1)/au)**1.5d0) then
          tflow = (xz(2)/(0.191d0*au*(xz(1)/au)**1.5d0))**(-f1o3) *
     $            tcloudm
      else
          tflow = 10d0 * tmax
      endif
c
c     Need to find the location of the point: CM or TSC
c     or is it starting in the disk
c
      rcen = 0.0625d0*cloud_asound*((m0*t)**3)*(cloud_omega**2)
      if(r0(ip)/(cloud_asound*t00).lt.(cloud_omega*t00)**2) then
          costh = cos(th0(ip))
          rrc = r0(ip) / rcen
          if(rrc.eq.1d0) then
              mu0 = costh**f1o3
          else if(costh.eq.0d0) then
              if(rrc.lt.1d0) then
                  mu0 = sqrt(1d0-rrc)
              else
                  mu0 = 0d0
              endif
          else
              acub(4) = 1d0
              acub(3) = 0d0
              acub(2) = rrc-1d0
              acub(1) = -rrc*costh
              call cbcrt(acub,xcub)
              mu0 = xcub(1)
          endif
          if(costh.gt.0d0) then
              check = mu0 * (1d0-(1d0-mu0**2)/rrc)
              if(abs(1d0-check/costh).gt.1d-6) stop 3362
          endif
          tcm81 = t00
          loc = 'cm__'
      else
          loc = 'tsc_'
      endif
c
c     Let's start the evolution of the points from here
c     Need to find the density and the time step
c
      r = sqrt(xz(1)*xz(1)+xz(2)*xz(2))
      theta = atan2(xz(1),xz(2))
      rinv = 1d0/r

 1000 continue

c
c     Find the next time step
c     Increase the time and steps
c     Set the time step
c
      time_dt = max(de_time(2),t*time_dtrel)
c
c     Check if we are about to cross the save time
c
      if(t+time_dt.ge.de_time(itime)) then
          time_dt = de_time(itime)-t
      endif
c
c     Check if we are about to cross the time-limit
c
      if(t+time_dt.gt.time_end) then
          time_dt = time_end - t
      endif
c      write(*,*) t, time_dt, step,itime
      t = t + time_dt + 0.3*year
      tinv = 1d0/t
      lnt = log(t)
c
c     Need to hunt for the time
c
      jtime = 1
      call hunt(de_lntime,de_ntime,lnt,jtime)
      if( (jtime.le.1).or.(jtime.ge.de_ntime)) then
          write(*,*) itime, time_dt
          write(*,*) t, lnt,t/year,jtime, de_ntime
          stop 5132
      endif

c
c     Need to find the star and disk masses
c

      deeps = (lnt-de_lntime(jtime)) / 
     $     (de_lntime(jtime+1)-de_lntime(jtime))
      
      deepsn = 1d0-deeps
      
      if((deeps.lt.0.d0).or.(deeps.gt.1.d0)) stop 1542
      mstar = deepsn * de_mstar(jtime) + deeps * de_mstar(jtime+1)
      mdisk = deepsn * de_mdisk(jtime) + deeps * de_mdisk(jtime+1)
c
c     Need to find the location of the point: CM or TSC
c
      if(t.lt.tcm81) then
          xx = r * cloud_asound_inv * tinv
          if(xx.lt.(cloud_omega*t)**2) then
              tcm81 = t
              rdisk = 0.0625d0 * cloud_asound * (m0*t)**3 *
     $                cloud_omega**2
              rrc = r / rdisk
              costh = xz(2) * rinv
              if(rrc.eq.1d0) then
                  mu0 = costh**f1o3
              else if(costh.eq.0d0) then
                  if(rrc.lt.1d0) then
                      mu0 = sqrt(1d0-rrc)
                  else
                      mu0 = 0d0
                  endif
              else
                  acub(4) = 1d0
                  acub(3) = 0d0
                  acub(2) = rrc-1d0
                  acub(1) = -rrc*costh
                  call cbcrt(acub,xcub)
                  mu0 = xcub(1)
              endif
              if(dummu0.lt.-2d0) dummu0 = mu0
              if(costh.gt.0d0) then
                  check = mu0 * (1d0-(1d0-mu0**2)/rrc)
                  if(abs(1d0-check/costh).gt.1d-6) then
                      write(*,*) r,rdisk,costh,mu0,check
                      stop 10003
                  endif
              endif
          endif
      endif
c
c     Calculate the densities following the fparcel way 
c     The method of computing the velocity depends on the type of regime
c     the parcel is currently in. If we're in the TSC84 or CM81 regime,
c     we also need to get the disk density to determine if the
c     parcel has entered the disk. The 'velden' routines overwrite the
c     value of the global parameter mden, so we have to store the
c     densities in other variables.
c
      if(t.gt.tdisk) then
          call get_diskvelden(xz(1),xz(2),t,jtime,deeps,deepsn,1)
          ndendisk = mden * mpeffinv
      elseif(t.gt.tcm81) then
          call get_cmvelden(r,theta,t,mu0,rinv,0)
          ndenenv = mden * mpeffinv
          ndenf = mdenf * mpeffinv
          prsenv =  mden * (ux*ux + uz*uz)
          call get_diskvelden(xz(1),xz(2),t,jtime,deeps,deepsn,0)
          ndendisk = mden * mpeffinv
          prsdisk = ndendisk * kk * dktemp
      else
          call get_tscvelden(r,theta,t,tinv,0)
          ndenenv = mden * mpeffinv
          ndenf = mdenf * mpeffinv
          prsenv = mden * (ux*ux + uz*uz)
          call get_diskvelden(xz(1),xz(2),t,jtime,deeps,deepsn,0)
          ndendisk = mden * mpeffinv
          prsdisk = ndendisk * kk * dktemp
      endif

c
c    Check if we are in the disk or envelope
c

      if( t.lt.tdisk) then
         
         if(prsdisk.gt.prsenv) then
             tdisk = t
             if(t.lt.tcm81) tcm81 = t
             nden = ndendisk
         else
           
             nden = ndenenv
c     
c     Note: DH, Sep 15, 2011
c     if xx E [0.01 and scalexx]: Need to scale the densities
c
c
c     Initialize some variables.
c
             xx = r / (cloud_asound * t)
             tau = cloud_omega * t
             tau2 = tau * tau
     
             if( (xx.lt.upperbnd*tau2*scalexx).and.
     $          (xx.gt.lowerbnd*tau2)) then
                call get_cmvelden(r,theta,t,mu0,1d0/r,0)
                nden81 = mden * mpeffinv
                call get_tscvelden(r,theta,t,1d0/t,0)
                nden84 = mden * mpeffinv
          
                if (cloud_omega.gt.8d-14) then
                    eps = (norm1) * (log10(xx/tau2) + offset)
                elseif (cloud_omega.lt.1.2d-14) then
                    eps = (norm2) * (log10(xx/tau2) + offset)
                endif

                if(eps.lt.0d0.or.eps.gt.1d0) then
                    write(*,151) t/year,r,xx,tau2,nden81,nden84,
     $              eps,scalexx
  151               format(8(1PE13.5))
                    stop 5123
                endif
                nden = eps * nden84 + (1d0-eps) * nden81
              endif
          endif
      else
          nden = ndendisk
      endif
c      if (ip.gt.65) then
c          write(*,*) step, nden
c      endif

      nden = max(1d-50,nden)

c
c     If we crossed a diskevol time step, compute and write the position
c     corresponding to this step. Also check if we perhaps crossed two
c     or more diskevol steps at once.
c
 1030 if(t.ge.de_time(itime)) then
          eps = (de_time(itime)-told) / (t-told)
          if((eps.lt.0.d0).or.(eps.gt.1.d0)) then
              write(*,*) ip,eps,de_time(itime),told,t,itime
              write(*,*) de_time(de_ntime)/year,tcloudm/year
              stop 8341
          endif
          rsave = r
          thsave = theta
          ndensave = (1.d0-eps) * ndenold + eps * nden
          temp = get_dust_temp(t,r0(ip),th0(ip))
          if(nparc.lt.25000.or.itime.eq.102.or.itime.eq.104) then
              write(10,100) 
     $         itime,de_time(itime),rsave,thsave,ndensave,temp        
          endif
          tstore(step) = de_time(itime)
          rstore(step) = rsave
          thstore(step) = thsave
          nstore(step) = ndensave
c          if (ip.gt.65) then
c              write(*,*) step, tstore(step), nstore(step)
c              write(*,*)
c          endif
c
c    Check if these are 0s...
c
           if( (tstore(step).eq.0).or.(nstore(step).eq.0).or.
     #         (tstore(step).eq.tstore(step-1)).or.
     #         (nstore(step).eq.nstore(step-1)) ) then
               write(*,*) ip,step,itime
               write(*,*) t, told
               write(*,*) nden,ndenold
               write(*,*) tstore(step),tstore(step-1)
               write(*,*) nstore(step),nstore(step-1)
               write(*,*) t, t/year, de_time(itime)
               write(*,*) eps
               stop 1974
           endif
c
c     update the steps and time
c
           itime = itime + 1
           step = step +1
           if(itime.le.de_ntime) goto 1030
      endif
c
c     Store the results from this step.
c
      tstore(step) = t
      rstore(step) = r
      thstore(step) = theta
      nstore(step) = nden

c      if (ip.gt.65) then
c         write(*,*) step, tstore(step), nstore(step)
c         write(*,*)
c         if (ip.eq.68) stop 5132
c      endif

c
c     Did we reach the star or cross the maximum time? If not, save the
c     current values and do another loop. If we did, stop the iteration
c     for this parcel.
c
      if(t.gt.tmax) then
          step = step-1
          goto 1020
          
      endif
      told = t
      rold = r
      thold = theta
      ndenold = nden - (nden*1d-3)
c
c     Did we reach the max number of steps? If so, report the problem
c     and stop the iteration. If not, go to the next step.
c
      if(step.eq.NSTORE) then
          write(*,*) 'Max. number of steps reached. You may want to ',
     $               'recompile with a larger NSTORE.'
          write(*,*) ip,t
          stop 54123
      endif
      step = step+1
      goto 1000
 1020 continue
c
c     Write values for the remaining steps (radial distance of 1 cm
c     and angle of pi/4).
c
      if(itime.le.de_ntime) then
          do j=itime,de_ntime
              if(nparc.lt.25000.or.j.eq.102.or.j.eq.104) then
                  write(10,100) j,de_time(j),5d0,0.78539816d0,1d0,5d0
              endif
          enddo
      endif
      if(write_traj.eq.1) then
          write(20,101) -1,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,
     $                  0d0,0d0,0d0,'end_',0
      endif
c
c     Create an input file for the chemistry with the temperature,
c     density and time derivative of the density at a number of time
c     steps.
c
      if(dowritechem.eq.1) then
          call write_chem_input(ip,do_tau,do_cdh)
      endif
c
      return
      end
c
c
c
c     ------------------------------------------------------------------
c     This subroutine computes the infall trajectory for a single
c     parcel using the DVODE solver.
c     ------------------------------------------------------------------
      subroutine make_trajvode(ip,dowritechem,do_tau,do_cdh)
c
      implicit none
      integer ip,dowritechem,do_tau,do_cdh
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
#include "common_traj.h"
c
      integer j,itol,ntask,iopt,mf,neq,leniwd,lenwd,nstate,ipar,jtime
      parameter(leniwd=32)
      parameter(lenwd=62) !48
      integer iwork(leniwd),itime,zmethod,i,nroot
c Added by DH Sep 15, 2011 -> scale the CM81 boundary
      doubleprecision scalexx
      doubleprecision xz(2),rtol,atol,t,t00,told,rpar,lnt,deeps,deepsn
      doubleprecision work(lenwd),dr,dum,erfi,scaleold,tflow,rold,thold
      doubleprecision check,costh,eps,nden,nden0,r,ndenold,ndensave
      doubleprecision rcen,rrc,rsave,temp,theta,thsave,get_dust_temp
      doubleprecision xcub(3),acub(4),xx,tau,tau2,nden81,nden84
      doubleprecision dummu0, mdendisk, mdenenv, prsdisk,prsenv
      doubleprecision rhom,rinv,ztemp,tcloud,get_lstar
      character*4 loc
      external fparcel,jacobn
c     Added by Daniel, Sep 15, 2011
c     Need variables to define the transition boundary between
c     TSC and CM solutions
      doubleprecision upperbnd, lowerbnd,norm1,norm2
      doubleprecision offset
      upperbnd = 1d1
      lowerbnd = 1d-2
      offset = 2.0d0
      norm1 = 0.3333333d0
      norm2 = 0.5d0

!       if (ip > 5999) then
! c        write(*,*) 'I am here...'
!       endif

c     Added by Daniel: Sep 15, 2011
c     need to scale the xx for the transition between solution
c     this depends on the rc -> smaller rotation lead to smaller rc
c     hence this means that tsc84 solution is valid into smaller scale
      
      if (cloud_omega.gt.8d-14) then
        scalexx = 1d0
      elseif (cloud_omega.lt.1.2d-14) then
        scalexx = 1d-1
      else
        scalexx = (cloud_omega/1d-13)*(cloud_omega/1d-13)
      endif
c
c     Initialize several parameters. The time is initialized to the
c     point where the expansion wave reaches the initial position
c     (r = c_s t).
c
      xz(1) = r0(ip)
      xz(2) = th0(ip)
      r0(ip) = sqrt(xz(1)**2+xz(2)**2)
      th0(ip) = atan2(xz(1),xz(2))
      t = t0(ip)
      if(t.lt.r0(ip)*cloud_asound_inv) then
          t = 1.00000001d0 * r0(ip) * cloud_asound_inv
      endif
      t00 = t
      rold = r0(ip)
      theta = th0(ip)
      thold = theta
      told = t
      nden0 = mpeffinv * 0.5d0 * piinv * GGinv *
     $        (cloud_asound/r0(ip))**2
      ndenold = nden0
      scaleold = -1d0
      tcm81 = 1d90
      tdisk = 1d90
c
c     Get the stellar properties.
c
      rcloud = GG * cloud_mass / (2d0*cloud_asound**2)
      tcloud = rcloud / cloud_asound
      if(t.lt.tcloudm) then
          ml_mcentr = (0.5d0*m0) * cloud_mass * (t/tcloud)
          mdot = m0 * cloud_asound**3 * GGinv
          ml_rcentr = 0.0625d0 * cloud_asound * ((m0*t)**3) *
     $                (cloud_omega**2)
      else
          ml_mcentr = cloud_mass
          mdot = 0d0
          ml_rcentr = 0.0625d0 * cloud_asound * ((m0*tcloudm)**3) *
     $                (cloud_omega**2)
      endif
      time_time = t
      global_mstar = mstar
      call get_rstar()
      global_lstar = get_lstar()
c
c     Initialize constant parameters for DVODE.
c
      neq = 2
      itol = 1
      rtol = 1d-8
      ntask = 2
      work(5) = 0d0
      work(6) = 1d12
      work(7) = 0d0
      iwork(5) = 0
      iwork(6) = 30000
      iwork(7) = 0
      iopt = 1
      mf = 12 !22
      nstate = 1
c
c     If the parcel's initial angle is less than the maximum outflow
c     angle, compute the time at which the outflow angle equals the
c     parcel's initial angle.
c
      if(xz(2).gt.0.191d0*au*(xz(1)/au)**1.5d0) then
          tflow = (xz(2)/(0.191d0*au*(xz(1)/au)**1.5d0))**(-f1o3) *
     $            tcloudm
      else
          tflow = 10d0 * tmax
      endif
c
c     If t00 exceeds the maximum output time step or if the outflow
c     will reach the parcel before it begins to fall in, write the
c     initial position for each time step and skip straight to the end.
c     Terminate if and when the outflow reaches the parcel.
c
      if(t00.gt.tmax.or.t00.gt.tflow) then
          step = 2
          tstore(1) = 0d0
          rstore(1) = r0(ip)
          thstore(1) = th0(ip)
          nstore(1) = nden0
          tstore(2) = min(tmax,tflow)
          rstore(2) = r0(ip)
          thstore(2) = th0(ip)
          nstore(2) = nden0
          do j=1,de_ntime
              if(nparc.lt.50000.or.j.eq.102.or.j.eq.104) then
                  if(de_time(j).lt.tflow) then
                      temp = get_dust_temp(t00,r0(ip),th0(ip))
                      write(10,100)
     #                 j,de_time(j),r0(ip),th0(ip),max(1d0,nden0),temp
                  else
                      write(10,100) j,de_time(j),1d0,0.25d0*pi,1d0,5d0
                  endif
              endif
          enddo
          itime = de_ntime + 1
          goto 1020
      endif
c
c     Locate t00 on the diskevol time grid.
c
      itime = 1
      call hunt(de_time,de_ntime,t00,itime)
      if((itime.le.0).or.(itime.ge.de_ntime)) stop 8350
c
c     Write the initial values at all time steps up to itime.
c
      do j=1,itime
          if(nparc.lt.50000.or.j.eq.102.or.j.eq.104) then
              temp = get_dust_temp(t00,r0(ip),th0(ip))
              write(10,100) 
     $         j,de_time(j),r0(ip),th0(ip),max(1d0,nden0),temp
          endif
      enddo
  100 format(I5,2(1PE16.8),0PF12.8,1PE16.8,1PE13.5)
  101 format(I5,4(1PE16.8),0PF12.8,6(1PE16.8),1X,A5,I2)
      itime = itime + 1
c
c     If the infall starts after the first hydrostatic core phase
c     [assumed to correspond to 2.01d4*year], add two steps to the
c     storage arrays to get the correct density and temperature profile.
c
      if(t.gt.2.01d4*year) then
          tstore(1) = 0.999d0*2.01d4*year
          rstore(1) = r0(ip)
          thstore(1) = th0(ip)
          nstore(1) = nden0
          if(write_traj.eq.1) then
              write(20,101) 0,0.999d0*2.01d4*year,xz(1),xz(2),r0(ip),
     $                      th0(ip),nden0,0d0,0d0,0d0,0d0,zerotemp,
     $                      'init',0
          endif
          tstore(2) = 2.01d4*year
          rstore(2) = r0(ip)
          thstore(2) = th0(ip)
          nstore(2) = nden0
          if(write_traj.eq.1) then
              temp = get_dust_temp(2.01d4*year,r0(ip),th0(ip))
              write(20,101) 1,2.01d4*year,xz(1),xz(2),r0(ip),th0(ip),
     $                      nden0,0d0,0d0,0d0,0d0,temp,'init',0
          endif
          step = 3
      else
          step = 1
      endif
      if(write_traj.eq.1) then
          temp = get_dust_temp(t00,r0(ip),th0(ip))
          write(20,101) step-1,t00,xz(1),xz(2),r0(ip),th0(ip),nden0,0d0,
     $                  0d0,0d0,0d0,temp,'init',0
      endif
c
c     For any reasonable set of initial conditions and position, the
c     parcel will start in the TSC84 or CM81 regime. Check whether that
c     is indeed the case.
c
      rcen = 0.0625d0*cloud_asound*((m0*t)**3)*(cloud_omega**2)
      if(r0(ip)/(cloud_asound*t00).lt.(cloud_omega*t00)**2) then
          costh = cos(th0(ip))
          rrc = r0(ip) / rcen
          if(rrc.eq.1d0) then
              mu0 = costh**f1o3
          else if(costh.eq.0d0) then
              if(rrc.lt.1d0) then
                  mu0 = sqrt(1d0-rrc)
              else
                  mu0 = 0d0
              endif
          else
              acub(4) = 1d0
              acub(3) = 0d0
              acub(2) = rrc-1d0
              acub(1) = -rrc*costh
              call cbcrt(acub,xcub)
              mu0 = xcub(1)
          endif
          if(costh.gt.0d0) then
              check = mu0 * (1d0-(1d0-mu0**2)/rrc)
              if(abs(1d0-check/costh).gt.1d-6) stop 3362
          endif
          tcm81 = t00
          loc = 'cm__'
      else
          loc = 'tsc_'
      endif
!       rpar = 1e-99
c
c     Start the infall.
c
 1000 continue
c
c     Set time-dependent parameters for DVODE.
c
      atol = 1d-8 * xz(1)
      work(1) = tmax
      global_ndt0 = 0
c
c     Call DVODE to advance one time step.
c
      call dvode(fparcel,neq,xz,t,tmax,itol,rtol,atol,ntask,
     $           nstate,iopt,work,lenwd,iwork,leniwd,jacobn,mf,
     $           rpar,ipar)
!       write(*,*) nstate, neq, t/year,tmax/year
!       write(*,*)  t/year, xz(1)/au, xz(2)/au
c
c     If inside the disk, get the new altitude.
c
      zmethod = 0
      if(t.ge.tdisk) then
c
c NOTE: Sep 2 2011
c After it enters the disk it just move away from it
c NOTE: Sep 5 2011
c Need to check if I am still in the disk...
c Above no longer is true: Sep 9 2011
c Found one line that needed to change in fparcel
c It was missing the disk boundary definition
c         
          if((loc.eq.'tsc_').or.(loc.eq.'cm__')) then
              write(12,*) ip,t,r,theta
              loc = 'disk'
          endif
          lnt = log(t)
          call hunt(de_lntime,de_ntime,lnt,jtime)
          if((jtime.le.0).or.(jtime.ge.de_ntime)) then
              write(*,*) jtime,t/year,tmax/year
              stop 11001
          endif
          deeps = (lnt-de_lntime(jtime)) /
     $            (de_lntime(jtime+1)-de_lntime(jtime))
          deepsn = 1d0 - deeps
          if((deeps.lt.0.d0).or.(deeps.gt.1.d0)) stop 11002
          call get_diskvelden(xz(1),xz(2),t,jtime,deeps,deepsn,1)
          nden = mden * mpeffinv
          dr = xepsn * grid_dr(jgrid) + xeps * grid_dr(jgrid+1)
          dum = surferf / surf2
          if(dum.lt.1d0) then
              xz(2) = scale2 * sqrt(2d0) * erfi(dum)
              zmethod = 1
          else
              if(scaleold.lt.0d0) then
                  write(104,*) ip,t/year
                  step = step - 1
                  goto 1020
                  stop 10442
              endif
              xz(2) = scale2 * xz(2) / scaleold
              zmethod = 2
          endif
          scaleold = scale2
      endif

c
c     Note: DH, Sep 13, 2011
c     need to recalculate r and th in case in the disk
c
      r = sqrt(xz(1)*xz(1)+xz(2)*xz(2))
      theta = atan2(xz(1),xz(2))
c
c     Initialize some variables.
c
      xx = r / (cloud_asound * t)
      tau = cloud_omega * t
      tau2 = tau * tau
c
c     If in the envelope, get the density.
c
      outflow = 0
      if(t.ge.tcm81.and.t.lt.tdisk) then
          loc = 'cm__'
          call get_cmvelden(r,theta,t,mu0,1d0/r,0)
          nden = max(1d-50,mdenf * mpeffinv)
      endif
      if(t.lt.tcm81) then
          loc = 'tsc_'
          call get_tscvelden(r,theta,t,1d0/t,0)
          nden = max(1d-50,mdenf * mpeffinv)
      endif

c     
c     Note: DH, Sep 15, 2011
c     if xx E [0.01 and scalexx]: Need to scale the densities
c

      if( (t.lt.tdisk)) then
        if((xx.lt.upperbnd*tau2*scalexx).and.(xx.gt.lowerbnd*tau2)) then
            call get_cmvelden(r,theta,t,mu0,1d0/r,0)
            nden81 = max(1d-50,mdenf * mpeffinv)
            call get_tscvelden(r,theta,t,1d0/t,0)
            nden84 = max(1d-50,mdenf * mpeffinv)
           
            if (cloud_omega.gt.8d-14) then
                eps = (norm1) * (log10(xx/tau2) + offset)
            elseif (cloud_omega.lt.1.2d-14) then
                eps = (norm2) * (log10(xx/tau2) + offset)
            endif

            if(eps.lt.0d0.or.eps.gt.1d0) then
                write(*,121) t/year,r,xx,tau2,nden81,nden84,
     $          eps,scalexx
  121           format(8(1PE13.5))
                stop 5243
            endif
            nden = eps * nden84 + (1d0-eps) * nden81
        endif
      endif

c
c     Check if the parcel entered the outflow. If so, delete the last
c     step and stop the calculation for this parcel.
c
      if(t.lt.tdisk.and.outflow.eq.1) then
          step = step - 1
          goto 1020
      endif
c
c     If we crossed a diskevol time step, compute and write the position
c     corresponding to this step. Also check if we perhaps crossed two
c     or more diskevol steps at once.
c
 1030 if(t.gt.de_time(itime)) then
          eps = (de_time(itime)-told) / (t-told)
          if((eps.lt.0.d0).or.(eps.gt.1.d0)) then
              write(*,*) ip,eps,de_time(itime),told,t,itime
              write(*,*) de_time(de_ntime)/year,tcloudm/year
              stop 8341
          endif
          rsave = (1.d0-eps) * rold + eps * r
          thsave = (1.d0-eps) * thold + eps * theta
          ndensave = (1.d0-eps) * ndenold + eps * nden
          if(nparc.lt.50000.or.itime.eq.102.or.itime.eq.104) then
              temp = get_dust_temp(de_time(itime),rsave,thsave)
              write(10,100)
     #     itime,de_time(itime),rsave,thsave,max(1d0,ndensave),temp
          endif
          tstore(step) = de_time(itime)
          rstore(step) = rsave
          thstore(step) = thsave
          nstore(step) = max(1d0,ndensave)
          step = step + 1
          itime = itime + 1
          if(itime.le.de_ntime) goto 1030
      endif
c
c     Store the results from this step.
c
      tstore(step) = t
      rstore(step) = r
      thstore(step) = theta
      nstore(step) = nden
      if(write_traj.eq.1) then
          temp = get_dust_temp(t,r,theta)
          write(20,101) step,t,xz(1),xz(2),r,theta,nden,scale2,surf2,
     $                  ux,uz,temp,loc,zmethod
      endif
c
c     Did we reach the star or cross the maximum time? If not, save the
c     current values and do another loop. If we did, stop the iteration
c     for this parcel.
c
      if(xz(1).lt.rin.or.xz(2).lt.0d0.or.t.gt.tmax) then
          goto 1020
      endif
      told = t
      rold = r
      thold = theta
      ndenold = nden
c
c     Did we reach the max number of steps? If so, report the problem
c     and stop the iteration. If not, go to the next step.
c
      if(step.eq.NSTORE) then
          write(*,*) 'Max. number of steps reached. You may want to ',
     $               'recompile with a larger NSTORE.'
          write(*,*) ip,loc,ux,uz, told/year, de_time(de_ntime)/year
          write(*,*) step, NSTORE, rold/au, thold
          write(*,*) rstore(step-1)/au, thstore(step-1)
!           stop 900
          goto 1020  
      endif
      step = step + 1
      goto 1000
 1020 continue
c
c     Write values for the remaining steps (radial distance of 1 cm
c     and angle of pi/4).
c
      if(itime.le.de_ntime) then
          do j=itime,de_ntime
              if(nparc.lt.50000.or.j.eq.102.or.j.eq.104) then
                  write(10,100) j,de_time(j),5d0,0.78539816d0,1d0,5d0
              endif
          enddo
      endif
      if(write_traj.eq.1) then
          write(20,101) -1,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,
     $                  0d0,0d0,0d0,'end_',0
      endif
c
c     Create an input file for the chemistry with the temperature,
c     density and time derivative of the density at a number of time
c     steps.
c
      if(dowritechem.eq.1) then
          call write_chem_input(ip,do_tau,do_cdh)
      endif
c
      return
      end
c
c
c
c     ------------------------------------------------------------------
c     This subroutine computes the velocities in two-dimensional
c     cylindrical coordinates as a function of time and position.
c     ------------------------------------------------------------------
      subroutine fparcel(n,t,y,ydot,rpar,ipar)
c
      implicit none
      doubleprecision t,y(n),ydot(n),rpar
      integer n,ipar
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
#include "common_traj.h"
c
      integer jtime,i,nroot
      doubleprecision lnt,tinv,r,theta,rinv,deeps,deepsn,xx,rrc
      doubleprecision nden,ndenenv,ndendisk,check
      doubleprecision prsdisk,prsenv,tcloud,get_lstar
c Added by DH Sep 4, 2011 -> scale the CM81 boundary
      doubleprecision scalexx
      doubleprecision costh,dr,dummu0,xcub(3),acub(4)
c     Added by Daniel, Jun 11, 2011
c     Need variables to define the transition boundary between
c     TSC and CM solutions
      doubleprecision upperbnd, lowerbnd,norm1,norm2
      doubleprecision offset
c     Attempt to fix DVODE error: DH 2013
      doubleprecision told, dt
      character*4 oldloc
      upperbnd = 1d1
      lowerbnd = 1d-2
      offset = 2.0d0
      norm1 = 0.3333333d0
      norm2 = 0.5d0
c
c     If the time .le. 0, return zero velocities.
c
      if(t.le.0d0) then
          ydot(1) = 0d0
          ydot(2) = 0d0
          return
      endif
c
c     Store the log and inverse of the time for later use. Convert
c     cylindrical to polar coordinates.
c
      lnt = log(t)
      tinv = 1d0 / t
      r = sqrt(y(1)**2+y(2)**2)
      theta = atan2(y(1),y(2))
      rinv = 1d0 / r

c     Added by Daniel: Sep 4, 2011
c     need to scale the xx for the transition between solution
c     this depends on the rc -> smaller rotation lead to smaller rc
c     hence this means that tsc84 solution is valid into smaller scale
      
      if (cloud_omega.gt.8d-14) then
        scalexx = 1d0
      elseif (cloud_omega.lt.1.2d-14) then
        scalexx = 1d-1
      else
        scalexx = (cloud_omega/1d-13)*(cloud_omega/1d-13)
      endif
c
c     Added 10-10,2011
c     parameters for get_vertden
c
      rcloud = GG * cloud_mass / (2d0*cloud_asound**2)
      tcloud = rcloud / cloud_asound
      if(t.lt.tcloudm) then
          ml_mcentr = (0.5d0*m0) * cloud_mass * (t/tcloud)
          mdot = m0 * cloud_asound**3 * GGinv
          ml_rcentr = 0.0625d0 * cloud_asound * ((m0*t)**3) *
     $                (cloud_omega**2)
      else
          ml_mcentr = cloud_mass
          mdot = 0d0
          ml_rcentr = 0.0625d0 * cloud_asound * ((m0*tcloudm)**3) *
     $                (cloud_omega**2)
      endif
      time_time = t
      global_mstar = mstar
      call get_rstar()
      global_lstar = get_lstar()
c
c     Find the current time step on the disk evolution grid. Set the
c     mass of the star and the disk accordingly.
c
      jtime = 1
      call hunt(de_lntime,de_ntime,lnt,jtime)
      if((jtime.le.0).or.(jtime.ge.de_ntime)) then
          write(*,*) t,lnt,t/year,jtime,de_ntime
          write(*,*) exp(de_lntime(1))/year, 
     $     exp(de_lntime(de_ntime))/year
          stop 10001
      endif
      deeps = (lnt-de_lntime(jtime)) /
     $        (de_lntime(jtime+1)-de_lntime(jtime))
      deepsn = 1d0 - deeps
      if((deeps.lt.0.d0).or.(deeps.gt.1.d0)) stop 10002
      mstar = deepsn * de_mstar(jtime) + deeps * de_mstar(jtime+1)
      mdisk = deepsn * de_mdisk(jtime) + deeps * de_mdisk(jtime+1)
c
c     If we were not yet in the CM81 regime, check if we just got
c     there.
c
      if(t.lt.tcm81) then
          xx = r * cloud_asound_inv * tinv
          if(xx.lt.lowerbnd*(cloud_omega*t)**2) then
              tcm81 = t
              rdisk = 0.0625d0 * cloud_asound * (m0*t)**3 *
     $                cloud_omega**2
              rrc = r / rdisk
              costh = y(2) * rinv
              if(rrc.eq.1d0) then
                  mu0 = costh**f1o3
              else if(costh.eq.0d0) then
                  if(rrc.lt.1d0) then
                      mu0 = sqrt(1d0-rrc)
                  else
                      mu0 = 0d0
                  endif
              else
                  acub(4) = 1d0
                  acub(3) = 0d0
                  acub(2) = rrc-1d0
                  acub(1) = -rrc*costh
                  call cbcrt(acub,xcub)
                  mu0 = xcub(1)
              endif
              if(dummu0.lt.-2d0) dummu0 = mu0
              if(costh.gt.0d0) then
                  check = mu0 * (1d0-(1d0-mu0**2)/rrc)
                  if(abs(1d0-check/costh).gt.1d-6) then
                      write(*,*) r,rdisk,costh,mu0,check
                      stop 10003
                  endif
              endif
          endif
      endif
c
c     The method of computing the velocity depends on the type of regime
c     the parcel is currently in. If we're in the TSC84 or CM81 regime,
c     we also need to get the disk density to determine if the
c     parcel has entered the disk. The 'velden' routines overwrite the
c     value of the global parameter mden, so we have to store the
c     densities in other variables.
c
      if(t.gt.tdisk) then
          call get_diskvelden(y(1),y(2),t,jtime,deeps,deepsn,1)
          uz = 0d0
c      elseif(t.gt.tcm81) then
c          call get_cmvelden(r,theta,t,mu0,rinv,0)
c          call get_vertden(y(1),y(2),r,theta,t,-1d0,0)
c          ndenenv = mden * mpeffinv
c          prsenv =  mden * (ux*ux + uz*uz) 
c          call get_diskvelden(y(1),y(2),t,jtime,deeps,deepsn,0)
c          ndendisk = mden * mpeffinv
c          prsdisk = ndendisk * kk * dktemp
      elseif(t.lt.tcm81) then
          call get_tscvelden(r,theta,t,tinv,0)
          ndenenv = mden * mpeffinv
          prsenv = mden * (ux*ux + uz*uz)
          call get_diskvelden(y(1),y(2),t,jtime,deeps,deepsn,0)
          ndendisk = mden * mpeffinv
          prsdisk = ndendisk * kk * dktemp
      else
c     Added by DH, 10-10-2011:
c     Need to do the same thing as in diskevol: interpolation between
c     TSC and CM solutins
          call get_vertden(y(1),abs(y(2)),r,theta,t,-1d0,0)
          ndenenv = max(1d0,mdenf*mpeffinv)
          prsenv = mden * (ux*ux + uz*uz)
          call get_diskvelden(y(1),y(2),t,jtime,deeps,deepsn,1)
          ndendisk = mden * mpeffinv
          prsdisk = ndendisk * kk *dktemp
      endif

c
c     If we were not yet in the disk regime, check if we just got
c     there.
c
      if(t.lt.tdisk) then
c
c Need to change this to pressure!!!
c Changed by Daniel Sep 2, 2011
c Modified by DH: Dec, 2013:
c  Make sure that it calls the functions again to get the right
c  velocities
c
c          if(ndendisk.gt.ndenenv) then
          if(prsdisk.gt.prsenv) then
              tdisk = t
              if(t.lt.tcm81) tcm81 = t
              dr = xepsn * grid_dr(jgrid) + xeps * grid_dr(jgrid+1)
              surferf = surf2 * erf(y(2)/(scale2*sqrt(2d0)))
              call get_diskvelden(y(1),y(2),t,jtime,deeps,deepsn,1) 
          else 
              call get_vertden(y(1),abs(y(2)),r,theta,t,-1d0,0)
          endif
      endif
c
c     If we are in the disk.... check if we get out of the disk
c     added by DH: Dec 2013
c
      if(t.gt.tdisk) then
          call get_vertden(y(1),abs(y(2)),r,theta,t,-1d0,0)
          ndenenv = max(1d0,mdenf*mpeffinv)
          prsenv = mden * (ux*ux + uz*uz)
          call get_diskvelden(y(1),y(2),t,jtime,deeps,deepsn,1)
          ndendisk = mden * mpeffinv
          prsdisk = ndendisk * kk *dktemp
          if(prsdisk.lt.prsenv) then
c                 We have left the disk
                call get_vertden(y(1),abs(y(2)),r,theta,t,-1d0,0)
          else
                call get_diskvelden(y(1),y(2),t,jtime,deeps,deepsn,1)
                uz = 0d0
          endif
      endif
c    DONE MODIFICATIONS DH: DEC 2013
c  
c     Store the velocities in the ydot array.
c
      ydot(1) = ux
      ydot(2) = uz
c
c    This is a hack...
      if (y(2)/au.lt.0.9) then
        ydot(2) = 0d0
      endif
      if (global_dt/year.le.1e-3) then
        global_ndt0 = global_ndt0+1
      endif
!       write(*,*) t/year, global_told/year, global_dt/year
!       write(*,*) y(1)/au, y(2)/au, outflow, ydot(1)/au, ydot(2)/au
!       write(*,*) global_loc, global_oldloc, global_ndt0
!       write(*,*)
c
c    DH 2013: attempt to fix DVODE? 
c
      if (global_ndt0.gt.30) then
          ydot(1) = 0d0
          ydot(2) = 0d0
      endif
      global_dt = abs(global_told-t)
      global_told = t
      global_oldloc = global_loc
c
      return
      end
c
c
c
c     ------------------------------------------------------------------
c     This subroutine computes the infall trajectory for a single
c     parcel using the DVODE solver.
c     ------------------------------------------------------------------
      subroutine make_trajvodeback(ip,write_chem,do_tau,do_cdh)
c
      implicit none
      integer ip,write_chem,do_tau,do_cdh,dowritechem
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
#include "common_traj.h"
c
      integer j,itol,ntask,iopt,mf,neq,leniwd,lenwd,nstate,ipar,jtime
      integer zmethod,itime,i,nroot
      parameter(leniwd=32)
      parameter(lenwd=62) !48
      integer iwork(leniwd)
      doubleprecision dummy(NSTORE)
      doubleprecision xz(2),rtol,atol,tnew,t,t00,rpar,lnt,deeps,deepsn
      doubleprecision work(lenwd),dr,dum,erfi,scaleold,rrc,costh,r,rinv
      doubleprecision check,dummu0,mdendisk,mdenenv
      doubleprecision rcen,rhom,theta,ztemp,eps,nden,nden0,temp,tcloud
      doubleprecision ndendisk,ndenenv,outangle,rold,rsave,tflow,thold
      doubleprecision told,get_dust_temp,tinv,aout,thsave,lnx,xcub(3)
      doubleprecision get_lstar,acub(4),prsenv,prsdisk,ndenold,ndensave
      character*4 loc
      external fparcelback,jacobn
c
c     Initialize several parameters. The time is initialized to the
c     point where the expansion wave reaches the initial position
c     (r = c_s t).
c
      xz(1) = r0(ip)
      xz(2) = th0(ip)
      t = t0(ip)
      t00 = 1d0
      if(t.gt.1.01d0*tmax) then
          write(*,*) 't > tmax for ip =',ip
          write(*,*) t/year,tmax/year
          stop
      endif
      if(t.gt.0.999999d0*tmax) t = 0.999999d0 * tmax
      r = sqrt(xz(1)**2+xz(2)**2)
      rinv = 1d0 / r
      theta = atan2(xz(1),xz(2))
      tcm81 = 0d0
      ttsc84 = 0d0
      told = t
      rold = r
      thold = theta
      zmethod = 0
      dowritechem = write_chem
c
c     Get the surface density for the initial position.
c
      lnt = log(t)
      call hunt(de_lntime,de_ntime,lnt,jtime)
      if((jtime.le.0).or.(jtime.ge.de_ntime)) then
          write(*,*) jtime,t/year,tmax/year
          stop 11001
      endif
      deeps = (lnt-de_lntime(jtime)) /
     $        (de_lntime(jtime+1)-de_lntime(jtime))
      deepsn = 1d0 - deeps
      if((deeps.lt.0.d0).or.(deeps.gt.1.d0)) stop 11002
      mstar = deepsn * de_mstar(jtime) + deeps * de_mstar(jtime+1)
      mdisk = deepsn * de_mdisk(jtime) + deeps * de_mdisk(jtime+1)
      call get_diskvelden(xz(1),xz(2),t,jtime,deeps,deepsn,0)
      surferf = surf2 * erf(xz(2)/(scale2*sqrt(2d0)))
c
c     Locate the current position on the disk evolution grid.
c
      lnx = log(xz(1))
      jgrid = 1
      call hunt(grid_lnr,de_ngrid,lnx,jgrid)
      if(jgrid.le.0.or.jgrid.ge.de_ngrid) then
          write(*,500) ip,jgrid,de_ngrid,xz(1),grid_r(1),
     #     grid_r(de_ngrid)
  500     format(I5,I5,I5,3(1PE13.5),I5)
          stop 8400
      endif
      xeps = (lnx-grid_lnr(jgrid)) / (grid_lnr(jgrid+1)-grid_lnr(jgrid))
      xepsn = 1d0 - xeps
      if((xeps.lt.0.d0).or.(xeps.gt.1.d0)) stop 8401
c
c     Compute the location of the envelope/disk boundary.
c
      surface = deepsn * (xepsn * de_surface(jtime,jgrid) +
     $                    xeps  * de_surface(jtime,jgrid+1)) +
     $          deeps  * (xepsn * de_surface(jtime+1,jgrid) +
     $                    xeps  * de_surface(jtime+1,jgrid+1))
c      write(*,47) jtime,jgrid,xeps,xepsn,deeps,deepsn
c   47 format(I5,I5,4(0PF9.4))
c      write(*,*) 'Surface:',surface/au
c
c     Get the stellar properties.
c
      rcloud = GG * cloud_mass / (2d0*cloud_asound**2)
      tcloud = rcloud / cloud_asound
      if(t.lt.tcloudm) then
          ml_mcentr = (0.5d0*m0) * cloud_mass * (t/tcloud)
          mdot = m0 * cloud_asound**3 * GGinv
          ml_rcentr = 0.0625d0 * cloud_asound * ((m0*t)**3) *
     $                (cloud_omega**2)
      else
          ml_mcentr = cloud_mass
          mdot = 0d0
          ml_rcentr = 0.0625d0 * cloud_asound * ((m0*tcloudm)**3) *
     $                (cloud_omega**2)
      endif
      time_time = t
      global_mstar = mstar
      call get_rstar()
      global_lstar = get_lstar()
c      write(*,*) global_mstar/ms,global_rstar/rs,global_lstar/ls
c
c     Initialize constant parameters for DVODE.
c
      neq = 2
      itol = 1
      rtol = 1d-8
      ntask = 2
      work(5) = 0d0
      work(6) = 1d12
      work(7) = 0d0
      iwork(5) = 0
      iwork(6) = 30000
      iwork(7) = 0
      iopt = 1
      mf = 12 !22
      nstate = 1
c
c     Are we starting in the disk or in the envelope?
c
      tinv = 1d0 / t
      lnt = log(t)
      call hunt(de_lntime,de_ntime,lnt,jtime)
      deeps = (lnt-de_lntime(jtime)) /
     $        (de_lntime(jtime+1)-de_lntime(jtime))
      deepsn = 1d0 - deeps
      if((deeps.lt.0.d0).or.(deeps.gt.1.d0)) stop 11002
      if((jtime.le.0).or.(jtime.ge.de_ntime)) then
          write(*,*) jtime,t/year,tmax/year
          stop 11001
      endif
c      write(*,48) xz(1),xz(2),t
c      write(*,48) ndendisk,ndenenv,mden*mpeffinv
      call get_vertden(xz(1),xz(2),r,theta,t,-1d0,0)
      ndenenv = mden * mpeffinv
      prsenv = prs
      mu0 = global_mu0
      call get_diskvelden(xz(1),xz(2),t,jtime,deeps,deepsn,0)
      ndendisk = mden * mpeffinv
      prsdisk = ndendisk * kk*dktemp
c
c     Sept 13: A bug here? disk pressure > env issue?
c     Oct 13: redefine the tcm81 and ttsc84
c
      if(global_loc.eq.'cm__') then
c          if(ndenenv.gt.ndendisk) then
          if(prsdisk.gt.prsenv) then
!               tcm81 = t
              loc = 'disk'
              nden0 = ndendisk
          else
              loc = 'cm__'
              tcm81 = t
              nden0 = ndenenv
          endif
      else
c          if(ndenenv.gt.ndendisk) then
          if(prsdisk.gt.prsenv) then
!               tcm81 = t
!               ttsc84 = t
              loc = 'disk'
              nden0 = ndendisk
          else
              tcm81 = t
              ttsc84 = t
              loc = 'disk'
              nden0 = ndenenv
          endif
      endif
      ndenold = nden0
c      write(*,48) ndendisk,ndenenv,mden
c   48 format(3(1PE16.9))
c      stop 48
c
c     If we're starting in the envelope, check if we're starting inside
c     the outflow.
c
      outflow = 0
      if(loc.eq.'tsc_'.or.loc.eq.'cm__') then
          outangle = 0.5519d0*(t/tcloudm)**0.26d0
          if(outangle.gt.pi/2) outangle = pi/2
          aout = 0.191d0 / (t/tcloudm)**3
          if(aout.lt.(0.191d0)) aout = 0.191d0
          if(xz(2).gt.aout*au*(xz(1)/au)**1.5d0) then
              outflow = 1
              tflow = (xz(2)/(0.191d0*au*(xz(1)/au)**1.5d0))**(-f1o3) *
     $                tcloudm
              t = 0.999999d0 * tflow
              nden0 = mpeffinv * 3d0 * mflow /
     $            (2d0*pi*rcloud**3*(1d0-cos(0.5d0*theta)))
              loc = 'flow'
          endif
      endif
!       write(*,*) 1, tmax/year, xz(1)/au, xz(2)/au, loc
c
c     Store the position at tmax and at the first time step.
c
      tstore(1) = tmax
      rstore(1) = r
      thstore(1) = theta
      nstore(1) = nden0
      if(write_traj.eq.1) then
          temp = get_dust_temp(tmax,r,theta)
          write(20,101) 1,tmax,xz(1),xz(2),r,theta,nden0,0d0,
     $                  0d0,0d0,0d0,temp,'init',0
      endif
      tstore(2) = t
      rstore(2) = r
      thstore(2) = theta
      nstore(2) = nden0
      if(write_traj.eq.1) then
          temp = get_dust_temp(t,r,theta)
          write(20,101) 2,t,xz(1),xz(2),r,theta,nden0,0d0,
     $                  0d0,0d0,0d0,temp,'init',0
      endif
      step = 3
c
c     If we start in the outflow, see if the parcel is not in the
c     outflow at an earlier time.
c
      if(loc.eq.'flow') then
          if(tflow.lt.de_time(2)) then
              itime = de_ntime
              tstore(3) = 0d0
              rstore(3) = r
              thstore(3) = theta
              nstore(3) = nden0
              if(write_traj.eq.1) then
                  write(20,101) 3,0d0,xz(1),xz(2),r,theta,nden0,0d0,
     $                          0d0,0d0,0d0,10d0,'flow',0
              endif
              goto 210
          endif
      endif
c
c     Locate t on the diskevol time grid.
c
      itime = 1
      call hunt(de_time,de_ntime,t,itime)
      if((itime.le.0).or.(itime.ge.de_ntime)) stop 8350
c
c     Write the initial values at all time steps up to itime.
c
      do j=de_ntime,itime+1,-1
c          if(nparc.lt.50000.or.j.eq.102.or.j.eq.104) then
c           if(
           temp = get_dust_temp(t00,r0(ip),th0(ip))
           write(10,100) j,de_time(j),r,theta,max(1d0,nden0),temp
c          endif
      enddo
  100 format(I5,2(1PE16.8),0PF12.8,1PE16.8,1PE13.5)
  101 format(I5,4(1PE16.8),0PF12.8,6(1PE16.8),1X,A5,I2)
c
c     Start the infall.
c
      scaleold = -1d0
 1000 continue
c
c     Set time-dependent parameters for DVODE.
c
      atol = 1d-8 * xz(1)
      work(1) = t00
c
c     Call DVODE.
c
c      if(t.gt.tcm81) nstate = 1
      call dvode(fparcelback,neq,xz,t,t00,itol,rtol,atol,ntask,
     $           nstate,iopt,work,lenwd,iwork,leniwd,jacobn,mf,
     $           rpar,ipar)
      r = sqrt(xz(1)**2+xz(2)**2)
      theta = atan2(xz(1),xz(2))
!       if (r/au.lt.200) write(*,*) step, t/year, xz(1)/au, xz(2)/au, loc
c
c     If inside the disk, get the new altitude.
c
      if(t.gt.tcm81) then
!           write(*,*) 'Inside the disk'
          lnt = log(t)
          call hunt(de_lntime,de_ntime,lnt,jtime)
          if((jtime.le.0).or.(jtime.ge.de_ntime)) then
              write(*,*) jtime,t/year,tmax/year
              stop 11001
          endif
          deeps = (lnt-de_lntime(jtime)) /
     $            (de_lntime(jtime+1)-de_lntime(jtime))
          deepsn = 1d0 - deeps
          if((deeps.lt.0.d0).or.(deeps.gt.1.d0)) stop 11002
          call get_diskvelden(xz(1),xz(2),t,jtime,deeps,deepsn,0)
          dr = xepsn * grid_dr(jgrid) + xeps * grid_dr(jgrid+1)
          dum = surferf / surf2
          if(dum.lt.1d0) then
              xz(2) = scale2 * sqrt(2d0) * erfi(dum)
              zmethod = 1
          else
              xz(2) = scale2 * xz(2) / scaleold
              zmethod = 2              
! c
! c     Scan upwards and find the point where the envelope density equals
! c     the disk density.
! c
!               rcen = 0.0625d0 * cloud_asound * (m0*t)**3 *
!      $               cloud_omega**2
!               rhom = surf2 / (scale2 * sqrt(2d0*pi))
!               ztemp = xz(2)
!   200         continue
!               r = sqrt(xz(1)**2+ztemp**2)
!               rinv = 1d0 / r
!               rrc = r / rcen
!               theta = atan2(xz(1),ztemp)
!               costh = ztemp * rinv
!               if(rrc.eq.1d0) then
!                   mu0 = costh**f1o3
!               else if(costh.eq.0d0) then
!                   if(rrc.lt.1d0) then
!                       mu0 = sqrt(1d0-rrc)
!                   else
!                       mu0 = 0d0
!                   endif
!               else
!                   acub(4) = 1d0
!                   acub(3) = 0d0
!                   acub(2) = rrc-1d0
!                   acub(1) = -rrc*costh
!                   call cbcrt(acub,xcub)
!                   mu0 = xcub(1)
!               endif
!               if(dummu0.lt.-2d0) dummu0 = mu0
!               if(costh.gt.0d0) then
!                   check = mu0 * (1d0-(1d0-mu0**2)/rrc)
!                   if(abs(1d0-check/costh).gt.1d-6) then
!                       write(*,*) r,ml_rcentr,costh,mu0,check
!                       stop 10005
!                   endif
!               endif
! !               call get_cmvelden(r,theta,t,mu0,rinv,0)
! !             DH, Dec 2013: vertden has the new interpolation
!               call get_vertden(xz(1),xz(2),r,theta,t,-1d0,0)
!               mdenenv = mden
!               prsenv = mden * (ux**2+uz**2)
!               mdendisk = rhom * exp(-0.5d0*(ztemp/scale2)**2)
!               prsdisk = mdendisk * mpeffinv * kk*dktemp
!               if(prsenv.gt.prsdisk) then
!                   ztemp = ztemp * 1.001d0
!                   if(ztemp.gt.1d3*au) then
!                       write(*,*) ip,r0(ip)/au,th0(ip)/au
!                       write(*,*) xz(1)/au,xz(2)/au,ztemp/au
!                       stop 9100
!                   endif
!                   goto 200
!               endif
!               xz(2) = ztemp
!               tcm81 = t
          endif
          scaleold = scale2
      endif
c
c     Get the density.
c     DH:modify Dec 2013: 
c     call get_vertden and disk to get disk surface
c
      outflow = 0
      call get_diskvelden(xz(1),xz(2), t, jtime, deeps, deepsn,0)
      ndendisk = mden*mpeffinv
      prsdisk = ndendisk*kk*dktemp
      call get_vertden(xz(1),xz(2),r,theta,t,-1d0,0)
      nden = mden*mpeffinv
      prsenv = mden*(ux**2+uz**2)
      if (prsdisk.gt.prsenv) then
        call get_diskvelden(xz(1),xz(2), t, jtime, deeps, deepsn,0)
        nden = mden*mpeffinv
        loc  = 'disk'
      else
        call get_vertden(xz(1),xz(2),r,theta,t,-1d0,0)
        loc = global_loc
        if (loc.eq.'cm__') then
            tcm81 = t
        else
            ttsc84 = t
        endif
        nden = (mden*mpeffinv) + ndendisk
      endif
! 
!       if(t.gt.tcm81) then
!           call get_diskvelden(xz(1),xz(2),t,jtime,deeps,deepsn,0)
!           nden = mden * mpeffinv
!           loc = 'disk'
!       elseif(t.le.tcm81.and.t.gt.ttsc84) then
!           call get_vertden(xz(1),xz(2),r,theta,t,-1d0,0)
!           nden = mden * mpeffinv
!           if(loc.eq.'disk') write(12,*) ip,t,r,theta
!           loc = 'cm__'
!       else
!           call get_vertden(xz(1),xz(2),r,theta,t,-1d0,0)
!           nden = mden * mpeffinv
!           if(loc.eq.'disk') write(12,*) ip,t,r,theta
!           loc = 'tsc_'
!       endif
!       nden = nden + ndendisk
c
c     Check if the parcel entered the outflow.
c
      if(t.le.tcm81.and.outflow.eq.1) then
c      if (outflow.eq.1) then
          tstore(step) = t
          rstore(step) = r
          thstore(step) = theta
          nstore(step) = mpeffinv * 3d0 * mflow /
     $                   (2d0*pi*rcloud**3*(1d0-cos(0.5d0*theta)))
          if(write_traj.eq.1) then
              write(20,101) step,t,xz(1),xz(2),r,theta,nstore(step),
     $                      0d0,0d0,0d0,0d0,10d0,'flow',0
          endif
          r = 1d0
          theta = 0.25d0 * pi
          nden = 1d0
          dowritechem = 0
          goto 210
      endif
c
c     If we crossed a diskevol time step, compute and write the position
c     corresponding to this step. Also check if we perhaps crossed two
c     or more diskevol steps at once.
c
 1030 if(t.lt.de_time(itime)) then
          eps = (told-de_time(itime)) / (told-t)
          if((eps.lt.0.d0).or.(eps.gt.1.d0)) then
              write(*,*) ip,eps,de_time(itime),told,t,itime
              write(*,*) de_time(de_ntime)/year,tcloudm/year
              stop 8341
          endif
          rsave = (1.d0-eps) * rold + eps * r
          thsave = (1.d0-eps) * thold + eps * theta
          ndensave = (1.d0-eps) * ndenold + eps * nden
c          ndensave = ndensave + ndendisk
c          if(nparc.lt.50000.or.itime.eq.102.or.itime.eq.104) then
          temp = get_dust_temp(de_time(itime),rsave,thsave)
          write(10,100) 
     $       itime,de_time(itime),rsave,thsave,max(1d0,nden),temp
c          endif
          tstore(step) = de_time(itime)
          rstore(step) = rsave
          thstore(step) = thsave
          nstore(step) = ndensave
          step = step + 1
          itime = itime - 1
          if(itime.gt.1) goto 1030
      endif
c
c     Store the results from this step.
c
      tstore(step) = t
      rstore(step) = r
      thstore(step) = theta
      nstore(step) = nden
      if(write_traj.eq.1) then
          temp = get_dust_temp(t,r,theta)
          write(20,101) step,t,xz(1),xz(2),r,theta,nden,scale2,surf2,
     $                  ux,uz,temp,loc,zmethod
      endif
c
c     Did we reach the point where the parcel is at the edge of the
c     expansion wave? If so, stop the iteration. If not, save the
c     current values and do another loop.
c
      if(t.lt.r*cloud_asound_inv) goto 210
      told = t
      rold = r
      thold = theta
      ndenold = nden
c
c     Did we reach the max number of steps? If so, report the problem
c     and stop the iteration. If not, go to the next step.
c
      if(step.eq.NSTORE) then
          write(*,*) 'Max. number of steps reached. You may want to ',
     $               'recompile with a larger NSTORE.'
          write(*,*) ip,t
c          stop 900
          goto 210
      endif
      step = step + 1
c
c     Write results.
c          
      goto 1000
  210 continue
c
c     Write values for the remaining time steps down to zero.
c
      if(itime.gt.0) then
          do j=itime,1,-1
c              if(nparc.lt.50000.or.j.eq.102.or.j.eq.104) then
               temp = get_dust_temp(de_time(j),r,theta)
c   Debugging segmentation faults
               write(10,100) 
     #     j,de_time(j),r,theta,max(1d-50,nden),temp
c              endif
          enddo
      endif
      if(write_traj.eq.1) then
          write(20,101) -1,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,
     $                  0d0,0d0,0d0,'end_',0
      endif
c
c     Reverse the storage arrays.
c
      do j=1,step
          dummy(j) = tstore(step+1-j)
      enddo
      do j=1,step
          tstore(j) = dummy(j)
      enddo
      do j=1,step
          dummy(j) = rstore(step+1-j)
      enddo
      do j=1,step
          rstore(j) = dummy(j)
      enddo
      do j=1,step
          dummy(j) = thstore(step+1-j)
      enddo
      do j=1,step
          thstore(j) = dummy(j)
      enddo
      do j=1,step
          dummy(j) = nstore(step+1-j)
      enddo
      do j=1,step
          nstore(j) = dummy(j)
      enddo
c
c     Create an input file for the chemistry with the temperature,
c     density and time derivative of the density at a number of time
c     steps.
c
      if(dowritechem.eq.1) then
c          stop
          call write_chem_input(ip,do_tau,do_cdh)
      endif
c
      return
      end
c
c
c
c     ------------------------------------------------------------------
c     This subroutine computes the velocities in two-dimensional
c     cylindrical coordinates as a function of time and position.
c     ------------------------------------------------------------------
      subroutine fparcelback(nxz,t,y,ydot,rpar,ipar)
c
      implicit none
      doubleprecision t,y(nxz),ydot(nxz),rpar
      integer nxz,ipar,i,nroot
c
#include "common_natconst.h"
#include "common_disk.h"
#include "common_infall.h"
#include "common_traj.h"
c
      integer jtime
      doubleprecision lnt,tinv,r,theta,rinv,deeps,deepsn,xx,rrc
      doubleprecision nden,ndenenv,ndendisk,check,prsenv,prsdisk
      doubleprecision costh,dr,xcub(3),acub(4),tcloud
C
C Type declarations for labeled COMMON block DVOD01 --------------------
C
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
c
c     If the time .le. 0, return zero velocities.
c
      if(t.le.0d0) then
          ydot(1) = 0d0
          ydot(2) = 0d0
          return
      endif
c
c     Store the log and inverse of the time for later use. Convert
c     cylindrical to polar coordinates.
c
      lnt = log(t)
      tinv = 1d0 / t
      r = sqrt(y(1)**2+y(2)**2)
      theta = atan2(y(1),y(2))
      rinv = 1d0 / r
      rcloud = GG * cloud_mass / (2d0*cloud_asound**2)
      tcloud = rcloud / cloud_asound
      if(t.lt.tcloudm) then
          ml_mcentr = (0.5d0*m0) * cloud_mass * (t/tcloud)
          mdot = m0 * cloud_asound**3 * GGinv
          ml_rcentr = 0.0625d0 * cloud_asound * ((m0*t)**3) *
     $                (cloud_omega**2)
      else
          ml_mcentr = cloud_mass
          mdot = 0d0
          ml_rcentr = 0.0625d0 * cloud_asound * ((m0*tcloudm)**3) *
     $                (cloud_omega**2)
      endif
c
c     Find the current time step on the disk evolution grid. Set the
c     mass of the star and the disk accordingly.
c
      jtime = 1
      call hunt(de_lntime,de_ntime,lnt,jtime)
      if((jtime.le.1).or.(jtime.ge.de_ntime)) stop 10011
      deeps = (lnt-de_lntime(jtime)) /
     $        (de_lntime(jtime+1)-de_lntime(jtime))
      deepsn = 1d0 - deeps
      if((deeps.lt.0.d0).or.(deeps.gt.1.d0)) stop 10012
      mstar = deepsn * de_mstar(jtime) + deeps * de_mstar(jtime+1)
      mdisk = deepsn * de_mdisk(jtime) + deeps * de_mdisk(jtime+1)
c
c     If we were in the disk, check if we just left it.
c
      if(t.gt.tcm81) then
          xx = r * cloud_asound_inv * tinv
          if(xx.lt.(cloud_omega*t)**2) then
              rdisk = 0.0625d0 * cloud_asound * (m0*t)**3 *
     $                cloud_omega**2
              rrc = r / rdisk
              costh = y(2) * rinv
              if(rrc.eq.1d0) then
                  mu0 = costh**f1o3
              else if(costh.eq.0d0) then
                  if(rrc.lt.1d0) then
                      mu0 = sqrt(1d0-rrc)
                  else
                      mu0 = 0d0
                  endif
              else
                  acub(4) = 1d0
                  acub(3) = 0d0
                  acub(2) = rrc-1d0
                  acub(1) = -rrc*costh
                  call cbcrt(acub,xcub)
                  mu0 = xcub(1)
              endif
              if(costh.gt.0d0) then
                  check = mu0 * (1d0-(1d0-mu0**2)/rrc)
                  if(abs(1d0-check/costh).gt.1d-6) then
                      write(*,*) r,rdisk,costh,mu0,check
                      stop 10006
                  endif
              endif
c              call get_cmvelden(r,theta,t,mu0,rinv,0)
              call get_vertden(y(1), y(2),r,theta,t,-1d0,0)
              ndenenv = mden * mpeffinv
              prsenv = mden * (ux**2 + uz**2)
              call get_diskvelden(y(1),y(2),t,jtime,deeps,deepsn,0)
              ndendisk = mden * mpeffinv
              prsdisk = ndendisk * kk*dktemp
c              if(ndenenv.gt.ndendisk) then
              if(prsenv.gt.prsdisk) then
                  tcm81 = t
              endif
          else
              call get_tscvelden(r,theta,t,tinv,0)
              ndenenv = mden * mpeffinv
              prsenv = mden * (ux**2+uz**2)
              call get_diskvelden(y(1),y(2),t,jtime,deeps,deepsn,0)
              ndendisk = mden * mpeffinv
              prsdisk = ndendisk * kk*dktemp
c              if(ndenenv.gt.ndendisk) then
              if(prsenv.gt.prsdisk) then
                  tcm81 = t
                  ttsc84 = t
              endif
          endif
      endif
c
c     If we were in the CM81 regime, check if we just left it.
c
      if(t.lt.tcm81.and.t.gt.ttsc84) then
          xx = r * cloud_asound_inv * tinv
          if(xx.gt.(cloud_omega*t)**2) then
              ttsc84 = t
          endif
      endif
c
c     The method of computing the velocity depends on the type of regime
c     the parcel is currently in. If we're in the TSC84 or CM81 regime,
c     we also need to get the disk density to determine if the
c     parcel has entered the disk. The 'velden' routines overwrite the
c     value of the global parameter mden, so we have to store the
c     densities in other variables.
c
      if(t.le.ttsc84) then
          call get_tscvelden(r,theta,t,tinv,0)
          ndenenv = mden * mpeffinv
      elseif(t.le.tcm81) then
c          call get_cmvelden(r,theta,t,mu0,rinv,0)
          call get_vertden(y(1),y(2),r,theta,t,-1d0,0)
          ndenenv = mden * mpeffinv
      else
          call get_diskvelden(y(1),y(2),t,jtime,deeps,deepsn,1)
          uz = 0d0
      endif
c
c     Store the velocities in the ydot array.
c
      ydot(1) = ux
      ydot(2) = uz
c
      return
      end
c
c
c
      subroutine dummy
      implicit doubleprecision(a-h,o-z),integer(i-n)
      entry jacobn(n,t,y,dfdy,matdim,ml,mu)
      return
      end
c
c
c
c     --------------------------------------------------------------
c                          INVERSE ERROR FUNCTION
c              For a known X, compute Y such that ERF(Y) = X
c     --------------------------------------------------------------
      function erfi(x)
      doubleprecision x,erfi,keps,kconst,derfi,y0,y1,dy0,dy1
      integer kmaxit
c
      kmaxit = 50
      keps = 1d-14
      kconst = 0.8862269254527579d0 !sqrt(pi)/2.0
c
      if(abs(x).le.keps) then
          erfi = kconst * x
          return
      endif
      if(abs(x).lt.1d0) then
          erfi = kconst * abs(x)
          y0 = erf(0.9d0*erfi)
          derfi = 0.1d0 * erfi
          do iter=0,kmaxit
              y1 = 1d0 - erfc(erfi)
              dy1 = abs(x) - y1
              if (abs(dy1).lt.keps) then
                  if(x.lt.0d0) erfi = -erfi
                  return
              endif
              dy0 = y1 - y0
              derfi = derfi * dy1/dy0
              y0 = y1
              erfi = erfi + derfi
              if(abs(derfi/erfi).lt.keps) then
                  if(x.lt.0) erfi = -erfi
                  return
              endif
          enddo
      endif
      erfi = 0d0
      return
      end
c
c
C From Leonard J. Moss of SLAC:

C Here's a hybrid QuickSort I wrote a number of years ago.  It's
C based on suggestions in Knuth, Volume 3, and performs much better
C than a pure QuickSort on short or partially ordered input arrays.  

      SUBROUTINE SORTRX(N,DATA,INDEX)
C===================================================================
C
C     SORTRX -- SORT, Real input, indeX output
C
C
C     Input:  N     INTEGER
C             DATA  REAL
C
C     Output: INDEX INTEGER (DIMENSION N)
C
C This routine performs an in-memory sort of the first N elements of
C array DATA, returning into array INDEX the indices of elements of
C DATA arranged in ascending order.  Thus,
C
C    DATA(INDEX(1)) will be the smallest number in array DATA;
C    DATA(INDEX(N)) will be the largest number in DATA.
C
C The original data is not physically rearranged.  The original order
C of equal input values is not necessarily preserved.
C
C===================================================================
C
C SORTRX uses a hybrid QuickSort algorithm, based on several
C suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
C "pivot key" [my term] for dividing each subsequence is chosen to be
C the median of the first, last, and middle values of the subsequence;
C and the QuickSort is cut off when a subsequence has 9 or fewer
C elements, and a straight insertion sort of the entire array is done
C at the end.  The result is comparable to a pure insertion sort for
C very short arrays, and very fast for very large arrays (of order 12
C micro-sec/element on the 3081K for arrays of 10K elements).  It is
C also not subject to the poor performance of the pure QuickSort on
C partially ordered data.
C
C Created:  15 Jul 1986  Len Moss
C
C===================================================================
 
      INTEGER   N,INDEX(N)
      DOUBLEPRECISION DATA(N)
 
      INTEGER   LSTK(31),RSTK(31),ISTK
      INTEGER   L,R,I,J,P,INDEXP,INDEXT
      DOUBLEPRECISION DATAP
 
C     QuickSort Cutoff
C
C     Quit QuickSort-ing when a subsequence contains M or fewer
C     elements and finish off at end with straight insertion sort.
C     According to Knuth, V.3, the optimum value of M is around 9.
 
      INTEGER   M
      PARAMETER (M=9)
 
C===================================================================
C
C     Make initial guess for INDEX
 
      DO 50 I=1,N
         INDEX(I)=I
   50    CONTINUE
 
C     If array is short, skip QuickSort and go directly to
C     the straight insertion sort.
 
      IF (N.LE.M) GOTO 900
 
C===================================================================
C
C     QuickSort
C
C     The "Qn:"s correspond roughly to steps in Algorithm Q,
C     Knuth, V.3, PP.116-117, modified to select the median
C     of the first, last, and middle elements as the "pivot
C     key" (in Knuth's notation, "K").  Also modified to leave
C     data in place and produce an INDEX array.  To simplify
C     comments, let DATA[I]=DATA(INDEX(I)).
 
C Q1: Initialize
      ISTK=0
      L=1
      R=N
 
  200 CONTINUE
 
C Q2: Sort the subsequence DATA[L]..DATA[R].
C
C     At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
C     r > R, and L <= m <= R.  (First time through, there is no
C     DATA for l < L or r > R.)
 
      I=L
      J=R
 
C Q2.5: Select pivot key
C
C     Let the pivot, P, be the midpoint of this subsequence,
C     P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
C     so the corresponding DATA values are in increasing order.
C     The pivot key, DATAP, is then DATA[P].
 
      P=(L+R)/2
      INDEXP=INDEX(P)
      DATAP=DATA(INDEXP)
 
      IF (DATA(INDEX(L)) .GT. DATAP) THEN
         INDEX(P)=INDEX(L)
         INDEX(L)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
      IF (DATAP .GT. DATA(INDEX(R))) THEN
         IF (DATA(INDEX(L)) .GT. DATA(INDEX(R))) THEN
            INDEX(P)=INDEX(L)
            INDEX(L)=INDEX(R)
         ELSE
            INDEX(P)=INDEX(R)
         ENDIF
         INDEX(R)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
C     Now we swap values between the right and left sides and/or
C     move DATAP until all smaller values are on the left and all
C     larger values are on the right.  Neither the left or right
C     side will be internally ordered yet; however, DATAP will be
C     in its final position.
 
  300 CONTINUE
 
C Q3: Search for datum on left >= DATAP
C
C     At this point, DATA[L] <= DATAP.  We can therefore start scanning
C     up from L, looking for a value >= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         I=I+1
         IF (DATA(INDEX(I)).LT.DATAP) GOTO 300
 
  400 CONTINUE
 
C Q4: Search for datum on right <= DATAP
C
C     At this point, DATA[R] >= DATAP.  We can therefore start scanning
C     down from R, looking for a value <= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         J=J-1
         IF (DATA(INDEX(J)).GT.DATAP) GOTO 400
 
C Q5: Have the two scans collided?
 
      IF (I.LT.J) THEN
 
C Q6: No, interchange DATA[I] <--> DATA[J] and continue
 
         INDEXT=INDEX(I)
         INDEX(I)=INDEX(J)
         INDEX(J)=INDEXT
         GOTO 300
      ELSE
 
C Q7: Yes, select next subsequence to sort
C
C     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
C     for all L <= l < I and J < r <= R.  If both subsequences are
C     more than M elements long, push the longer one on the stack and
C     go back to QuickSort the shorter; if only one is more than M
C     elements long, go back and QuickSort it; otherwise, pop a
C     subsequence off the stack and QuickSort it.
 
         IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=J+1
            RSTK(ISTK)=R
            R=I-1
         ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=L
            RSTK(ISTK)=I-1
            L=J+1
         ELSE IF (R-J .GT. M) THEN
            L=J+1
         ELSE IF (I-L .GT. M) THEN
            R=I-1
         ELSE
C Q8: Pop the stack, or terminate QuickSort if empty
            IF (ISTK.LT.1) GOTO 900
            L=LSTK(ISTK)
            R=RSTK(ISTK)
            ISTK=ISTK-1
         ENDIF
         GOTO 200
      ENDIF
 
  900 CONTINUE
 
C===================================================================
C
C Q9: Straight Insertion sort
 
      DO 950 I=2,N
         IF (DATA(INDEX(I-1)) .GT. DATA(INDEX(I))) THEN
            INDEXP=INDEX(I)
            DATAP=DATA(INDEXP)
            P=I-1
  920       CONTINUE
               INDEX(P+1) = INDEX(P)
               P=P-1
               IF (P.GT.0) THEN
                  IF (DATA(INDEX(P)).GT.DATAP) GOTO 920
               ENDIF
            INDEX(P+1) = INDEXP
         ENDIF
  950    CONTINUE
 
C===================================================================
C
C     All done
 
      END
c
c     The Hull implementation of the Douglas and Peucker algorithm.
c     (c) Phil Wade, 1984
c     http://www2.dcs.hull.ac.uk/CISRG/publications/DPs/DP4/DP4.html
c
      SUBROUTINE GENERALIZE(XCOORD,YCOORD,PTSINCHAIN,GENCODE)
      IMPLICIT NONE
      INTEGER PTSINCHAIN,I,J,CORDSTART,CORDEND,POINT
      DOUBLEPRECISION XCOORD(PTSINCHAIN),YCOORD(PTSINCHAIN)
      DOUBLEPRECISION GENCODE(PTSINCHAIN)
      DOUBLEPRECISION DIS,MAXDIS
      
C     Set the generalization code of the first and last points to 999999
C     and set the others to -1.
      GENCODE(1)=1d90
      DO 30,I=2,PTSINCHAIN-1
          GENCODE(I)=-1
   30 CONTINUE
      GENCODE(PTSINCHAIN)=1d90
      CORDSTART=1
      
C     Main loop to give a generalization code to each point.
      DO 50,I=1,PTSINCHAIN-2
          
C     Find the first point on the line to have a generalization code.
          DO WHILE (GENCODE(CORDSTART+1) .GT. -1)
              CORDSTART=CORDSTART+1
          END DO
          
C     Find the first point after this with a generalization code.
          CORDEND=CORDSTART+2
          DO WHILE(GENCODE(CORDEND) .EQ. -1)
              CORDEND=CORDEND+1
          END DO 
          
C     Find the point with the maximum distance from the "Anchor-Floater"
C     line
          MAXDIS=-1.0
          DO 40,J=CORDSTART+1,CORDEND-1
              CALL CLCDIS(DIS,XCOORD(CORDSTART),YCOORD(CORDSTART),
     -            XCOORD(CORDEND),YCOORD(CORDEND),
     -            XCOORD(J),YCOORD(J))
              IF (DIS .GT. MAXDIS) THEN
                  POINT=J
                  MAXDIS=DIS
              END IF
   40     CONTINUE
          
C     Give the point which is furthest away from the "Anchor-Floater"
C     line the generalization code of this distance away.
          GENCODE(POINT)=MIN(GENCODE(CORDSTART),
     -        GENCODE(CORDEND),
     -        SQRT(MAXDIS))
   50 CONTINUE
      END
c
      SUBROUTINE CLCDIS (DIS,XCDST,YCDST,XCDEND,YCDEND,XPT,YPT)
      IMPLICIT NONE
      DOUBLEPRECISION XCDST,YCDST,XCDEND,YCDEND,XPT,YPT
      DOUBLEPRECISION X1,Y1,X2,Y2,X3,Y3,LAMBDA,X,Y,DIS
      
      X1=XCDST
      Y1=YCDST
      X2=XCDEND
      Y2=YCDEND
      X3=XPT
      Y3=YPT
      
      IF    (XCDST .EQ. XCDEND .AND. YCDST .EQ. YCDEND) THEN
          LAMBDA=0.0
      ELSE
          LAMBDA=(X1*(X1-X2-X3)+X2*X3+
     -        Y1*(Y1-Y2-Y3)+Y2*Y3)/
     -        ((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1))
      END IF
      
C     Calculate the values of x and y, this is the co-ordinate of the
C     point where the perpendicular from the point in question joins
C     the line from the anchor to the floater 
      
      IF      (LAMBDA .LT. 0.0) THEN
          X=X1
          Y=Y1
      ELSE IF (LAMBDA .GT. 1.0) THEN
          X=X2
          Y=Y2
      ELSE
          X=X1+LAMBDA*(X2-X1)
          Y=Y1+LAMBDA*(Y2-Y1)
      END IF 
      
C     Calculate the square of the distance from (X3,Y3) to (X,Y)
      DIS=(X3-X)*(X3-X)+(Y3-Y)*(Y3-Y)
      
      RETURN
      END
c
c     begin dvode
c
*DECK DVODE
      SUBROUTINE DVODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF,
     2            RPAR, IPAR)
      EXTERNAL F, JAC
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK, RPAR
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW,
     1        MF, IPAR
      DIMENSION Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW),
     1          RPAR(*), IPAR(*)
C-----------------------------------------------------------------------
C DVODE: Variable-coefficient Ordinary Differential Equation solver,
C with fixed-leading-coefficient implementation.
C This version is in double precision.
C
C DVODE solves the initial value problem for stiff or nonstiff
C systems of first order ODEs,
C     dy/dt = f(t,y) ,  or, in component form,
C     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
C DVODE is a package based on the EPISODE and EPISODEB packages, and
C on the ODEPACK user interface standard, with minor modifications.
C-----------------------------------------------------------------------
C Authors:
C               Peter N. Brown and Alan C. Hindmarsh
C               Center for Applied Scientific Computing, L-561
C               Lawrence Livermore National Laboratory
C               Livermore, CA 94551
C and
C               George D. Byrne
C               Illinois Institute of Technology
C               Chicago, IL 60616
C-----------------------------------------------------------------------
C References:
C
C 1. P. N. Brown, G. D. Byrne, and A. C. Hindmarsh, "VODE: A Variable
C    Coefficient ODE Solver," SIAM J. Sci. Stat. Comput., 10 (1989),
C    pp. 1038-1051.  Also, LLNL Report UCRL-98412, June 1988.
C 2. G. D. Byrne and A. C. Hindmarsh, "A Polyalgorithm for the
C    Numerical Solution of Ordinary Differential Equations,"
C    ACM Trans. Math. Software, 1 (1975), pp. 71-96.
C 3. A. C. Hindmarsh and G. D. Byrne, "EPISODE: An Effective Package
C    for the Integration of Systems of Ordinary Differential
C    Equations," LLNL Report UCID-30112, Rev. 1, April 1977.
C 4. G. D. Byrne and A. C. Hindmarsh, "EPISODEB: An Experimental
C    Package for the Integration of Systems of Ordinary Differential
C    Equations with Banded Jacobians," LLNL Report UCID-30132, April
C    1976.
C 5. A. C. Hindmarsh, "ODEPACK, a Systematized Collection of ODE
C    Solvers," in Scientific Computing, R. S. Stepleman et al., eds.,
C    North-Holland, Amsterdam, 1983, pp. 55-64.
C 6. K. R. Jackson and R. Sacks-Davis, "An Alternative Implementation
C    of Variable Step-Size Multistep Formulas for Stiff ODEs," ACM
C    Trans. Math. Software, 6 (1980), pp. 295-318.
C-----------------------------------------------------------------------
C Summary of usage.
C
C Communication between the user and the DVODE package, for normal
C situations, is summarized here.  This summary describes only a subset
C of the full set of options available.  See the full description for
C details, including optional communication, nonstandard options,
C and instructions for special situations.  See also the example
C problem (with program and output) following this summary.
C
C A. First provide a subroutine of the form:
C           SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR)
C           DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ), RPAR
C which supplies the vector function f by loading YDOT(i) with f(i).
C
C B. Next determine (or guess) whether or not the problem is stiff.
C Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
C whose real part is negative and large in magnitude, compared to the
C reciprocal of the t span of interest.  If the problem is nonstiff,
C use a method flag MF = 10.  If it is stiff, there are four standard
C choices for MF (21, 22, 24, 25), and DVODE requires the Jacobian
C matrix in some form.  In these cases (MF .gt. 0), DVODE will use a
C saved copy of the Jacobian matrix.  If this is undesirable because of
C storage limitations, set MF to the corresponding negative value
C (-21, -22, -24, -25).  (See full description of MF below.)
C The Jacobian matrix is regarded either as full (MF = 21 or 22),
C or banded (MF = 24 or 25).  In the banded case, DVODE requires two
C half-bandwidth parameters ML and MU.  These are, respectively, the
C widths of the lower and upper parts of the band, excluding the main
C diagonal.  Thus the band consists of the locations (i,j) with
C i-ML .le. j .le. i+MU, and the full bandwidth is ML+MU+1.
C
C C. If the problem is stiff, you are encouraged to supply the Jacobian
C directly (MF = 21 or 24), but if this is not feasible, DVODE will
C compute it internally by difference quotients (MF = 22 or 25).
C If you are supplying the Jacobian, provide a subroutine of the form:
C           SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
C           DOUBLE PRECISION T, Y(NEQ), PD(NROWPD,NEQ), RPAR
C which supplies df/dy by loading PD as follows:
C     For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),
C the partial derivative of f(i) with respect to y(j).  (Ignore the
C ML and MU arguments in this case.)
C     For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with
C df(i)/dy(j), i.e. load the diagonal lines of df/dy into the rows of
C PD from the top down.
C     In either case, only nonzero elements need be loaded.
C
C D. Write a main program which calls subroutine DVODE once for
C each point at which answers are desired.  This should also provide
C for possible use of logical unit 6 for output of error messages
C by DVODE.  On the first call to DVODE, supply arguments as follows:
C F      = Name of subroutine for right-hand side vector f.
C          This name must be declared external in calling program.
C NEQ    = Number of first order ODEs.
C Y      = Array of initial values, of length NEQ.
C T      = The initial value of the independent variable.
C TOUT   = First point where output is desired (.ne. T).
C ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
C RTOL   = Relative tolerance parameter (scalar).
C ATOL   = Absolute tolerance parameter (scalar or array).
C          The estimated local error in Y(i) will be controlled so as
C          to be roughly less (in magnitude) than
C             EWT(i) = RTOL*abs(Y(i)) + ATOL     if ITOL = 1, or
C             EWT(i) = RTOL*abs(Y(i)) + ATOL(i)  if ITOL = 2.
C          Thus the local error test passes if, in each component,
C          either the absolute error is less than ATOL (or ATOL(i)),
C          or the relative error is less than RTOL.
C          Use RTOL = 0.0 for pure absolute error control, and
C          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
C          control.  Caution: Actual (global) errors may exceed these
C          local tolerances, so choose them conservatively.
C ITASK  = 1 for normal computation of output values of Y at t = TOUT.
C ISTATE = Integer flag (input and output).  Set ISTATE = 1.
C IOPT   = 0 to indicate no optional input used.
C RWORK  = Real work array of length at least:
C             20 + 16*NEQ                      for MF = 10,
C             22 +  9*NEQ + 2*NEQ**2           for MF = 21 or 22,
C             22 + 11*NEQ + (3*ML + 2*MU)*NEQ  for MF = 24 or 25.
C LRW    = Declared length of RWORK (in user's DIMENSION statement).
C IWORK  = Integer work array of length at least:
C             30        for MF = 10,
C             30 + NEQ  for MF = 21, 22, 24, or 25.
C          If MF = 24 or 25, input in IWORK(1),IWORK(2) the lower
C          and upper half-bandwidths ML,MU.
C LIW    = Declared length of IWORK (in user's DIMENSION statement).
C JAC    = Name of subroutine for Jacobian matrix (MF = 21 or 24).
C          If used, this name must be declared external in calling
C          program.  If not used, pass a dummy name.
C MF     = Method flag.  Standard values are:
C          10 for nonstiff (Adams) method, no Jacobian used.
C          21 for stiff (BDF) method, user-supplied full Jacobian.
C          22 for stiff method, internally generated full Jacobian.
C          24 for stiff method, user-supplied banded Jacobian.
C          25 for stiff method, internally generated banded Jacobian.
C RPAR,IPAR = user-defined real and integer arrays passed to F and JAC.
C Note that the main program must declare arrays Y, RWORK, IWORK,
C and possibly ATOL, RPAR, and IPAR.
C
C E. The output from the first call (or any call) is:
C      Y = Array of computed values of y(t) vector.
C      T = Corresponding value of independent variable (normally TOUT).
C ISTATE = 2  if DVODE was successful, negative otherwise.
C          -1 means excess work done on this call. (Perhaps wrong MF.)
C          -2 means excess accuracy requested. (Tolerances too small.)
C          -3 means illegal input detected. (See printed message.)
C          -4 means repeated error test failures. (Check all input.)
C          -5 means repeated convergence failures. (Perhaps bad
C             Jacobian supplied or wrong choice of MF or tolerances.)
C          -6 means error weight became zero during problem. (Solution
C             component i vanished, and ATOL or ATOL(i) = 0.)
C
C F. To continue the integration after a successful return, simply
C reset TOUT and call DVODE again.  No other parameters need be reset.
C
C-----------------------------------------------------------------------
C EXAMPLE PROBLEM
C
C The following is a simple example problem, with the coding
C needed for its solution by DVODE.  The problem is from chemical
C kinetics, and consists of the following three rate equations:
C     dy1/dt = -.04*y1 + 1.e4*y2*y3
C     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
C     dy3/dt = 3.e7*y2**2
C on the interval from t = 0.0 to t = 4.e10, with initial conditions
C y1 = 1.0, y2 = y3 = 0.  The problem is stiff.
C
C The following coding solves this problem with DVODE, using MF = 21
C and printing results at t = .4, 4., ..., 4.e10.  It uses
C ITOL = 2 and ATOL much smaller for y2 than y1 or y3 because
C y2 has much smaller values.
C At the end of the run, statistical quantities of interest are
C printed. (See optional output in the full description below.)
C To generate Fortran source code, replace C in column 1 with a blank
C in the coding below.
C
C     EXTERNAL FEX, JEX
C     DOUBLE PRECISION ATOL, RPAR, RTOL, RWORK, T, TOUT, Y
C     DIMENSION Y(3), ATOL(3), RWORK(67), IWORK(33)
C     NEQ = 3
C     Y(1) = 1.0D0
C     Y(2) = 0.0D0
C     Y(3) = 0.0D0
C     T = 0.0D0
C     TOUT = 0.4D0
C     ITOL = 2
C     RTOL = 1.D-4
C     ATOL(1) = 1.D-8
C     ATOL(2) = 1.D-14
C     ATOL(3) = 1.D-6
C     ITASK = 1
C     ISTATE = 1
C     IOPT = 0
C     LRW = 67
C     LIW = 33
C     MF = 21
C     DO 40 IOUT = 1,12
C       CALL DVODE(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
C    1            IOPT,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR)
C       WRITE(6,20)T,Y(1),Y(2),Y(3)
C 20    FORMAT(' At t =',D12.4,'   y =',3D14.6)
C       IF (ISTATE .LT. 0) GO TO 80
C 40    TOUT = TOUT*10.
C     WRITE(6,60) IWORK(11),IWORK(12),IWORK(13),IWORK(19),
C    1            IWORK(20),IWORK(21),IWORK(22)
C 60  FORMAT(/' No. steps =',I4,'   No. f-s =',I4,
C    1       '   No. J-s =',I4,'   No. LU-s =',I4/
C    2       '  No. nonlinear iterations =',I4/
C    3       '  No. nonlinear convergence failures =',I4/
C    4       '  No. error test failures =',I4/)
C     STOP
C 80  WRITE(6,90)ISTATE
C 90  FORMAT(///' Error halt: ISTATE =',I3)
C     STOP
C     END
C
C     SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
C     DOUBLE PRECISION RPAR, T, Y, YDOT
C     DIMENSION Y(NEQ), YDOT(NEQ)
C     YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
C     YDOT(3) = 3.D7*Y(2)*Y(2)
C     YDOT(2) = -YDOT(1) - YDOT(3)
C     RETURN
C     END
C
C     SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
C     DOUBLE PRECISION PD, RPAR, T, Y
C     DIMENSION Y(NEQ), PD(NRPD,NEQ)
C     PD(1,1) = -.04D0
C     PD(1,2) = 1.D4*Y(3)
C     PD(1,3) = 1.D4*Y(2)
C     PD(2,1) = .04D0
C     PD(2,3) = -PD(1,3)
C     PD(3,2) = 6.D7*Y(2)
C     PD(2,2) = -PD(1,2) - PD(3,2)
C     RETURN
C     END
C
C The following output was obtained from the above program on a
C Cray-1 computer with the CFT compiler.
C
C At t =  4.0000e-01   y =  9.851680e-01  3.386314e-05  1.479817e-02
C At t =  4.0000e+00   y =  9.055255e-01  2.240539e-05  9.445214e-02
C At t =  4.0000e+01   y =  7.158108e-01  9.184883e-06  2.841800e-01
C At t =  4.0000e+02   y =  4.505032e-01  3.222940e-06  5.494936e-01
C At t =  4.0000e+03   y =  1.832053e-01  8.942690e-07  8.167938e-01
C At t =  4.0000e+04   y =  3.898560e-02  1.621875e-07  9.610142e-01
C At t =  4.0000e+05   y =  4.935882e-03  1.984013e-08  9.950641e-01
C At t =  4.0000e+06   y =  5.166183e-04  2.067528e-09  9.994834e-01
C At t =  4.0000e+07   y =  5.201214e-05  2.080593e-10  9.999480e-01
C At t =  4.0000e+08   y =  5.213149e-06  2.085271e-11  9.999948e-01
C At t =  4.0000e+09   y =  5.183495e-07  2.073399e-12  9.999995e-01
C At t =  4.0000e+10   y =  5.450996e-08  2.180399e-13  9.999999e-01
C
C No. steps = 595   No. f-s = 832   No. J-s =  13   No. LU-s = 112
C  No. nonlinear iterations = 831
C  No. nonlinear convergence failures =   0
C  No. error test failures =  22
C-----------------------------------------------------------------------
C Full description of user interface to DVODE.
C
C The user interface to DVODE consists of the following parts.
C
C i.   The call sequence to subroutine DVODE, which is a driver
C      routine for the solver.  This includes descriptions of both
C      the call sequence arguments and of user-supplied routines.
C      Following these descriptions is
C        * a description of optional input available through the
C          call sequence,
C        * a description of optional output (in the work arrays), and
C        * instructions for interrupting and restarting a solution.
C
C ii.  Descriptions of other routines in the DVODE package that may be
C      (optionally) called by the user.  These provide the ability to
C      alter error message handling, save and restore the internal
C      COMMON, and obtain specified derivatives of the solution y(t).
C
C iii. Descriptions of COMMON blocks to be declared in overlay
C      or similar environments.
C
C iv.  Description of two routines in the DVODE package, either of
C      which the user may replace with his own version, if desired.
C      these relate to the measurement of errors.
C
C-----------------------------------------------------------------------
C Part i.  Call Sequence.
C
C The call sequence parameters used for input only are
C     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF,
C and those used for both input and output are
C     Y, T, ISTATE.
C The work arrays RWORK and IWORK are also used for conditional and
C optional input and optional output.  (The term output here refers
C to the return from subroutine DVODE to the user's calling program.)
C
C The legality of input parameters will be thoroughly checked on the
C initial call for the problem, but not checked thereafter unless a
C change in input parameters is flagged by ISTATE = 3 in the input.
C
C The descriptions of the call arguments are as follows.
C
C F      = The name of the user-supplied subroutine defining the
C          ODE system.  The system must be put in the first-order
C          form dy/dt = f(t,y), where f is a vector-valued function
C          of the scalar t and the vector y.  Subroutine F is to
C          compute the function f.  It is to have the form
C               SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR)
C               DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ), RPAR
C          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
C          is output.  Y and YDOT are arrays of length NEQ.
C          Subroutine F should not alter Y(1),...,Y(NEQ).
C          F must be declared EXTERNAL in the calling program.
C
C          Subroutine F may access user-defined real and integer
C          work arrays RPAR and IPAR, which are to be dimensioned
C          in the main program.
C
C          If quantities computed in the F routine are needed
C          externally to DVODE, an extra call to F should be made
C          for this purpose, for consistent and accurate results.
C          If only the derivative dy/dt is needed, use DVINDY instead.
C
C NEQ    = The size of the ODE system (number of first order
C          ordinary differential equations).  Used only for input.
C          NEQ may not be increased during the problem, but
C          can be decreased (with ISTATE = 3 in the input).
C
C Y      = A real array for the vector of dependent variables, of
C          length NEQ or more.  Used for both input and output on the
C          first call (ISTATE = 1), and only for output on other calls.
C          On the first call, Y must contain the vector of initial
C          values.  In the output, Y contains the computed solution
C          evaluated at T.  If desired, the Y array may be used
C          for other purposes between calls to the solver.
C
C          This array is passed as the Y argument in all calls to
C          F and JAC.
C
C T      = The independent variable.  In the input, T is used only on
C          the first call, as the initial point of the integration.
C          In the output, after each call, T is the value at which a
C          computed solution Y is evaluated (usually the same as TOUT).
C          On an error return, T is the farthest point reached.
C
C TOUT   = The next value of t at which a computed solution is desired.
C          Used only for input.
C
C          When starting the problem (ISTATE = 1), TOUT may be equal
C          to T for one call, then should .ne. T for the next call.
C          For the initial T, an input value of TOUT .ne. T is used
C          in order to determine the direction of the integration
C          (i.e. the algebraic sign of the step sizes) and the rough
C          scale of the problem.  Integration in either direction
C          (forward or backward in t) is permitted.
C
C          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
C          the first call (i.e. the first call with TOUT .ne. T).
C          Otherwise, TOUT is required on every call.
C
C          If ITASK = 1, 3, or 4, the values of TOUT need not be
C          monotone, but a value of TOUT which backs up is limited
C          to the current internal t interval, whose endpoints are
C          TCUR - HU and TCUR.  (See optional output, below, for
C          TCUR and HU.)
C
C ITOL   = An indicator for the type of error control.  See
C          description below under ATOL.  Used only for input.
C
C RTOL   = A relative error tolerance parameter, either a scalar or
C          an array of length NEQ.  See description below under ATOL.
C          Input only.
C
C ATOL   = An absolute error tolerance parameter, either a scalar or
C          an array of length NEQ.  Input only.
C
C          The input parameters ITOL, RTOL, and ATOL determine
C          the error control performed by the solver.  The solver will
C          control the vector e = (e(i)) of estimated local errors
C          in Y, according to an inequality of the form
C                      rms-norm of ( e(i)/EWT(i) )   .le.   1,
C          where       EWT(i) = RTOL(i)*abs(Y(i)) + ATOL(i),
C          and the rms-norm (root-mean-square norm) here is
C          rms-norm(v) = sqrt(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
C          is a vector of weights which must always be positive, and
C          the values of RTOL and ATOL should all be non-negative.
C          The following table gives the types (scalar/array) of
C          RTOL and ATOL, and the corresponding form of EWT(i).
C
C             ITOL    RTOL       ATOL          EWT(i)
C              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
C              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
C              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
C              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
C
C          When either of these parameters is a scalar, it need not
C          be dimensioned in the user's calling program.
C
C          If none of the above choices (with ITOL, RTOL, and ATOL
C          fixed throughout the problem) is suitable, more general
C          error controls can be obtained by substituting
C          user-supplied routines for the setting of EWT and/or for
C          the norm calculation.  See Part iv below.
C
C          If global errors are to be estimated by making a repeated
C          run on the same problem with smaller tolerances, then all
C          components of RTOL and ATOL (i.e. of EWT) should be scaled
C          down uniformly.
C
C ITASK  = An index specifying the task to be performed.
C          Input only.  ITASK has the following values and meanings.
C          1  means normal computation of output values of y(t) at
C             t = TOUT (by overshooting and interpolating).
C          2  means take one step only and return.
C          3  means stop at the first internal mesh point at or
C             beyond t = TOUT and return.
C          4  means normal computation of output values of y(t) at
C             t = TOUT but without overshooting t = TCRIT.
C             TCRIT must be input as RWORK(1).  TCRIT may be equal to
C             or beyond TOUT, but not behind it in the direction of
C             integration.  This option is useful if the problem
C             has a singularity at or beyond t = TCRIT.
C          5  means take one step, without passing TCRIT, and return.
C             TCRIT must be input as RWORK(1).
C
C          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C          (within roundoff), it will return T = TCRIT (exactly) to
C          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
C          in which case answers at T = TOUT are returned first).
C
C ISTATE = an index used for input and output to specify the
C          the state of the calculation.
C
C          In the input, the values of ISTATE are as follows.
C          1  means this is the first call for the problem
C             (initializations will be done).  See note below.
C          2  means this is not the first call, and the calculation
C             is to continue normally, with no change in any input
C             parameters except possibly TOUT and ITASK.
C             (If ITOL, RTOL, and/or ATOL are changed between calls
C             with ISTATE = 2, the new values will be used but not
C             tested for legality.)
C          3  means this is not the first call, and the
C             calculation is to continue normally, but with
C             a change in input parameters other than
C             TOUT and ITASK.  Changes are allowed in
C             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF, ML, MU,
C             and any of the optional input except H0.
C             (See IWORK description for ML and MU.)
C          Note:  A preliminary call with TOUT = T is not counted
C          as a first call here, as no initialization or checking of
C          input is done.  (Such a call is sometimes useful to include
C          the initial conditions in the output.)
C          Thus the first call for which TOUT .ne. T requires
C          ISTATE = 1 in the input.
C
C          In the output, ISTATE has the following values and meanings.
C           1  means nothing was done, as TOUT was equal to T with
C              ISTATE = 1 in the input.
C           2  means the integration was performed successfully.
C          -1  means an excessive amount of work (more than MXSTEP
C              steps) was done on this call, before completing the
C              requested task, but the integration was otherwise
C              successful as far as T.  (MXSTEP is an optional input
C              and is normally 500.)  To continue, the user may
C              simply reset ISTATE to a value .gt. 1 and call again.
C              (The excess work step counter will be reset to 0.)
C              In addition, the user may increase MXSTEP to avoid
C              this error return.  (See optional input below.)
C          -2  means too much accuracy was requested for the precision
C              of the machine being used.  This was detected before
C              completing the requested task, but the integration
C              was successful as far as T.  To continue, the tolerance
C              parameters must be reset, and ISTATE must be set
C              to 3.  The optional output TOLSF may be used for this
C              purpose.  (Note: If this condition is detected before
C              taking any steps, then an illegal input return
C              (ISTATE = -3) occurs instead.)
C          -3  means illegal input was detected, before taking any
C              integration steps.  See written message for details.
C              Note:  If the solver detects an infinite loop of calls
C              to the solver with illegal input, it will cause
C              the run to stop.
C          -4  means there were repeated error test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              The problem may have a singularity, or the input
C              may be inappropriate.
C          -5  means there were repeated convergence test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              This may be caused by an inaccurate Jacobian matrix,
C              if one is being used.
C          -6  means EWT(i) became zero for some i during the
C              integration.  Pure relative error control (ATOL(i)=0.0)
C              was requested on a variable which has now vanished.
C              The integration was successful as far as T.
C
C          Note:  Since the normal output value of ISTATE is 2,
C          it does not need to be reset for normal continuation.
C          Also, since a negative input value of ISTATE will be
C          regarded as illegal, a negative output value requires the
C          user to change it, and possibly other input, before
C          calling the solver again.
C
C IOPT   = An integer flag to specify whether or not any optional
C          input is being used on this call.  Input only.
C          The optional input is listed separately below.
C          IOPT = 0 means no optional input is being used.
C                   Default values will be used in all cases.
C          IOPT = 1 means optional input is being used.
C
C RWORK  = A real working array (double precision).
C          The length of RWORK must be at least
C             20 + NYH*(MAXORD + 1) + 3*NEQ + LWM    where
C          NYH    = the initial value of NEQ,
C          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
C                   smaller value is given as an optional input),
C          LWM = length of work space for matrix-related data:
C          LWM = 0             if MITER = 0,
C          LWM = 2*NEQ**2 + 2  if MITER = 1 or 2, and MF.gt.0,
C          LWM = NEQ**2 + 2    if MITER = 1 or 2, and MF.lt.0,
C          LWM = NEQ + 2       if MITER = 3,
C          LWM = (3*ML+2*MU+2)*NEQ + 2 if MITER = 4 or 5, and MF.gt.0,
C          LWM = (2*ML+MU+1)*NEQ + 2   if MITER = 4 or 5, and MF.lt.0.
C          (See the MF description for METH and MITER.)
C          Thus if MAXORD has its default value and NEQ is constant,
C          this length is:
C             20 + 16*NEQ                    for MF = 10,
C             22 + 16*NEQ + 2*NEQ**2         for MF = 11 or 12,
C             22 + 16*NEQ + NEQ**2           for MF = -11 or -12,
C             22 + 17*NEQ                    for MF = 13,
C             22 + 18*NEQ + (3*ML+2*MU)*NEQ  for MF = 14 or 15,
C             22 + 17*NEQ + (2*ML+MU)*NEQ    for MF = -14 or -15,
C             20 +  9*NEQ                    for MF = 20,
C             22 +  9*NEQ + 2*NEQ**2         for MF = 21 or 22,
C             22 +  9*NEQ + NEQ**2           for MF = -21 or -22,
C             22 + 10*NEQ                    for MF = 23,
C             22 + 11*NEQ + (3*ML+2*MU)*NEQ  for MF = 24 or 25.
C             22 + 10*NEQ + (2*ML+MU)*NEQ    for MF = -24 or -25.
C          The first 20 words of RWORK are reserved for conditional
C          and optional input and optional output.
C
C          The following word in RWORK is a conditional input:
C            RWORK(1) = TCRIT = critical value of t which the solver
C                       is not to overshoot.  Required if ITASK is
C                       4 or 5, and ignored otherwise.  (See ITASK.)
C
C LRW    = The length of the array RWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C IWORK  = An integer work array.  The length of IWORK must be at least
C             30        if MITER = 0 or 3 (MF = 10, 13, 20, 23), or
C             30 + NEQ  otherwise (abs(MF) = 11,12,14,15,21,22,24,25).
C          The first 30 words of IWORK are reserved for conditional and
C          optional input and optional output.
C
C          The following 2 words in IWORK are conditional input:
C            IWORK(1) = ML     These are the lower and upper
C            IWORK(2) = MU     half-bandwidths, respectively, of the
C                       banded Jacobian, excluding the main diagonal.
C                       The band is defined by the matrix locations
C                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU
C                       must satisfy  0 .le.  ML,MU  .le. NEQ-1.
C                       These are required if MITER is 4 or 5, and
C                       ignored otherwise.  ML and MU may in fact be
C                       the band parameters for a matrix to which
C                       df/dy is only approximately equal.
C
C LIW    = the length of the array IWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C Note:  The work arrays must not be altered between calls to DVODE
C for the same problem, except possibly for the conditional and
C optional input, and except for the last 3*NEQ words of RWORK.
C The latter space is used for internal scratch space, and so is
C available for use by the user outside DVODE between calls, if
C desired (but not for use by F or JAC).
C
C JAC    = The name of the user-supplied routine (MITER = 1 or 4) to
C          compute the Jacobian matrix, df/dy, as a function of
C          the scalar t and the vector y.  It is to have the form
C               SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD,
C                               RPAR, IPAR)
C               DOUBLE PRECISION T, Y(NEQ), PD(NROWPD,NEQ), RPAR
C          where NEQ, T, Y, ML, MU, and NROWPD are input and the array
C          PD is to be loaded with partial derivatives (elements of the
C          Jacobian matrix) in the output.  PD must be given a first
C          dimension of NROWPD.  T and Y have the same meaning as in
C          Subroutine F.
C               In the full matrix case (MITER = 1), ML and MU are
C          ignored, and the Jacobian is to be loaded into PD in
C          columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
C               In the band matrix case (MITER = 4), the elements
C          within the band are to be loaded into PD in columnwise
C          manner, with diagonal lines of df/dy loaded into the rows
C          of PD. Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
C          ML and MU are the half-bandwidth parameters. (See IWORK).
C          The locations in PD in the two triangular areas which
C          correspond to nonexistent matrix elements can be ignored
C          or loaded arbitrarily, as they are overwritten by DVODE.
C               JAC need not provide df/dy exactly.  A crude
C          approximation (possibly with a smaller bandwidth) will do.
C               In either case, PD is preset to zero by the solver,
C          so that only the nonzero elements need be loaded by JAC.
C          Each call to JAC is preceded by a call to F with the same
C          arguments NEQ, T, and Y.  Thus to gain some efficiency,
C          intermediate quantities shared by both calculations may be
C          saved in a user COMMON block by F and not recomputed by JAC,
C          if desired.  Also, JAC may alter the Y array, if desired.
C          JAC must be declared external in the calling program.
C               Subroutine JAC may access user-defined real and integer
C          work arrays, RPAR and IPAR, whose dimensions are set by the
C          user in the main program.
C
C MF     = The method flag.  Used only for input.  The legal values of
C          MF are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, 25,
C          -11, -12, -14, -15, -21, -22, -24, -25.
C          MF is a signed two-digit integer, MF = JSV*(10*METH + MITER).
C          JSV = SIGN(MF) indicates the Jacobian-saving strategy:
C            JSV =  1 means a copy of the Jacobian is saved for reuse
C                     in the corrector iteration algorithm.
C            JSV = -1 means a copy of the Jacobian is not saved
C                     (valid only for MITER = 1, 2, 4, or 5).
C          METH indicates the basic linear multistep method:
C            METH = 1 means the implicit Adams method.
C            METH = 2 means the method based on backward
C                     differentiation formulas (BDF-s).
C          MITER indicates the corrector iteration method:
C            MITER = 0 means functional iteration (no Jacobian matrix
C                      is involved).
C            MITER = 1 means chord iteration with a user-supplied
C                      full (NEQ by NEQ) Jacobian.
C            MITER = 2 means chord iteration with an internally
C                      generated (difference quotient) full Jacobian
C                      (using NEQ extra calls to F per df/dy value).
C            MITER = 3 means chord iteration with an internally
C                      generated diagonal Jacobian approximation
C                      (using 1 extra call to F per df/dy evaluation).
C            MITER = 4 means chord iteration with a user-supplied
C                      banded Jacobian.
C            MITER = 5 means chord iteration with an internally
C                      generated banded Jacobian (using ML+MU+1 extra
C                      calls to F per df/dy evaluation).
C          If MITER = 1 or 4, the user must supply a subroutine JAC
C          (the name is arbitrary) as described above under JAC.
C          For other values of MITER, a dummy argument can be used.
C
C RPAR     User-specified array used to communicate real parameters
C          to user-supplied subroutines.  If RPAR is a vector, then
C          it must be dimensioned in the user's main program.  If it
C          is unused or it is a scalar, then it need not be
C          dimensioned.
C
C IPAR     User-specified array used to communicate integer parameter
C          to user-supplied subroutines.  The comments on dimensioning
C          RPAR apply to IPAR.
C-----------------------------------------------------------------------
C Optional Input.
C
C The following is a list of the optional input provided for in the
C call sequence.  (See also Part ii.)  For each such input variable,
C this table lists its name as used in this documentation, its
C location in the call sequence, its meaning, and the default value.
C The use of any of this input requires IOPT = 1, and in that
C case all of this input is examined.  A value of zero for any
C of these optional input variables will cause the default value to be
C used.  Thus to use a subset of the optional input, simply preload
C locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
C then set those of interest to nonzero values.
C
C NAME    LOCATION      MEANING AND DEFAULT VALUE
C
C H0      RWORK(5)  The step size to be attempted on the first step.
C                   The default value is determined by the solver.
C
C HMAX    RWORK(6)  The maximum absolute step size allowed.
C                   The default value is infinite.
C
C HMIN    RWORK(7)  The minimum absolute step size allowed.
C                   The default value is 0.  (This lower bound is not
C                   enforced on the final step before reaching TCRIT
C                   when ITASK = 4 or 5.)
C
C MAXORD  IWORK(5)  The maximum order to be allowed.  The default
C                   value is 12 if METH = 1, and 5 if METH = 2.
C                   If MAXORD exceeds the default value, it will
C                   be reduced to the default value.
C                   If MAXORD is changed during the problem, it may
C                   cause the current order to be reduced.
C
C MXSTEP  IWORK(6)  Maximum number of (internally defined) steps
C                   allowed during one call to the solver.
C                   The default value is 500.
C
C MXHNIL  IWORK(7)  Maximum number of messages printed (per problem)
C                   warning that T + H = T on a step (H = step size).
C                   This must be positive to result in a non-default
C                   value.  The default value is 10.
C
C-----------------------------------------------------------------------
C Optional Output.
C
C As optional additional output from DVODE, the variables listed
C below are quantities related to the performance of DVODE
C which are available to the user.  These are communicated by way of
C the work arrays, but also have internal mnemonic names as shown.
C Except where stated otherwise, all of this output is defined
C on any successful return from DVODE, and on any return with
C ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return
C (ISTATE = -3), they will be unchanged from their existing values
C (if any), except possibly for TOLSF, LENRW, and LENIW.
C On any error return, output relevant to the error will be defined,
C as noted below.
C
C NAME    LOCATION      MEANING
C
C HU      RWORK(11) The step size in t last used (successfully).
C
C HCUR    RWORK(12) The step size to be attempted on the next step.
C
C TCUR    RWORK(13) The current value of the independent variable
C                   which the solver has actually reached, i.e. the
C                   current internal mesh point in t.  In the output,
C                   TCUR will always be at least as far from the
C                   initial value of t as the current argument T,
C                   but may be farther (if interpolation was done).
C
C TOLSF   RWORK(14) A tolerance scale factor, greater than 1.0,
C                   computed when a request for too much accuracy was
C                   detected (ISTATE = -3 if detected at the start of
C                   the problem, ISTATE = -2 otherwise).  If ITOL is
C                   left unaltered but RTOL and ATOL are uniformly
C                   scaled up by a factor of TOLSF for the next call,
C                   then the solver is deemed likely to succeed.
C                   (The user may also ignore TOLSF and alter the
C                   tolerance parameters in any other way appropriate.)
C
C NST     IWORK(11) The number of steps taken for the problem so far.
C
C NFE     IWORK(12) The number of f evaluations for the problem so far.
C
C NJE     IWORK(13) The number of Jacobian evaluations so far.
C
C NQU     IWORK(14) The method order last used (successfully).
C
C NQCUR   IWORK(15) The order to be attempted on the next step.
C
C IMXER   IWORK(16) The index of the component of largest magnitude in
C                   the weighted local error vector ( e(i)/EWT(i) ),
C                   on an error return with ISTATE = -4 or -5.
C
C LENRW   IWORK(17) The length of RWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C LENIW   IWORK(18) The length of IWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C NLU     IWORK(19) The number of matrix LU decompositions so far.
C
C NNI     IWORK(20) The number of nonlinear (Newton) iterations so far.
C
C NCFN    IWORK(21) The number of convergence failures of the nonlinear
C                   solver so far.
C
C NETF    IWORK(22) The number of error test failures of the integrator
C                   so far.
C
C The following two arrays are segments of the RWORK array which
C may also be of interest to the user as optional output.
C For each array, the table below gives its internal name,
C its base address in RWORK, and its description.
C
C NAME    BASE ADDRESS      DESCRIPTION
C
C YH      21             The Nordsieck history array, of size NYH by
C                        (NQCUR + 1), where NYH is the initial value
C                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
C                        of YH contains HCUR**j/factorial(j) times
C                        the j-th derivative of the interpolating
C                        polynomial currently representing the
C                        solution, evaluated at t = TCUR.
C
C ACOR     LENRW-NEQ+1   Array of size NEQ used for the accumulated
C                        corrections on each step, scaled in the output
C                        to represent the estimated local error in Y
C                        on the last step.  This is the vector e in
C                        the description of the error control.  It is
C                        defined only on a successful return from DVODE.
C
C-----------------------------------------------------------------------
C Interrupting and Restarting
C
C If the integration of a given problem by DVODE is to be
C interrrupted and then later continued, such as when restarting
C an interrupted run or alternating between two or more ODE problems,
C the user should save, following the return from the last DVODE call
C prior to the interruption, the contents of the call sequence
C variables and internal COMMON blocks, and later restore these
C values before the next DVODE call for that problem.  To save
C and restore the COMMON blocks, use subroutine DVSRCO, as
C described below in part ii.
C
C In addition, if non-default values for either LUN or MFLAG are
C desired, an extra call to XSETUN and/or XSETF should be made just
C before continuing the integration.  See Part ii below for details.
C
C-----------------------------------------------------------------------
C Part ii.  Other Routines Callable.
C
C The following are optional calls which the user may make to
C gain additional capabilities in conjunction with DVODE.
C (The routines XSETUN and XSETF are designed to conform to the
C SLATEC error handling package.)
C
C     FORM OF CALL                  FUNCTION
C  CALL XSETUN(LUN)           Set the logical unit number, LUN, for
C                             output of messages from DVODE, if
C                             the default is not desired.
C                             The default value of LUN is 6.
C
C  CALL XSETF(MFLAG)          Set a flag to control the printing of
C                             messages by DVODE.
C                             MFLAG = 0 means do not print. (Danger:
C                             This risks losing valuable information.)
C                             MFLAG = 1 means print (the default).
C
C                             Either of the above calls may be made at
C                             any time and will take effect immediately.
C
C  CALL DVSRCO(RSAV,ISAV,JOB) Saves and restores the contents of
C                             the internal COMMON blocks used by
C                             DVODE. (See Part iii below.)
C                             RSAV must be a real array of length 49
C                             or more, and ISAV must be an integer
C                             array of length 40 or more.
C                             JOB=1 means save COMMON into RSAV/ISAV.
C                             JOB=2 means restore COMMON from RSAV/ISAV.
C                                DVSRCO is useful if one is
C                             interrupting a run and restarting
C                             later, or alternating between two or
C                             more problems solved with DVODE.
C
C  CALL DVINDY(,,,,,)         Provide derivatives of y, of various
C        (See below.)         orders, at a specified point T, if
C                             desired.  It may be called only after
C                             a successful return from DVODE.
C
C The detailed instructions for using DVINDY are as follows.
C The form of the call is:
C
C  CALL DVINDY (T, K, RWORK(21), NYH, DKY, IFLAG)
C
C The input parameters are:
C
C T         = Value of independent variable where answers are desired
C             (normally the same as the T last returned by DVODE).
C             For valid results, T must lie between TCUR - HU and TCUR.
C             (See optional output for TCUR and HU.)
C K         = Integer order of the derivative desired.  K must satisfy
C             0 .le. K .le. NQCUR, where NQCUR is the current order
C             (see optional output).  The capability corresponding
C             to K = 0, i.e. computing y(T), is already provided
C             by DVODE directly.  Since NQCUR .ge. 1, the first
C             derivative dy/dt is always available with DVINDY.
C RWORK(21) = The base address of the history array YH.
C NYH       = Column length of YH, equal to the initial value of NEQ.
C
C The output parameters are:
C
C DKY       = A real array of length NEQ containing the computed value
C             of the K-th derivative of y(t).
C IFLAG     = Integer flag, returned as 0 if K and T were legal,
C             -1 if K was illegal, and -2 if T was illegal.
C             On an error return, a message is also written.
C-----------------------------------------------------------------------
C Part iii.  COMMON Blocks.
C If DVODE is to be used in an overlay situation, the user
C must declare, in the primary overlay, the variables in:
C   (1) the call sequence to DVODE,
C   (2) the two internal COMMON blocks
C         /DVOD01/  of length  81  (48 double precision words
C                         followed by 33 integer words),
C         /DVOD02/  of length  9  (1 double precision word
C                         followed by 8 integer words),
C
C If DVODE is used on a system in which the contents of internal
C COMMON blocks are not preserved between calls, the user should
C declare the above two COMMON blocks in his main program to insure
C that their contents are preserved.
C
C-----------------------------------------------------------------------
C Part iv.  Optionally Replaceable Solver Routines.
C
C Below are descriptions of two routines in the DVODE package which
C relate to the measurement of errors.  Either routine can be
C replaced by a user-supplied version, if desired.  However, since such
C a replacement may have a major impact on performance, it should be
C done only when absolutely necessary, and only with great caution.
C (Note: The means by which the package version of a routine is
C superseded by the user's version may be system-dependent.)
C
C (a) DEWSET.
C The following subroutine is called just before each internal
C integration step, and sets the array of error weights, EWT, as
C described under ITOL/RTOL/ATOL above:
C     SUBROUTINE DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
C where NEQ, ITOL, RTOL, and ATOL are as in the DVODE call sequence,
C YCUR contains the current dependent variable vector, and
C EWT is the array of weights set by DEWSET.
C
C If the user supplies this subroutine, it must return in EWT(i)
C (i = 1,...,NEQ) a positive quantity suitable for comparison with
C errors in Y(i).  The EWT array returned by DEWSET is passed to the
C DVNORM routine (See below.), and also used by DVODE in the computation
C of the optional output IMXER, the diagonal Jacobian approximation,
C and the increments for difference quotient Jacobians.
C
C In the user-supplied version of DEWSET, it may be desirable to use
C the current values of derivatives of y.  Derivatives up to order NQ
C are available from the history array YH, described above under
C Optional Output.  In DEWSET, YH is identical to the YCUR array,
C extended to NQ + 1 columns with a column length of NYH and scale
C factors of h**j/factorial(j).  On the first call for the problem,
C given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
C NYH is the initial value of NEQ.  The quantities NQ, H, and NST
C can be obtained by including in DEWSET the statements:
C     DOUBLE PRECISION RVOD, H, HU
C     COMMON /DVOD01/ RVOD(48), IVOD(33)
C     COMMON /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C     NQ = IVOD(28)
C     H = RVOD(21)
C Thus, for example, the current value of dy/dt can be obtained as
C YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
C unnecessary when NST = 0).
C
C (b) DVNORM.
C The following is a real function routine which computes the weighted
C root-mean-square norm of a vector v:
C     D = DVNORM (N, V, W)
C where:
C   N = the length of the vector,
C   V = real array of length N containing the vector,
C   W = real array of length N containing weights,
C   D = sqrt( (1/N) * sum(V(i)*W(i))**2 ).
C DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
C EWT is as set by subroutine DEWSET.
C
C If the user supplies this function, it should return a non-negative
C value of DVNORM suitable for use in the error control in DVODE.
C None of the arguments should be altered by DVNORM.
C For example, a user-supplied DVNORM routine might:
C   -substitute a max-norm of (V(i)*W(i)) for the rms-norm, or
C   -ignore some components of V in the norm, with the effect of
C    suppressing the error control on those components of Y.
C-----------------------------------------------------------------------
C REVISION HISTORY (YYYYMMDD)
C  19890615  Date Written.  Initial release.
C  19890922  Added interrupt/restart ability, minor changes throughout.
C  19910228  Minor revisions in line format,  prologue, etc.
C  19920227  Modifications by D. Pang:
C            (1) Applied subgennam to get generic intrinsic names.
C            (2) Changed intrinsic names to generic in comments.
C            (3) Added *DECK lines before each routine.
C  19920721  Names of routines and labeled Common blocks changed, so as
C            to be unique in combined single/double precision code (ACH).
C  19920722  Minor revisions to prologue (ACH).
C  19920831  Conversion to double precision done (ACH).
C  19921106  Fixed minor bug: ETAQ,ETAQM1 in DVSTEP SAVE statement (ACH).
C  19921118  Changed LUNSAV/MFLGSV to IXSAV (ACH).
C  19941222  Removed MF overwrite; attached sign to H in estimated second 
C            deriv. in DVHIN; misc. comment changes throughout (ACH).
C  19970515  Minor corrections to comments in prologue, DVJAC (ACH).
C  19981111  Corrected Block B by adding final line, GO TO 200 (ACH).
C  20020430  Various upgrades (ACH): Use ODEPACK error handler package.
C            Replaced D1MACH by DUMACH.  Various changes to main
C            prologue and other routine prologues.
C-----------------------------------------------------------------------
C Other Routines in the DVODE Package.
C
C In addition to subroutine DVODE, the DVODE package includes the
C following subroutines and function routines:
C  DVHIN     computes an approximate step size for the initial step.
C  DVINDY    computes an interpolated value of the y vector at t = TOUT.
C  DVSTEP    is the core integrator, which does one step of the
C            integration and the associated error control.
C  DVSET     sets all method coefficients and test constants.
C  DVNLSD    solves the underlying nonlinear system -- the corrector.
C  DVJAC     computes and preprocesses the Jacobian matrix J = df/dy
C            and the Newton iteration matrix P = I - (h/l1)*J.
C  DVSOL     manages solution of linear system in chord iteration.
C  DVJUST    adjusts the history array on a change of order.
C  DEWSET    sets the error weight vector EWT before each step.
C  DVNORM    computes the weighted r.m.s. norm of a vector.
C  DVSRCO    is a user-callable routine to save and restore
C            the contents of the internal COMMON blocks.
C  DACOPY    is a routine to copy one two-dimensional array to another.
C  DGEFA and DGESL   are routines from LINPACK for solving full
C            systems of linear algebraic equations.
C  DGBFA and DGBSL   are routines from LINPACK for solving banded
C            linear systems.
C  DAXPY, DSCAL, and DCOPY are basic linear algebra modules (BLAS).
C  DUMACH    sets the unit roundoff of the machine.
C  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH handle the printing of all
C            error messages and warnings.  XERRWD is machine-dependent.
C Note:  DVNORM, DUMACH, IXSAV, and IUMACH are function routines.
C All the others are subroutines.
C
C-----------------------------------------------------------------------
C
C Type declarations for labeled COMMON block DVOD01 --------------------
C
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for labeled COMMON block DVOD02 --------------------
C
      DOUBLE PRECISION HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C
C Type declarations for local variables --------------------------------
C
      EXTERNAL DVNLSD
      LOGICAL IHIT
      DOUBLE PRECISION ATOLI, BIG, EWTI, FOUR, H0, HMAX, HMX, HUN, ONE,
     1   PT2, RH, RTOLI, SIZE, TCRIT, TNEXT, TOLSF, TP, TWO, ZERO
      INTEGER I, IER, IFLAG, IMXER, JCO, KGO, LENIW, LENJ, LENP, LENRW,
     1   LENWM, LF0, MBAND, MFA, ML, MORD, MU, MXHNL0, MXSTP0, NITER, 
     2   NSLAST
      CHARACTER*80 MSG
C
C Type declaration for function subroutines called ---------------------
C
      DOUBLE PRECISION DUMACH, DVNORM
C
      DIMENSION MORD(2)
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to DVODE.
C-----------------------------------------------------------------------
      SAVE MORD, MXHNL0, MXSTP0
      SAVE ZERO, ONE, TWO, FOUR, PT2, HUN
C-----------------------------------------------------------------------
C The following internal COMMON blocks contain variables which are
C communicated between subroutines in the DVODE package, or which are
C to be saved between calls to DVODE.
C In each block, real variables precede integers.
C The block /DVOD01/ appears in subroutines DVODE, DVINDY, DVSTEP,
C DVSET, DVNLSD, DVJAC, DVSOL, DVJUST and DVSRCO.
C The block /DVOD02/ appears in subroutines DVODE, DVINDY, DVSTEP,
C DVNLSD, DVJAC, and DVSRCO.
C
C The variables stored in the internal COMMON blocks are as follows:
C
C ACNRM  = Weighted r.m.s. norm of accumulated correction vectors.
C CCMXJ  = Threshhold on DRC for updating the Jacobian. (See DRC.)
C CONP   = The saved value of TQ(5).
C CRATE  = Estimated corrector convergence rate constant.
C DRC    = Relative change in H*RL1 since last DVJAC call.
C EL     = Real array of integration coefficients.  See DVSET.
C ETA    = Saved tentative ratio of new to old H.
C ETAMAX = Saved maximum value of ETA to be allowed.
C H      = The step size.
C HMIN   = The minimum absolute value of the step size H to be used.
C HMXI   = Inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C HNEW   = The step size to be attempted on the next step.
C HSCAL  = Stepsize in scaling of YH array.
C PRL1   = The saved value of RL1.
C RC     = Ratio of current H*RL1 to value on last DVJAC call.
C RL1    = The reciprocal of the coefficient EL(1).
C TAU    = Real vector of past NQ step sizes, length 13.
C TQ     = A real vector of length 5 in which DVSET stores constants
C          used for the convergence test, the error test, and the
C          selection of H at a new order.
C TN     = The independent variable, updated on each step taken.
C UROUND = The machine unit roundoff.  The smallest positive real number
C          such that  1.0 + UROUND .ne. 1.0
C ICF    = Integer flag for convergence failure in DVNLSD:
C            0 means no failures.
C            1 means convergence failure with out of date Jacobian
C                   (recoverable error).
C            2 means convergence failure with current Jacobian or
C                   singular matrix (unrecoverable error).
C INIT   = Saved integer flag indicating whether initialization of the
C          problem has been done (INIT = 1) or not.
C IPUP   = Saved flag to signal updating of Newton matrix.
C JCUR   = Output flag from DVJAC showing Jacobian status:
C            JCUR = 0 means J is not current.
C            JCUR = 1 means J is current.
C JSTART = Integer flag used as input to DVSTEP:
C            0  means perform the first step.
C            1  means take a new step continuing from the last.
C            -1 means take the next step with a new value of MAXORD,
C                  HMIN, HMXI, N, METH, MITER, and/or matrix parameters.
C          On return, DVSTEP sets JSTART = 1.
C JSV    = Integer flag for Jacobian saving, = sign(MF).
C KFLAG  = A completion code from DVSTEP with the following meanings:
C               0      the step was succesful.
C              -1      the requested error could not be achieved.
C              -2      corrector convergence could not be achieved.
C              -3, -4  fatal error in VNLS (can not occur here).
C KUTH   = Input flag to DVSTEP showing whether H was reduced by the
C          driver.  KUTH = 1 if H was reduced, = 0 otherwise.
C L      = Integer variable, NQ + 1, current order plus one.
C LMAX   = MAXORD + 1 (used for dimensioning).
C LOCJS  = A pointer to the saved Jacobian, whose storage starts at
C          WM(LOCJS), if JSV = 1.
C LYH, LEWT, LACOR, LSAVF, LWM, LIWM = Saved integer pointers
C          to segments of RWORK and IWORK.
C MAXORD = The maximum order of integration method to be allowed.
C METH/MITER = The method flags.  See MF.
C MSBJ   = The maximum number of steps between J evaluations, = 50.
C MXHNIL = Saved value of optional input MXHNIL.
C MXSTEP = Saved value of optional input MXSTEP.
C N      = The number of first-order ODEs, = NEQ.
C NEWH   = Saved integer to flag change of H.
C NEWQ   = The method order to be used on the next step.
C NHNIL  = Saved counter for occurrences of T + H = T.
C NQ     = Integer variable, the current integration method order.
C NQNYH  = Saved value of NQ*NYH.
C NQWAIT = A counter controlling the frequency of order changes.
C          An order change is about to be considered if NQWAIT = 1.
C NSLJ   = The number of steps taken as of the last Jacobian update.
C NSLP   = Saved value of NST as of last Newton matrix update.
C NYH    = Saved value of the initial value of NEQ.
C HU     = The step size in t last used.
C NCFN   = Number of nonlinear convergence failures so far.
C NETF   = The number of error test failures of the integrator so far.
C NFE    = The number of f evaluations for the problem so far.
C NJE    = The number of Jacobian evaluations so far.
C NLU    = The number of matrix LU decompositions so far.
C NNI    = Number of nonlinear iterations so far.
C NQU    = The method order last used.
C NST    = The number of steps taken for the problem so far.
C-----------------------------------------------------------------------
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C
      DATA  MORD(1) /12/, MORD(2) /5/, MXSTP0 /500/, MXHNL0 /10/
      DATA ZERO /0.0D0/, ONE /1.0D0/, TWO /2.0D0/, FOUR /4.0D0/,
     1     PT2 /0.2D0/, HUN /100.0D0/
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .gt. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .NE. 1) GO TO 603
c
c     Next line added by RV on 13-Feb-2009 to allow for the possibility
c     of changing Y outside DVODE.
c
      CALL DCOPY (N, Y, 1, RWORK(LYH), 1)
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 1),
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all input and various initializations.
C
C First check legality of the non-optional input NEQ, ITOL, IOPT,
C MF, ML, and MU.
C-----------------------------------------------------------------------
 20   IF (NEQ .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ .GT. N) GO TO 605
 25   N = NEQ
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      JSV = SIGN(1,MF)
      MFA = ABS(MF)
      METH = MFA/10
      MITER = MFA - 10*METH
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 0 .OR. MITER .GT. 5) GO TO 608
      IF (MITER .LE. 3) GO TO 30
      ML = IWORK(1)
      MU = IWORK(2)
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610
 30   CONTINUE
C Next process and check the optional input. ---------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .EQ. 1) H0 = ZERO
      HMXI = ZERO
      HMIN = ZERO
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. ZERO) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. ZERO) GO TO 615
      HMXI = ZERO
      IF (HMAX .GT. ZERO) HMXI = ONE/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. ZERO) GO TO 616
C-----------------------------------------------------------------------
C Set work array pointers and check lengths LRW and LIW.
C Pointers to segments of RWORK and IWORK are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
C Within WM, LOCJS is the location of the saved Jacobian (JSV .gt. 0).
C-----------------------------------------------------------------------
 60   LYH = 21
c      write(*,*) 'hmin: ',HMIN
      IF (ISTATE .EQ. 1) NYH = N
      LWM = LYH + (MAXORD + 1)*NYH
      JCO = MAX(0,JSV)
      IF (MITER .EQ. 0) LENWM = 0
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        LENWM = 2 + (1 + JCO)*N*N
        LOCJS = N*N + 3
      ENDIF
      IF (MITER .EQ. 3) LENWM = 2 + N
      IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        MBAND = ML + MU + 1
        LENP = (MBAND + ML)*N
        LENJ = MBAND*N
        LENWM = 2 + LENP + JCO*LENJ
        LOCJS = LENP + 3
        ENDIF
      LEWT = LWM + LENWM
      LSAVF = LEWT + N
      LACOR = LSAVF + N
      LENRW = LACOR + N - 1
      IWORK(17) = LENRW
      LIWM = 1
      LENIW = 30 + N
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) LENIW = 30
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
C Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 70 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. ZERO) GO TO 619
        IF (ATOLI .LT. ZERO) GO TO 620
 70     CONTINUE
      IF (ISTATE .EQ. 1) GO TO 100
C If ISTATE = 3, set flag to signal parameter changes to DVSTEP. -------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 90
C MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
      CALL DCOPY (N, RWORK(LWM), 1, RWORK(LSAVF), 1)
C Reload WM(1) = RWORK(LWM), since LWM may have changed. ---------------
 90   IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 1).
C It contains all remaining initializations, the initial call to F,
C and the calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  UROUND = DUMACH()
      TN = T
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. ZERO) GO TO 625
      IF (H0 .NE. ZERO .AND. (T + H0 - TCRIT)*H0 .GT. ZERO)
     1   H0 = TCRIT - T
 110  JSTART = 0
      IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
      CCMXJ = PT2
      MSBJ = 50
      NHNIL = 0
      NST = 0
      NJE = 0
      NNI = 0
      NCFN = 0
      NETF = 0
      NLU = 0
      NSLJ = 0
      NSLAST = 0
      HU = ZERO
      NQU = 0
C Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      CALL F (N, T, Y, RWORK(LF0), RPAR, IPAR)
c      write(*,*) 'y(13) 21: ',y(13)
      NFE = 1
C Load the initial value vector in YH. ---------------------------------
      CALL DCOPY (N, Y, 1, RWORK(LYH), 1)
c      write(*,*) 'y(13) 22: ',y(13)
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1
      H = ONE
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 120 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. ZERO) GO TO 621
 120    RWORK(I+LEWT-1) = ONE/RWORK(I+LEWT-1)
      IF (H0 .NE. ZERO) GO TO 180
C Call DVHIN to set initial step size H0 to be attempted. --------------
      CALL DVHIN (N, T, RWORK(LYH), RWORK(LF0), F, RPAR, IPAR, TOUT,
     1   UROUND, RWORK(LEWT), ITOL, ATOL, Y, RWORK(LACOR), H0,
     2   NITER, IER)
c      write(*,*) 'y(13) 23: ',y(13)
      NFE = NFE + NITER
      IF (IER .NE. 0) GO TO 622
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. ONE) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      CALL DSCAL (N, H0, RWORK(LF0), 1)
      GO TO 270
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C-----------------------------------------------------------------------
 200  NSLAST = NST
      KUTH = 0
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      CALL DVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
c      write(*,*) 'y(13) 24: ',y(13)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(ONE + HUN*UROUND)
      IF ((TP - TOUT)*H .GT. ZERO) GO TO 623
      IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. ZERO) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. ZERO) GO TO 625
      IF ((TN - TOUT)*H .LT. ZERO) GO TO 245
      CALL DVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
c      write(*,*) 'y(13) 25: ',y(13)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. ZERO) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + HNEW*(ONE + FOUR*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. ZERO) GO TO 250
      H = (TCRIT - TN)*(ONE - FOUR*UROUND)
      KUTH = 1
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator DVSTEP.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken, update EWT (if not at
C start of problem), check for too much accuracy being requested, and
C check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. ZERO) GO TO 510
 260    RWORK(I+LEWT-1) = ONE/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*DVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. ONE) GO TO 280
      TOLSF = TOLSF*TWO
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'DVODE--  Warning: internal T (=R1) and H (=R2) are'
      CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWD (MSG, 60, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      (H = step size). solver will continue anyway'
      CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'DVODE--  Above warning has been issued I1 times.  '
      CALL XERRWD (MSG, 50, 102, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      it will not be issued again for this problem'
      CALL XERRWD (MSG, 50, 102, 1, 1, MXHNIL, 0, 0, ZERO, ZERO)
 290  CONTINUE
C-----------------------------------------------------------------------
C CALL DVSTEP (Y, YH, NYH, YH, EWT, SAVF, VSAV, ACOR,
C              WM, IWM, F, JAC, F, DVNLSD, RPAR, IPAR)
C-----------------------------------------------------------------------
      CALL DVSTEP (Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   RWORK(LSAVF), Y, RWORK(LACOR), RWORK(LWM), IWORK(LIWM),
     2   F, JAC, F, DVNLSD, RPAR, IPAR)
c      write(*,*) 'y(13) 26: ',y(13)
      KGO = 1 - KFLAG
c      write(*,*) 'kflag 3',kflag
C Branch on KFLAG.  Note: In this version, KFLAG can not be set to -3.
C  KFLAG .eq. 0,   -1,  -2
c      write(*,*) 'y(13) 530: ',y(13)
      GO TO (300, 530, 540), KGO
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).  Test for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
      KUTH = 0
      GO TO (310, 400, 330, 340, 350), ITASK
C ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      CALL DVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
c      write(*,*) 'y(13) 27: ',y(13)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. ZERO) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. ZERO) GO TO 345
      CALL DVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
c      write(*,*) 'y(13) 28: ',y(13)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + HNEW*(ONE + FOUR*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. ZERO) GO TO 250
      H = (TCRIT - TN)*(ONE - FOUR*UROUND)
      KUTH = 1
      GO TO 250
C ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from DVODE.
C If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional output is loaded into the work
C arrays before returning.
C-----------------------------------------------------------------------
 400  CONTINUE
      CALL DCOPY (N, RWORK(LYH), 1, Y, 1)
c      write(*,*) 'y(13) 29: ',y(13)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      RWORK(11) = HU
      RWORK(12) = HNEW
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NEWQ
      IWORK(19) = NLU
      IWORK(20) = NNI
      IWORK(21) = NCFN
      IWORK(22) = NETF
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C if there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH, and T is set to TN.
C The optional output is loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'DVODE--  At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWD (MSG, 50, 201, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWD (MSG, 50, 201, 1, 1, MXSTEP, 0, 1, TN, ZERO)
      ISTATE = -1
      GO TO 580
C EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'DVODE--  At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWD (MSG, 50, 202, 1, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'DVODE--  At T (=R1), too much accuracy requested  '
      CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      for precision of machine:   see TOLSF (=R2) '
      CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'DVODE--  At T(=R1) and step size H(=R2), the error'
      CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      test failed repeatedly or with abs(H) = HMIN'
      CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'DVODE--  At T (=R1) and step size H (=R2), the    '
      CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      or with abs(H) = HMIN   '
      CALL XERRWD (MSG, 30, 205, 1, 0, 0, 0, 2, TN, H)
      ISTATE = -5
C Compute IMXER if relevant. -------------------------------------------
 560  BIG = ZERO
      IMXER = 1
      DO 570 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
C Set Y vector, T, and optional output. --------------------------------
 580  CONTINUE
      CALL DCOPY (N, RWORK(LYH), 1, Y, 1)
c      write(*,*) 'y(13) 30: ',y(13)
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = NLU
      IWORK(20) = NNI
      IWORK(21) = NCFN
      IWORK(22) = NETF
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.   If the illegal input
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'DVODE--  ISTATE (=I1) illegal '
      CALL XERRWD (MSG, 30, 1, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'DVODE--  ITASK (=I1) illegal  '
      CALL XERRWD (MSG, 30, 2, 1, 1, ITASK, 0, 0, ZERO, ZERO)
      GO TO 700
 603  MSG='DVODE--  ISTATE (=I1) .gt. 1 but DVODE not initialized      '
      CALL XERRWD (MSG, 60, 3, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
      GO TO 700
 604  MSG = 'DVODE--  NEQ (=I1) .lt. 1     '
      CALL XERRWD (MSG, 30, 4, 1, 1, NEQ, 0, 0, ZERO, ZERO)
      GO TO 700
 605  MSG = 'DVODE--  ISTATE = 3 and NEQ increased (I1 to I2)  '
      CALL XERRWD (MSG, 50, 5, 1, 2, N, NEQ, 0, ZERO, ZERO)
      GO TO 700
 606  MSG = 'DVODE--  ITOL (=I1) illegal   '
      CALL XERRWD (MSG, 30, 6, 1, 1, ITOL, 0, 0, ZERO, ZERO)
      GO TO 700
 607  MSG = 'DVODE--  IOPT (=I1) illegal   '
      CALL XERRWD (MSG, 30, 7, 1, 1, IOPT, 0, 0, ZERO, ZERO)
      GO TO 700
 608  MSG = 'DVODE--  MF (=I1) illegal     '
      CALL XERRWD (MSG, 30, 8, 1, 1, MF, 0, 0, ZERO, ZERO)
      GO TO 700
 609  MSG = 'DVODE--  ML (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)'
      CALL XERRWD (MSG, 50, 9, 1, 2, ML, NEQ, 0, ZERO, ZERO)
      GO TO 700
 610  MSG = 'DVODE--  MU (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)'
      CALL XERRWD (MSG, 50, 10, 1, 2, MU, NEQ, 0, ZERO, ZERO)
      GO TO 700
 611  MSG = 'DVODE--  MAXORD (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 11, 1, 1, MAXORD, 0, 0, ZERO, ZERO)
      GO TO 700
 612  MSG = 'DVODE--  MXSTEP (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 12, 1, 1, MXSTEP, 0, 0, ZERO, ZERO)
      GO TO 700
 613  MSG = 'DVODE--  MXHNIL (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 13, 1, 1, MXHNIL, 0, 0, ZERO, ZERO)
      GO TO 700
 614  MSG = 'DVODE--  TOUT (=R1) behind T (=R2)      '
      CALL XERRWD (MSG, 40, 14, 1, 0, 0, 0, 2, TOUT, T)
      MSG = '      integration direction is given by H0 (=R1)  '
      CALL XERRWD (MSG, 50, 14, 1, 0, 0, 0, 1, H0, ZERO)
      GO TO 700
 615  MSG = 'DVODE--  HMAX (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 15, 1, 0, 0, 0, 1, HMAX, ZERO)
      GO TO 700
 616  MSG = 'DVODE--  HMIN (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 16, 1, 0, 0, 0, 1, HMIN, ZERO)
      GO TO 700
 617  CONTINUE
      MSG='DVODE--  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWD (MSG, 60, 17, 1, 2, LENRW, LRW, 0, ZERO, ZERO)
      GO TO 700
 618  CONTINUE
      MSG='DVODE--  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWD (MSG, 60, 18, 1, 2, LENIW, LIW, 0, ZERO, ZERO)
      GO TO 700
 619  MSG = 'DVODE--  RTOL(I1) is R1 .lt. 0.0        '
      CALL XERRWD (MSG, 40, 19, 1, 1, I, 0, 1, RTOLI, ZERO)
      GO TO 700
 620  MSG = 'DVODE--  ATOL(I1) is R1 .lt. 0.0        '
      CALL XERRWD (MSG, 40, 20, 1, 1, I, 0, 1, ATOLI, ZERO)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'DVODE--  EWT(I1) is R1 .le. 0.0         '
      CALL XERRWD (MSG, 40, 21, 1, 1, I, 0, 1, EWTI, ZERO)
      GO TO 700
 622  CONTINUE
      MSG='DVODE--  TOUT (=R1) too close to T(=R2) to start integration'
      CALL XERRWD (MSG, 60, 22, 1, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  CONTINUE
      MSG='DVODE--  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWD (MSG, 60, 23, 1, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  CONTINUE
      MSG='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWD (MSG, 60, 24, 1, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  CONTINUE
      MSG='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWD (MSG, 60, 25, 1, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'DVODE--  At start of problem, too much accuracy   '
      CALL XERRWD (MSG, 50, 26, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG='      requested for precision of machine:   see TOLSF (=R1) '
      CALL XERRWD (MSG, 60, 26, 1, 0, 0, 0, 1, TOLSF, ZERO)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG='DVODE--  Trouble from DVINDY.  ITASK = I1, TOUT = R1.       '
      CALL XERRWD (MSG, 60, 27, 1, 1, ITASK, 0, 1, TOUT, ZERO)
C
 700  CONTINUE
      ISTATE = -3
      RETURN
C
 800  MSG = 'DVODE--  Run aborted:  apparent infinite loop     '
      CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, ZERO, ZERO)
      RETURN
C----------------------- End of Subroutine DVODE -----------------------
      END
*DECK DVHIN
      SUBROUTINE DVHIN (N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND,
     1   EWT, ITOL, ATOL, Y, TEMP, H0, NITER, IER)
      EXTERNAL F
      DOUBLE PRECISION T0, Y0, YDOT, RPAR, TOUT, UROUND, EWT, ATOL, Y,
     1   TEMP, H0
      INTEGER N, IPAR, ITOL, NITER, IER
      DIMENSION Y0(*), YDOT(*), EWT(*), ATOL(*), Y(*),
     1   TEMP(*), RPAR(*), IPAR(*)
C-----------------------------------------------------------------------
C Call sequence input -- N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND,
C                        EWT, ITOL, ATOL, Y, TEMP
C Call sequence output -- H0, NITER, IER
C COMMON block variables accessed -- None
C
C Subroutines called by DVHIN:  F
C Function routines called by DVHI: DVNORM
C-----------------------------------------------------------------------
C This routine computes the step size, H0, to be attempted on the
C first step, when the user has not supplied a value for this.
C
C First we check that TOUT - T0 differs significantly from zero.  Then
C an iteration is done to approximate the initial second derivative
C and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.
C A bias factor of 1/2 is applied to the resulting h.
C The sign of H0 is inferred from the initial values of TOUT and T0.
C
C Communication with DVHIN is done with the following variables:
C
C N      = Size of ODE system, input.
C T0     = Initial value of independent variable, input.
C Y0     = Vector of initial conditions, input.
C YDOT   = Vector of initial first derivatives, input.
C F      = Name of subroutine for right-hand side f(t,y), input.
C RPAR, IPAR = Dummy names for user's real and integer work arrays.
C TOUT   = First output value of independent variable
C UROUND = Machine unit roundoff
C EWT, ITOL, ATOL = Error weights and tolerance parameters
C                   as described in the driver routine, input.
C Y, TEMP = Work arrays of length N.
C H0     = Step size to be attempted, output.
C NITER  = Number of iterations (and of f evaluations) to compute H0,
C          output.
C IER    = The error flag, returned with the value
C          IER = 0  if no trouble occurred, or
C          IER = -1 if TOUT and T0 are considered too close to proceed.
C-----------------------------------------------------------------------
C
C Type declarations for local variables --------------------------------
C
      DOUBLE PRECISION AFI, ATOLI, DELYI, H, HALF, HG, HLB, HNEW, HRAT,
     1     HUB, HUN, PT1, T1, TDIST, TROUND, TWO, YDDNRM
      INTEGER I, ITER
C
C Type declaration for function subroutines called ---------------------
C
      DOUBLE PRECISION DVNORM
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE HALF, HUN, PT1, TWO
      DATA HALF /0.5D0/, HUN /100.0D0/, PT1 /0.1D0/, TWO /2.0D0/
C
      NITER = 0
      TDIST = ABS(TOUT - T0)
      TROUND = UROUND*MAX(ABS(T0),ABS(TOUT))
      IF (TDIST .LT. TWO*TROUND) GO TO 100
C
C Set a lower bound on h based on the roundoff level in T0 and TOUT. ---
      HLB = HUN*TROUND
C Set an upper bound on h based on TOUT-T0 and the initial Y and YDOT. -
      HUB = PT1*TDIST
      ATOLI = ATOL(1)
      DO 10 I = 1, N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        DELYI = PT1*ABS(Y0(I)) + ATOLI
        AFI = ABS(YDOT(I))
        IF (AFI*HUB .GT. DELYI) HUB = DELYI/AFI
 10     CONTINUE
C
C Set initial guess for h as geometric mean of upper and lower bounds. -
      ITER = 0
      HG = SQRT(HLB*HUB)
C If the bounds have crossed, exit with the mean value. ----------------
      IF (HUB .LT. HLB) THEN
        H0 = HG
        GO TO 90
      ENDIF
C
C Looping point for iteration. -----------------------------------------
 50   CONTINUE
C Estimate the second derivative as a difference quotient in f. --------
      H = SIGN (HG, TOUT - T0)
      T1 = T0 + H
      DO 60 I = 1, N
 60     Y(I) = Y0(I) + H*YDOT(I)
c      write(*,*) 'y(13): ',y(13)
      CALL F (N, T1, Y, TEMP, RPAR, IPAR)
      DO 70 I = 1, N
 70     TEMP(I) = (TEMP(I) - YDOT(I))/H
      YDDNRM = DVNORM (N, TEMP, EWT)
C Get the corresponding new value of h. --------------------------------
      IF (YDDNRM*HUB*HUB .GT. TWO) THEN
        HNEW = SQRT(TWO/YDDNRM)
      ELSE
        HNEW = SQRT(HG*HUB)
      ENDIF
      ITER = ITER + 1
C-----------------------------------------------------------------------
C Test the stopping conditions.
C Stop if the new and previous h values differ by a factor of .lt. 2.
C Stop if four iterations have been done.  Also, stop with previous h
C if HNEW/HG .gt. 2 after first iteration, as this probably means that
C the second derivative value is bad because of cancellation error.
C-----------------------------------------------------------------------
      IF (ITER .GE. 4) GO TO 80
      HRAT = HNEW/HG
      IF ( (HRAT .GT. HALF) .AND. (HRAT .LT. TWO) ) GO TO 80
      IF ( (ITER .GE. 2) .AND. (HNEW .GT. TWO*HG) ) THEN
        HNEW = HG
        GO TO 80
      ENDIF
      HG = HNEW
      GO TO 50
C
C Iteration done.  Apply bounds, bias factor, and sign.  Then exit. ----
 80   H0 = HNEW*HALF
      IF (H0 .LT. HLB) H0 = HLB
      IF (H0 .GT. HUB) H0 = HUB
 90   H0 = SIGN(H0, TOUT - T0)
      NITER = ITER
      IER = 0
      RETURN
C Error return for TOUT - T0 too small. --------------------------------
 100  IER = -1
      RETURN
C----------------------- End of Subroutine DVHIN -----------------------
      END
*DECK DVINDY
      SUBROUTINE DVINDY (T, K, YH, LDYH, DKY, IFLAG)
      DOUBLE PRECISION T, YH, DKY
      INTEGER K, LDYH, IFLAG
      DIMENSION YH(LDYH,*), DKY(*)
C-----------------------------------------------------------------------
C Call sequence input -- T, K, YH, LDYH
C Call sequence output -- DKY, IFLAG
C COMMON block variables accessed:
C     /DVOD01/ --  H, TN, UROUND, L, N, NQ
C     /DVOD02/ --  HU
C
C Subroutines called by DVINDY: DSCAL, XERRWD
C Function routines called by DVINDY: None
C-----------------------------------------------------------------------
C DVINDY computes interpolated values of the K-th derivative of the
C dependent variable vector y, and stores it in DKY.  This routine
C is called within the package with K = 0 and T = TOUT, but may
C also be called by the user for any K up to the current order.
C (See detailed instructions in the usage documentation.)
C-----------------------------------------------------------------------
C The computed values in DKY are gotten by interpolation using the
C Nordsieck history array YH.  This array corresponds uniquely to a
C vector-valued polynomial of degree NQCUR or less, and DKY is set
C to the K-th derivative of this polynomial at T.
C The formula for DKY is:
C              q
C  DKY(i)  =  sum  c(j,K) * (T - TN)**(j-K) * H**(-j) * YH(i,j+1)
C             j=K
C where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, TN = TCUR, H = HCUR.
C The quantities  NQ = NQCUR, L = NQ+1, N, TN, and H are
C communicated by COMMON.  The above sum is done in reverse order.
C IFLAG is returned negative if either K or T is out of bounds.
C
C Discussion above and comments in driver explain all variables.
C-----------------------------------------------------------------------
C
C Type declarations for labeled COMMON block DVOD01 --------------------
C
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for labeled COMMON block DVOD02 --------------------
C
      DOUBLE PRECISION HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C
C Type declarations for local variables --------------------------------
C
      DOUBLE PRECISION C, HUN, R, S, TFUZZ, TN1, TP, ZERO
      INTEGER I, IC, J, JB, JB2, JJ, JJ1, JP1
      CHARACTER*80 MSG
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE HUN, ZERO
C
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C
      DATA HUN /100.0D0/, ZERO /0.0D0/
C
      IFLAG = 0
      IF (K .LT. 0 .OR. K .GT. NQ) GO TO 80
      TFUZZ = HUN*UROUND*SIGN(ABS(TN) + ABS(HU), HU)
      TP = TN - HU - TFUZZ
      TN1 = TN + TFUZZ
      IF ((T-TP)*(T-TN1) .GT. ZERO) GO TO 90
C
      S = (T - TN)/H
      IC = 1
      IF (K .EQ. 0) GO TO 15
      JJ1 = L - K
      DO 10 JJ = JJ1, NQ
 10     IC = IC*JJ
 15   C = REAL(IC)
      DO 20 I = 1, N
 20     DKY(I) = C*YH(I,L)
      IF (K .EQ. NQ) GO TO 55
      JB2 = NQ - K
      DO 50 JB = 1, JB2
        J = NQ - JB
        JP1 = J + 1
        IC = 1
        IF (K .EQ. 0) GO TO 35
        JJ1 = JP1 - K
        DO 30 JJ = JJ1, J
 30       IC = IC*JJ
 35     C = REAL(IC)
        DO 40 I = 1, N
 40       DKY(I) = C*YH(I,JP1) + S*DKY(I)
 50     CONTINUE
      IF (K .EQ. 0) RETURN
 55   R = H**(-K)
      CALL DSCAL (N, R, DKY, 1)
      RETURN
C
 80   MSG = 'DVINDY-- K (=I1) illegal      '
      CALL XERRWD (MSG, 30, 51, 1, 1, K, 0, 0, ZERO, ZERO)
      IFLAG = -1
      RETURN
 90   MSG = 'DVINDY-- T (=R1) illegal      '
      CALL XERRWD (MSG, 30, 52, 1, 0, 0, 0, 1, T, ZERO)
      MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
      CALL XERRWD (MSG, 60, 52, 1, 0, 0, 0, 2, TP, TN)
      IFLAG = -2
      RETURN
C----------------------- End of Subroutine DVINDY ----------------------
      END
*DECK DVSTEP
      SUBROUTINE DVSTEP (Y, YH, LDYH, YH1, EWT, SAVF, VSAV, ACOR,
     1                  WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR)
      EXTERNAL F, JAC, PSOL, VNLS
      DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, VSAV, ACOR, WM, RPAR
      INTEGER LDYH, IWM, IPAR
      DIMENSION Y(*), YH(LDYH,*), YH1(*), EWT(*), SAVF(*), VSAV(*),
     1   ACOR(*), WM(*), IWM(*), RPAR(*), IPAR(*)
C-----------------------------------------------------------------------
C Call sequence input -- Y, YH, LDYH, YH1, EWT, SAVF, VSAV,
C                        ACOR, WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR
C Call sequence output -- YH, ACOR, WM, IWM
C COMMON block variables accessed:
C     /DVOD01/  ACNRM, EL(13), H, HMIN, HMXI, HNEW, HSCAL, RC, TAU(13),
C               TQ(5), TN, JCUR, JSTART, KFLAG, KUTH,
C               L, LMAX, MAXORD, N, NEWQ, NQ, NQWAIT
C     /DVOD02/  HU, NCFN, NETF, NFE, NQU, NST
C
C Subroutines called by DVSTEP: F, DAXPY, DCOPY, DSCAL,
C                               DVJUST, VNLS, DVSET
C Function routines called by DVSTEP: DVNORM
C-----------------------------------------------------------------------
C DVSTEP performs one step of the integration of an initial value
C problem for a system of ordinary differential equations.
C DVSTEP calls subroutine VNLS for the solution of the nonlinear system
C arising in the time step.  Thus it is independent of the problem
C Jacobian structure and the type of nonlinear system solution method.
C DVSTEP returns a completion flag KFLAG (in COMMON).
C A return with KFLAG = -1 or -2 means either ABS(H) = HMIN or 10
C consecutive failures occurred.  On a return with KFLAG negative,
C the values of TN and the YH array are as of the beginning of the last
C step, and H is the last step size attempted.
C
C Communication with DVSTEP is done with the following variables:
C
C Y      = An array of length N used for the dependent variable vector.
C YH     = An LDYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C          j-th derivative of y(i), scaled by H**j/factorial(j)
C          (j = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C LDYH   = A constant integer .ge. N, the first dimension of YH.
C          N is the number of ODEs in the system.
C YH1    = A one-dimensional array occupying the same space as YH.
C EWT    = An array of length N containing multiplicative weights
C          for local error measurements.  Local errors in y(i) are
C          compared to 1.0/EWT(i) in various error tests.
C SAVF   = An array of working storage, of length N.
C          also used for input of YH(*,MAXORD+2) when JSTART = -1
C          and MAXORD .lt. the current order NQ.
C VSAV   = A work array of length N passed to subroutine VNLS.
C ACOR   = A work array of length N, used for the accumulated
C          corrections.  On a successful return, ACOR(i) contains
C          the estimated one-step local error in y(i).
C WM,IWM = Real and integer work arrays associated with matrix
C          operations in VNLS.
C F      = Dummy name for the user supplied subroutine for f.
C JAC    = Dummy name for the user supplied Jacobian subroutine.
C PSOL   = Dummy name for the subroutine passed to VNLS, for
C          possible use there.
C VNLS   = Dummy name for the nonlinear system solving subroutine,
C          whose real name is dependent on the method used.
C RPAR, IPAR = Dummy names for user's real and integer work arrays.
C-----------------------------------------------------------------------
C
C Type declarations for labeled COMMON block DVOD01 --------------------
C
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for labeled COMMON block DVOD02 --------------------
C
      DOUBLE PRECISION HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C
C Type declarations for local variables --------------------------------
C
      DOUBLE PRECISION ADDON, BIAS1,BIAS2,BIAS3, CNQUOT, DDN, DSM, DUP,
     1     ETACF, ETAMIN, ETAMX1, ETAMX2, ETAMX3, ETAMXF,
     2     ETAQ, ETAQM1, ETAQP1, FLOTL, ONE, ONEPSM,
     3     R, THRESH, TOLD, ZERO
      INTEGER I, I1, I2, IBACK, J, JB, KFC, KFH, MXNCF, NCF, NFLAG
C
C Type declaration for function subroutines called ---------------------
C
      DOUBLE PRECISION DVNORM
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE ADDON, BIAS1, BIAS2, BIAS3,
     1     ETACF, ETAMIN, ETAMX1, ETAMX2, ETAMX3, ETAMXF, ETAQ, ETAQM1,
     2     KFC, KFH, MXNCF, ONEPSM, THRESH, ONE, ZERO
C-----------------------------------------------------------------------
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C
      DATA KFC/-3/, KFH/-70/, MXNCF/10/
      DATA ADDON  /1.0D-6/,    BIAS1  /6.0D0/,     BIAS2  /6.0D0/,
     1     BIAS3  /10.0D0/,    ETACF  /0.25D0/,    ETAMIN /0.1D0/,
     2     ETAMXF /0.2D0/,     ETAMX1 /1.0D4/,     ETAMX2 /10.0D0/,
     3     ETAMX3 /10.0D0/,    ONEPSM /1.00001D0/, THRESH /1.5D0/
      DATA ONE/1.0D0/, ZERO/0.0D0/
C
      KFLAG = 0
      TOLD = TN
      NCF = 0
      JCUR = 0
      NFLAG = 0
      IF (JSTART .GT. 0) GO TO 20
      IF (JSTART .EQ. -1) GO TO 100
C-----------------------------------------------------------------------
C On the first call, the order is set to 1, and other variables are
C initialized.  ETAMAX is the maximum ratio by which H can be increased
C in a single step.  It is normally 10, but is larger during the
C first step to compensate for the small initial H.  If a failure
C occurs (in corrector convergence or error test), ETAMAX is set to 1
C for the next increase.
C-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      NQNYH = NQ*LDYH
      TAU(1) = H
      PRL1 = ONE
      RC = ZERO
      ETAMAX = ETAMX1
      NQWAIT = 2
      HSCAL = H
      GO TO 200
C-----------------------------------------------------------------------
C Take preliminary actions on a normal continuation step (JSTART.GT.0).
C If the driver changed H, then ETA must be reset and NEWH set to 1.
C If a change of order was dictated on the previous step, then
C it is done here and appropriate adjustments in the history are made.
C On an order decrease, the history array is adjusted by DVJUST.
C On an order increase, the history array is augmented by a column.
C On a change of step size H, the history array YH is rescaled.
C-----------------------------------------------------------------------
 20   CONTINUE
      IF (KUTH .EQ. 1) THEN
        ETA = MIN(ETA,H/HSCAL)
        NEWH = 1
        ENDIF
 50   IF (NEWH .EQ. 0) GO TO 200
      IF (NEWQ .EQ. NQ) GO TO 150
      IF (NEWQ .LT. NQ) THEN
        CALL DVJUST (YH, LDYH, -1)
        NQ = NEWQ
        L = NQ + 1
        NQWAIT = L
        GO TO 150
        ENDIF
      IF (NEWQ .GT. NQ) THEN
        CALL DVJUST (YH, LDYH, 1)
        NQ = NEWQ
        L = NQ + 1
        NQWAIT = L
        GO TO 150
      ENDIF
C-----------------------------------------------------------------------
C The following block handles preliminaries needed when JSTART = -1.
C If N was reduced, zero out part of YH to avoid undefined references.
C If MAXORD was reduced to a value less than the tentative order NEWQ,
C then NQ is set to MAXORD, and a new H ratio ETA is chosen.
C Otherwise, we take the same preliminary actions as for JSTART .gt. 0.
C In any case, NQWAIT is reset to L = NQ + 1 to prevent further
C changes in order for that many steps.
C The new H ratio ETA is limited by the input H if KUTH = 1,
C by HMIN if KUTH = 0, and by HMXI in any case.
C Finally, the history array YH is rescaled.
C-----------------------------------------------------------------------
 100  CONTINUE
      LMAX = MAXORD + 1
      IF (N .EQ. LDYH) GO TO 120
      I1 = 1 + (NEWQ + 1)*LDYH
      I2 = (MAXORD + 1)*LDYH
      IF (I1 .GT. I2) GO TO 120
      DO 110 I = I1, I2
 110    YH1(I) = ZERO
 120  IF (NEWQ .LE. MAXORD) GO TO 140
      FLOTL = REAL(LMAX)
      IF (MAXORD .LT. NQ-1) THEN
        DDN = DVNORM (N, SAVF, EWT)/TQ(1)
        ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
        ENDIF
      IF (MAXORD .EQ. NQ .AND. NEWQ .EQ. NQ+1) ETA = ETAQ
      IF (MAXORD .EQ. NQ-1 .AND. NEWQ .EQ. NQ+1) THEN
        ETA = ETAQM1
        CALL DVJUST (YH, LDYH, -1)
        ENDIF
      IF (MAXORD .EQ. NQ-1 .AND. NEWQ .EQ. NQ) THEN
        DDN = DVNORM (N, SAVF, EWT)/TQ(1)
        ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
        CALL DVJUST (YH, LDYH, -1)
        ENDIF
      ETA = MIN(ETA,ONE)
      NQ = MAXORD
      L = LMAX
 140  IF (KUTH .EQ. 1) ETA = MIN(ETA,ABS(H/HSCAL))
      IF (KUTH .EQ. 0) ETA = MAX(ETA,HMIN/ABS(HSCAL))
      ETA = ETA/MAX(ONE,ABS(HSCAL)*HMXI*ETA)
      NEWH = 1
      NQWAIT = L
      IF (NEWQ .LE. MAXORD) GO TO 50
C Rescale the history array for a change in H by a factor of ETA. ------
 150  R = ONE
      DO 180 J = 2, L
        R = R*ETA
        CALL DSCAL (N, R, YH(1,J), 1 )
 180    CONTINUE
      H = HSCAL*ETA
      HSCAL = H
      RC = RC*ETA
      NQNYH = NQ*LDYH
C-----------------------------------------------------------------------
C This section computes the predicted values by effectively
C multiplying the YH array by the Pascal triangle matrix.
C DVSET is called to calculate all integration coefficients.
C RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
C-----------------------------------------------------------------------
 200  TN = TN + H
      I1 = NQNYH + 1
      DO 220 JB = 1, NQ
        I1 = I1 - LDYH
        DO 210 I = I1, NQNYH
 210      YH1(I) = YH1(I) + YH1(I+LDYH)
 220  CONTINUE
      CALL DVSET
      RL1 = ONE/EL(2)
      RC = RC*(RL1/PRL1)
      PRL1 = RL1
C
C Call the nonlinear system solver. ------------------------------------
C
      CALL VNLS (Y, YH, LDYH, VSAV, SAVF, EWT, ACOR, IWM, WM,
     1           F, JAC, PSOL, NFLAG, RPAR, IPAR)
c      write(*,*) 'y(13) 31: ',y(13)
C
      IF (NFLAG .EQ. 0) GO TO 450
C-----------------------------------------------------------------------
C The VNLS routine failed to achieve convergence (NFLAG .NE. 0).
C The YH array is retracted to its values before prediction.
C The step size H is reduced and the step is retried, if possible.
C Otherwise, an error exit is taken.
C-----------------------------------------------------------------------
        NCF = NCF + 1
        NCFN = NCFN + 1
        ETAMAX = ONE
        TN = TOLD
        I1 = NQNYH + 1
        DO 430 JB = 1, NQ
          I1 = I1 - LDYH
          DO 420 I = I1, NQNYH
 420        YH1(I) = YH1(I) - YH1(I+LDYH)
 430      CONTINUE
        IF (NFLAG .LT. -1) GO TO 680
        IF (ABS(H) .LE. HMIN*ONEPSM) GO TO 670
        IF (NCF .EQ. MXNCF) GO TO 670
        ETA = ETACF
        ETA = MAX(ETA,HMIN/ABS(H))
        NFLAG = -1
        GO TO 150
C-----------------------------------------------------------------------
C The corrector has converged (NFLAG = 0).  The local error test is
C made and control passes to statement 500 if it fails.
C-----------------------------------------------------------------------
 450  CONTINUE
      DSM = ACNRM/TQ(2)
      IF (DSM .GT. ONE) GO TO 500
C-----------------------------------------------------------------------
C After a successful step, update the YH and TAU arrays and decrement
C NQWAIT.  If NQWAIT is then 1 and NQ .lt. MAXORD, then ACOR is saved
C for use in a possible order increase on the next step.
C If ETAMAX = 1 (a failure occurred this step), keep NQWAIT .ge. 2.
C-----------------------------------------------------------------------
      KFLAG = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      DO 470 IBACK = 1, NQ
        I = L - IBACK
 470    TAU(I+1) = TAU(I)
      TAU(1) = H
      DO 480 J = 1, L
        CALL DAXPY (N, EL(J), ACOR, 1, YH(1,J), 1 )
 480    CONTINUE
      NQWAIT = NQWAIT - 1
      IF ((L .EQ. LMAX) .OR. (NQWAIT .NE. 1)) GO TO 490
      CALL DCOPY (N, ACOR, 1, YH(1,LMAX), 1 )
      CONP = TQ(5)
 490  IF (ETAMAX .NE. ONE) GO TO 560
      IF (NQWAIT .LT. 2) NQWAIT = 2
      NEWQ = NQ
      NEWH = 0
      ETA = ONE
      HNEW = H
      GO TO 690
C-----------------------------------------------------------------------
C The error test failed.  KFLAG keeps track of multiple failures.
C Restore TN and the YH array to their previous values, and prepare
C to try the step again.  Compute the optimum step size for the
C same order.  After repeated failures, H is forced to decrease
C more rapidly.
C-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
c      write(*,*) 'kflag 1',kflag
      NETF = NETF + 1
      NFLAG = -2
      TN = TOLD
      I1 = NQNYH + 1
      DO 520 JB = 1, NQ
        I1 = I1 - LDYH
        DO 510 I = I1, NQNYH
 510      YH1(I) = YH1(I) - YH1(I+LDYH)
 520  CONTINUE
      IF (ABS(H) .LE. HMIN*ONEPSM) GO TO 660
      ETAMAX = ONE
      IF (KFLAG .LE. KFC) GO TO 530
C Compute ratio of new H to current H at the current order. ------------
      FLOTL = REAL(L)
      ETA = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
      ETA = MAX(ETA,HMIN/ABS(H),ETAMIN)
      IF ((KFLAG .LE. -2) .AND. (ETA .GT. ETAMXF)) ETA = ETAMXF
      GO TO 150
C-----------------------------------------------------------------------
C Control reaches this section if 3 or more consecutive failures
C have occurred.  It is assumed that the elements of the YH array
C have accumulated errors of the wrong order.  The order is reduced
C by one, if possible.  Then H is reduced by a factor of 0.1 and
C the step is retried.  After a total of 70 consecutive failures,
C an exit is taken with KFLAG = -1.
C-----------------------------------------------------------------------
 530  IF (KFLAG .EQ. KFH) GO TO 660
c      write(*,*) 'kflag 4',kflag
      IF (NQ .EQ. 1) GO TO 540
      ETA = MAX(ETAMIN,HMIN/ABS(H))
      CALL DVJUST (YH, LDYH, -1)
      L = NQ
      NQ = NQ - 1
      NQWAIT = L
      GO TO 150
 540  ETA = MAX(ETAMIN,HMIN/ABS(H))
      H = H*ETA
      HSCAL = H
      TAU(1) = H
      CALL F (N, TN, Y, SAVF, RPAR, IPAR)
c      write(*,*) 'y(13) 32: ',y(13)
      NFE = NFE + 1
      DO 550 I = 1, N
 550    YH(I,2) = H*SAVF(I)
      NQWAIT = 10
      GO TO 200
C-----------------------------------------------------------------------
C If NQWAIT = 0, an increase or decrease in order by one is considered.
C Factors ETAQ, ETAQM1, ETAQP1 are computed by which H could
C be multiplied at order q, q-1, or q+1, respectively.
C The largest of these is determined, and the new order and
C step size set accordingly.
C A change of H or NQ is made only if H increases by at least a
C factor of THRESH.  If an order change is considered and rejected,
C then NQWAIT is set to 2 (reconsider it after 2 steps).
C-----------------------------------------------------------------------
C Compute ratio of new H to current H at the current order. ------------
 560  FLOTL = REAL(L)
      ETAQ = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
      IF (NQWAIT .NE. 0) GO TO 600
      NQWAIT = 2
      ETAQM1 = ZERO
      IF (NQ .EQ. 1) GO TO 570
C Compute ratio of new H to current H at the current order less one. ---
      DDN = DVNORM (N, YH(1,L), EWT)/TQ(1)
      ETAQM1 = ONE/((BIAS1*DDN)**(ONE/(FLOTL - ONE)) + ADDON)
 570  ETAQP1 = ZERO
      IF (L .EQ. LMAX) GO TO 580
C Compute ratio of new H to current H at current order plus one. -------
      CNQUOT = (TQ(5)/CONP)*(H/TAU(2))**L
      DO 575 I = 1, N
 575    SAVF(I) = ACOR(I) - CNQUOT*YH(I,LMAX)
      DUP = DVNORM (N, SAVF, EWT)/TQ(3)
      ETAQP1 = ONE/((BIAS3*DUP)**(ONE/(FLOTL + ONE)) + ADDON)
 580  IF (ETAQ .GE. ETAQP1) GO TO 590
      IF (ETAQP1 .GT. ETAQM1) GO TO 620
      GO TO 610
 590  IF (ETAQ .LT. ETAQM1) GO TO 610
 600  ETA = ETAQ
      NEWQ = NQ
      GO TO 630
 610  ETA = ETAQM1
      NEWQ = NQ - 1
      GO TO 630
 620  ETA = ETAQP1
      NEWQ = NQ + 1
      CALL DCOPY (N, ACOR, 1, YH(1,LMAX), 1)
C Test tentative new H against THRESH, ETAMAX, and HMXI, then exit. ----
 630  IF (ETA .LT. THRESH .OR. ETAMAX .EQ. ONE) GO TO 640
      ETA = MIN(ETA,ETAMAX)
      ETA = ETA/MAX(ONE,ABS(H)*HMXI*ETA)
      NEWH = 1
      HNEW = H*ETA
      GO TO 690
 640  NEWQ = NQ
      NEWH = 0
      ETA = ONE
      HNEW = H
      GO TO 690
C-----------------------------------------------------------------------
C All returns are made through this section.
C On a successful return, ETAMAX is reset and ACOR is scaled.
C-----------------------------------------------------------------------
 660  KFLAG = -1
c      write(*,*) 'kflag 2',kflag
      GO TO 720
 670  KFLAG = -2
c      write(*,*) 'kflag 5',kflag
      GO TO 720
 680  IF (NFLAG .EQ. -2) KFLAG = -3
      IF (NFLAG .EQ. -3) KFLAG = -4
c      write(*,*) 'kflag 6',kflag
      GO TO 720
 690  ETAMAX = ETAMX3
      IF (NST .LE. 10) ETAMAX = ETAMX2
 700  R = ONE/TQ(2)
      CALL DSCAL (N, R, ACOR, 1)
 720  JSTART = 1
      RETURN
C----------------------- End of Subroutine DVSTEP ----------------------
      END
*DECK DVSET
      SUBROUTINE DVSET
C-----------------------------------------------------------------------
C Call sequence communication: None
C COMMON block variables accessed:
C     /DVOD01/ -- EL(13), H, TAU(13), TQ(5), L(= NQ + 1),
C                 METH, NQ, NQWAIT
C
C Subroutines called by DVSET: None
C Function routines called by DVSET: None
C-----------------------------------------------------------------------
C DVSET is called by DVSTEP and sets coefficients for use there.
C
C For each order NQ, the coefficients in EL are calculated by use of
C  the generating polynomial lambda(x), with coefficients EL(i).
C      lambda(x) = EL(1) + EL(2)*x + ... + EL(NQ+1)*(x**NQ).
C For the backward differentiation formulas,
C                                     NQ-1
C      lambda(x) = (1 + x/xi*(NQ)) * product (1 + x/xi(i) ) .
C                                     i = 1
C For the Adams formulas,
C                              NQ-1
C      (d/dx) lambda(x) = c * product (1 + x/xi(i) ) ,
C                              i = 1
C      lambda(-1) = 0,    lambda(0) = 1,
C where c is a normalization constant.
C In both cases, xi(i) is defined by
C      H*xi(i) = t sub n  -  t sub (n-i)
C              = H + TAU(1) + TAU(2) + ... TAU(i-1).
C
C
C In addition to variables described previously, communication
C with DVSET uses the following:
C   TAU    = A vector of length 13 containing the past NQ values
C            of H.
C   EL     = A vector of length 13 in which vset stores the
C            coefficients for the corrector formula.
C   TQ     = A vector of length 5 in which vset stores constants
C            used for the convergence test, the error test, and the
C            selection of H at a new order.
C   METH   = The basic method indicator.
C   NQ     = The current order.
C   L      = NQ + 1, the length of the vector stored in EL, and
C            the number of columns of the YH array being used.
C   NQWAIT = A counter controlling the frequency of order changes.
C            An order change is about to be considered if NQWAIT = 1.
C-----------------------------------------------------------------------
C
C Type declarations for labeled COMMON block DVOD01 --------------------
C
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for local variables --------------------------------
C
      DOUBLE PRECISION AHATN0, ALPH0, CNQM1, CORTES, CSUM, ELP, EM,
     1     EM0, FLOTI, FLOTL, FLOTNQ, HSUM, ONE, RXI, RXIS, S, SIX,
     2     T1, T2, T3, T4, T5, T6, TWO, XI, ZERO
      INTEGER I, IBACK, J, JP1, NQM1, NQM2
C
      DIMENSION EM(13)
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE CORTES, ONE, SIX, TWO, ZERO
C
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
C
      DATA CORTES /0.1D0/
      DATA ONE  /1.0D0/, SIX /6.0D0/, TWO /2.0D0/, ZERO /0.0D0/
C
      FLOTL = REAL(L)
      NQM1 = NQ - 1
      NQM2 = NQ - 2
      GO TO (100, 200), METH
C
C Set coefficients for Adams methods. ----------------------------------
 100  IF (NQ .NE. 1) GO TO 110
      EL(1) = ONE
      EL(2) = ONE
      TQ(1) = ONE
      TQ(2) = TWO
      TQ(3) = SIX*TQ(2)
      TQ(5) = ONE
      GO TO 300
 110  HSUM = H
      EM(1) = ONE
      FLOTNQ = FLOTL - ONE
      DO 115 I = 2, L
 115    EM(I) = ZERO
      DO 150 J = 1, NQM1
        IF ((J .NE. NQM1) .OR. (NQWAIT .NE. 1)) GO TO 130
        S = ONE
        CSUM = ZERO
        DO 120 I = 1, NQM1
          CSUM = CSUM + S*EM(I)/REAL(I+1)
 120      S = -S
        TQ(1) = EM(NQM1)/(FLOTNQ*CSUM)
 130    RXI = H/HSUM
        DO 140 IBACK = 1, J
          I = (J + 2) - IBACK
 140      EM(I) = EM(I) + EM(I-1)*RXI
        HSUM = HSUM + TAU(J)
 150    CONTINUE
C Compute integral from -1 to 0 of polynomial and of x times it. -------
      S = ONE
      EM0 = ZERO
      CSUM = ZERO
      DO 160 I = 1, NQ
        FLOTI = REAL(I)
        EM0 = EM0 + S*EM(I)/FLOTI
        CSUM = CSUM + S*EM(I)/(FLOTI+ONE)
 160    S = -S
C In EL, form coefficients of normalized integrated polynomial. --------
      S = ONE/EM0
      EL(1) = ONE
      DO 170 I = 1, NQ
 170    EL(I+1) = S*EM(I)/REAL(I)
      XI = HSUM/H
      TQ(2) = XI*EM0/CSUM
      TQ(5) = XI/EL(L)
      IF (NQWAIT .NE. 1) GO TO 300
C For higher order control constant, multiply polynomial by 1+x/xi(q). -
      RXI = ONE/XI
      DO 180 IBACK = 1, NQ
        I = (L + 1) - IBACK
 180    EM(I) = EM(I) + EM(I-1)*RXI
C Compute integral of polynomial. --------------------------------------
      S = ONE
      CSUM = ZERO
      DO 190 I = 1, L
        CSUM = CSUM + S*EM(I)/REAL(I+1)
 190    S = -S
      TQ(3) = FLOTL*EM0/CSUM
      GO TO 300
C
C Set coefficients for BDF methods. ------------------------------------
 200  DO 210 I = 3, L
 210    EL(I) = ZERO
      EL(1) = ONE
      EL(2) = ONE
      ALPH0 = -ONE
      AHATN0 = -ONE
      HSUM = H
      RXI = ONE
      RXIS = ONE
      IF (NQ .EQ. 1) GO TO 240
      DO 230 J = 1, NQM2
C In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------
        HSUM = HSUM + TAU(J)
        RXI = H/HSUM
        JP1 = J + 1
        ALPH0 = ALPH0 - ONE/REAL(JP1)
        DO 220 IBACK = 1, JP1
          I = (J + 3) - IBACK
 220      EL(I) = EL(I) + EL(I-1)*RXI
 230    CONTINUE
      ALPH0 = ALPH0 - ONE/REAL(NQ)
      RXIS = -EL(2) - ALPH0
      HSUM = HSUM + TAU(NQM1)
      RXI = H/HSUM
      AHATN0 = -EL(2) - RXI
      DO 235 IBACK = 1, NQ
        I = (NQ + 2) - IBACK
 235    EL(I) = EL(I) + EL(I-1)*RXIS
 240  T1 = ONE - AHATN0 + ALPH0
      T2 = ONE + REAL(NQ)*T1
      TQ(2) = ABS(ALPH0*T2/T1)
      TQ(5) = ABS(T2/(EL(L)*RXI/RXIS))
      IF (NQWAIT .NE. 1) GO TO 300
      CNQM1 = RXIS/EL(L)
      T3 = ALPH0 + ONE/REAL(NQ)
      T4 = AHATN0 + RXI
      ELP = T3/(ONE - T4 + T3)
      TQ(1) = ABS(ELP/CNQM1)
      HSUM = HSUM + TAU(NQ)
      RXI = H/HSUM
      T5 = ALPH0 - ONE/REAL(NQ+1)
      T6 = AHATN0 - RXI
      ELP = T2/(ONE - T6 + T5)
      TQ(3) = ABS(ELP*RXI*(FLOTL + ONE)*T5)
 300  TQ(4) = CORTES*TQ(2)
      RETURN
C----------------------- End of Subroutine DVSET -----------------------
      END
*DECK DVJUST
      SUBROUTINE DVJUST (YH, LDYH, IORD)
      DOUBLE PRECISION YH
      INTEGER LDYH, IORD
      DIMENSION YH(LDYH,*)
C-----------------------------------------------------------------------
C Call sequence input -- YH, LDYH, IORD
C Call sequence output -- YH
C COMMON block input -- NQ, METH, LMAX, HSCAL, TAU(13), N
C COMMON block variables accessed:
C     /DVOD01/ -- HSCAL, TAU(13), LMAX, METH, N, NQ,
C
C Subroutines called by DVJUST: DAXPY
C Function routines called by DVJUST: None
C-----------------------------------------------------------------------
C This subroutine adjusts the YH array on reduction of order,
C and also when the order is increased for the stiff option (METH = 2).
C Communication with DVJUST uses the following:
C IORD  = An integer flag used when METH = 2 to indicate an order
C         increase (IORD = +1) or an order decrease (IORD = -1).
C HSCAL = Step size H used in scaling of Nordsieck array YH.
C         (If IORD = +1, DVJUST assumes that HSCAL = TAU(1).)
C See References 1 and 2 for details.
C-----------------------------------------------------------------------
C
C Type declarations for labeled COMMON block DVOD01 --------------------
C
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for local variables --------------------------------
C
      DOUBLE PRECISION ALPH0, ALPH1, HSUM, ONE, PROD, T1, XI,XIOLD, ZERO
      INTEGER I, IBACK, J, JP1, LP1, NQM1, NQM2, NQP1
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE ONE, ZERO
C
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
C
      DATA ONE /1.0D0/, ZERO /0.0D0/
C
      IF ((NQ .EQ. 2) .AND. (IORD .NE. 1)) RETURN
      NQM1 = NQ - 1
      NQM2 = NQ - 2
      GO TO (100, 200), METH
C-----------------------------------------------------------------------
C Nonstiff option...
C Check to see if the order is being increased or decreased.
C-----------------------------------------------------------------------
 100  CONTINUE
      IF (IORD .EQ. 1) GO TO 180
C Order decrease. ------------------------------------------------------
      DO 110 J = 1, LMAX
 110    EL(J) = ZERO
      EL(2) = ONE
      HSUM = ZERO
      DO 130 J = 1, NQM2
C Construct coefficients of x*(x+xi(1))*...*(x+xi(j)). -----------------
        HSUM = HSUM + TAU(J)
        XI = HSUM/HSCAL
        JP1 = J + 1
        DO 120 IBACK = 1, JP1
          I = (J + 3) - IBACK
 120      EL(I) = EL(I)*XI + EL(I-1)
 130    CONTINUE
C Construct coefficients of integrated polynomial. ---------------------
      DO 140 J = 2, NQM1
 140    EL(J+1) = REAL(NQ)*EL(J)/REAL(J)
C Subtract correction terms from YH array. -----------------------------
      DO 170 J = 3, NQ
        DO 160 I = 1, N
 160      YH(I,J) = YH(I,J) - YH(I,L)*EL(J)
 170    CONTINUE
      RETURN
C Order increase. ------------------------------------------------------
C Zero out next column in YH array. ------------------------------------
 180  CONTINUE
      LP1 = L + 1
      DO 190 I = 1, N
 190    YH(I,LP1) = ZERO
      RETURN
C-----------------------------------------------------------------------
C Stiff option...
C Check to see if the order is being increased or decreased.
C-----------------------------------------------------------------------
 200  CONTINUE
      IF (IORD .EQ. 1) GO TO 300
C Order decrease. ------------------------------------------------------
      DO 210 J = 1, LMAX
 210    EL(J) = ZERO
      EL(3) = ONE
      HSUM = ZERO
      DO 230 J = 1,NQM2
C Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
        HSUM = HSUM + TAU(J)
        XI = HSUM/HSCAL
        JP1 = J + 1
        DO 220 IBACK = 1, JP1
          I = (J + 4) - IBACK
 220      EL(I) = EL(I)*XI + EL(I-1)
 230    CONTINUE
C Subtract correction terms from YH array. -----------------------------
      DO 250 J = 3,NQ
        DO 240 I = 1, N
 240      YH(I,J) = YH(I,J) - YH(I,L)*EL(J)
 250    CONTINUE
      RETURN
C Order increase. ------------------------------------------------------
 300  DO 310 J = 1, LMAX
 310    EL(J) = ZERO
      EL(3) = ONE
      ALPH0 = -ONE
      ALPH1 = ONE
      PROD = ONE
      XIOLD = ONE
      HSUM = HSCAL
      IF (NQ .EQ. 1) GO TO 340
      DO 330 J = 1, NQM1
C Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
        JP1 = J + 1
        HSUM = HSUM + TAU(JP1)
        XI = HSUM/HSCAL
        PROD = PROD*XI
        ALPH0 = ALPH0 - ONE/REAL(JP1)
        ALPH1 = ALPH1 + ONE/XI
        DO 320 IBACK = 1, JP1
          I = (J + 4) - IBACK
 320      EL(I) = EL(I)*XIOLD + EL(I-1)
        XIOLD = XI
 330    CONTINUE
 340  CONTINUE
      T1 = (-ALPH0 - ALPH1)/PROD
C Load column L + 1 in YH array. ---------------------------------------
      LP1 = L + 1
      DO 350 I = 1, N
 350    YH(I,LP1) = T1*YH(I,LMAX)
C Add correction terms to YH array. ------------------------------------
      NQP1 = NQ + 1
      DO 370 J = 3, NQP1
        CALL DAXPY (N, EL(J), YH(1,LP1), 1, YH(1,J), 1 )
 370  CONTINUE
      RETURN
C----------------------- End of Subroutine DVJUST ----------------------
      END
*DECK DVNLSD
      SUBROUTINE DVNLSD (Y, YH, LDYH, VSAV, SAVF, EWT, ACOR, IWM, WM,
     1                 F, JAC, PDUM, NFLAG, RPAR, IPAR)
      EXTERNAL F, JAC, PDUM
      DOUBLE PRECISION Y, YH, VSAV, SAVF, EWT, ACOR, WM, RPAR
      INTEGER LDYH, IWM, NFLAG, IPAR
      DIMENSION Y(*), YH(LDYH,*), VSAV(*), SAVF(*), EWT(*), ACOR(*),
     1          IWM(*), WM(*), RPAR(*), IPAR(*)
C-----------------------------------------------------------------------
C Call sequence input -- Y, YH, LDYH, SAVF, EWT, ACOR, IWM, WM,
C                        F, JAC, NFLAG, RPAR, IPAR
C Call sequence output -- YH, ACOR, WM, IWM, NFLAG
C COMMON block variables accessed:
C     /DVOD01/ ACNRM, CRATE, DRC, H, RC, RL1, TQ(5), TN, ICF,
C                JCUR, METH, MITER, N, NSLP
C     /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C
C Subroutines called by DVNLSD: F, DAXPY, DCOPY, DSCAL, DVJAC, DVSOL
C Function routines called by DVNLSD: DVNORM
C-----------------------------------------------------------------------
C Subroutine DVNLSD is a nonlinear system solver, which uses functional
C iteration or a chord (modified Newton) method.  For the chord method
C direct linear algebraic system solvers are used.  Subroutine DVNLSD
C then handles the corrector phase of this integration package.
C
C Communication with DVNLSD is done with the following variables. (For
C more details, please see the comments in the driver subroutine.)
C
C Y          = The dependent variable, a vector of length N, input.
C YH         = The Nordsieck (Taylor) array, LDYH by LMAX, input
C              and output.  On input, it contains predicted values.
C LDYH       = A constant .ge. N, the first dimension of YH, input.
C VSAV       = Unused work array.
C SAVF       = A work array of length N.
C EWT        = An error weight vector of length N, input.
C ACOR       = A work array of length N, used for the accumulated
C              corrections to the predicted y vector.
C WM,IWM     = Real and integer work arrays associated with matrix
C              operations in chord iteration (MITER .ne. 0).
C F          = Dummy name for user supplied routine for f.
C JAC        = Dummy name for user supplied Jacobian routine.
C PDUM       = Unused dummy subroutine name.  Included for uniformity
C              over collection of integrators.
C NFLAG      = Input/output flag, with values and meanings as follows:
C              INPUT
C                  0 first call for this time step.
C                 -1 convergence failure in previous call to DVNLSD.
C                 -2 error test failure in DVSTEP.
C              OUTPUT
C                  0 successful completion of nonlinear solver.
C                 -1 convergence failure or singular matrix.
C                 -2 unrecoverable error in matrix preprocessing
C                    (cannot occur here).
C                 -3 unrecoverable error in solution (cannot occur
C                    here).
C RPAR, IPAR = Dummy names for user's real and integer work arrays.
C
C IPUP       = Own variable flag with values and meanings as follows:
C              0,            do not update the Newton matrix.
C              MITER .ne. 0, update Newton matrix, because it is the
C                            initial step, order was changed, the error
C                            test failed, or an update is indicated by
C                            the scalar RC or step counter NST.
C
C For more details, see comments in driver subroutine.
C-----------------------------------------------------------------------
C Type declarations for labeled COMMON block DVOD01 --------------------
C
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for labeled COMMON block DVOD02 --------------------
C
      DOUBLE PRECISION HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C
C Type declarations for local variables --------------------------------
C
      DOUBLE PRECISION CCMAX, CRDOWN, CSCALE, DCON, DEL, DELP, ONE,
     1     RDIV, TWO, ZERO
      INTEGER I, IERPJ, IERSL, M, MAXCOR, MSBP
C
C Type declaration for function subroutines called ---------------------
C
      DOUBLE PRECISION DVNORM
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE CCMAX, CRDOWN, MAXCOR, MSBP, RDIV, ONE, TWO, ZERO
C
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C
      DATA CCMAX /0.3D0/, CRDOWN /0.3D0/, MAXCOR /3/, MSBP /20/,
     1     RDIV  /2.0D0/
      DATA ONE /1.0D0/, TWO /2.0D0/, ZERO /0.0D0/
C-----------------------------------------------------------------------
C On the first step, on a change of method order, or after a
C nonlinear convergence failure with NFLAG = -2, set IPUP = MITER
C to force a Jacobian update when MITER .ne. 0.
C-----------------------------------------------------------------------
      IF (JSTART .EQ. 0) NSLP = 0
      IF (NFLAG .EQ. 0) ICF = 0
      IF (NFLAG .EQ. -2) IPUP = MITER
      IF ( (JSTART .EQ. 0) .OR. (JSTART .EQ. -1) ) IPUP = MITER
C If this is functional iteration, set CRATE .eq. 1 and drop to 220
      IF (MITER .EQ. 0) THEN
        CRATE = ONE
        GO TO 220
      ENDIF
C-----------------------------------------------------------------------
C RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
C When RC differs from 1 by more than CCMAX, IPUP is set to MITER
C to force DVJAC to be called, if a Jacobian is involved.
C In any case, DVJAC is called at least every MSBP steps.
C-----------------------------------------------------------------------
      DRC = ABS(RC-ONE)
      IF (DRC .GT. CCMAX .OR. NST .GE. NSLP+MSBP) IPUP = MITER
C-----------------------------------------------------------------------
C Up to MAXCOR corrector iterations are taken.  A convergence test is
C made on the r.m.s. norm of each correction, weighted by the error
C weight vector EWT.  The sum of the corrections is accumulated in the
C vector ACOR(i).  The YH array is not altered in the corrector loop.
C-----------------------------------------------------------------------
 220  M = 0
      DELP = ZERO
      CALL DCOPY (N, YH(1,1), 1, Y, 1 )
      CALL F (N, TN, Y, SAVF, RPAR, IPAR)
      NFE = NFE + 1
      IF (IPUP .LE. 0) GO TO 250
C-----------------------------------------------------------------------
C If indicated, the matrix P = I - h*rl1*J is reevaluated and
C preprocessed before starting the corrector iteration.  IPUP is set
C to 0 as an indicator that this has been done.
C-----------------------------------------------------------------------
      CALL DVJAC (Y, YH, LDYH, EWT, ACOR, SAVF, WM, IWM, F, JAC, IERPJ,
     1           RPAR, IPAR)
      IPUP = 0
      RC = ONE
      DRC = ZERO
      CRATE = ONE
      NSLP = NST
C If matrix is singular, take error return to force cut in step size. --
      IF (IERPJ .NE. 0) GO TO 430
 250  DO 260 I = 1,N
 260    ACOR(I) = ZERO
C This is a looping point for the corrector iteration. -----------------
 270  IF (MITER .NE. 0) GO TO 350
C-----------------------------------------------------------------------
C In the case of functional iteration, update Y directly from
C the result of the last function evaluation.
C-----------------------------------------------------------------------
      DO 280 I = 1,N
 280    SAVF(I) = RL1*(H*SAVF(I) - YH(I,2))
      DO 290 I = 1,N
 290    Y(I) = SAVF(I) - ACOR(I)
c      write(*,*) 'y(13) 2: ',y(13)
      DEL = DVNORM (N, Y, EWT)
      DO 300 I = 1,N
 300    Y(I) = YH(I,1) + SAVF(I)
c      write(*,*) 'y(13) 3: ',y(13)
      CALL DCOPY (N, SAVF, 1, ACOR, 1)
      GO TO 400
C-----------------------------------------------------------------------
C In the case of the chord method, compute the corrector error,
C and solve the linear system with that as right-hand side and
C P as coefficient matrix.  The correction is scaled by the factor
C 2/(1+RC) to account for changes in h*rl1 since the last DVJAC call.
C-----------------------------------------------------------------------
 350  DO 360 I = 1,N
 360    Y(I) = (RL1*H)*SAVF(I) - (RL1*YH(I,2) + ACOR(I))
c      write(*,*) 'y(13) 4: ',y(13)
      CALL DVSOL (WM, IWM, Y, IERSL)
      NNI = NNI + 1
      IF (IERSL .GT. 0) GO TO 410
      IF (METH .EQ. 2 .AND. RC .NE. ONE) THEN
        CSCALE = TWO/(ONE + RC)
        CALL DSCAL (N, CSCALE, Y, 1)
      ENDIF
      DEL = DVNORM (N, Y, EWT)
      CALL DAXPY (N, ONE, Y, 1, ACOR, 1)
      DO 380 I = 1,N
 380    Y(I) = YH(I,1) + ACOR(I)
c      write(*,*) 'y(13) 5: ',y(13)
C-----------------------------------------------------------------------
C Test for convergence.  If M .gt. 0, an estimate of the convergence
C rate constant is stored in CRATE, and this is used in the test.
C-----------------------------------------------------------------------
 400  IF (M .NE. 0) CRATE = MAX(CRDOWN*CRATE,DEL/DELP)
      DCON = DEL*MIN(ONE,CRATE)/TQ(4)
      IF (DCON .LE. ONE) GO TO 450
      M = M + 1
      IF (M .EQ. MAXCOR) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. RDIV*DELP) GO TO 410
      DELP = DEL
      CALL F (N, TN, Y, SAVF, RPAR, IPAR)
      NFE = NFE + 1
      GO TO 270
C
 410  IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GO TO 430
      ICF = 1
      IPUP = MITER
      GO TO 220
C
 430  CONTINUE
      NFLAG = -1
      ICF = 2
      IPUP = MITER
      RETURN
C
C Return for successful step. ------------------------------------------
 450  NFLAG = 0
      JCUR = 0
      ICF = 0
      IF (M .EQ. 0) ACNRM = DEL
      IF (M .GT. 0) ACNRM = DVNORM (N, ACOR, EWT)
      RETURN
C----------------------- End of Subroutine DVNLSD ----------------------
      END
*DECK DVJAC
      SUBROUTINE DVJAC (Y, YH, LDYH, EWT, FTEM, SAVF, WM, IWM, F, JAC,
     1                 IERPJ, RPAR, IPAR)
      EXTERNAL F, JAC
      DOUBLE PRECISION Y, YH, EWT, FTEM, SAVF, WM, RPAR
      INTEGER LDYH, IWM, IERPJ, IPAR
      DIMENSION Y(*), YH(LDYH,*), EWT(*), FTEM(*), SAVF(*),
     1   WM(*), IWM(*), RPAR(*), IPAR(*)
C-----------------------------------------------------------------------
C Call sequence input -- Y, YH, LDYH, EWT, FTEM, SAVF, WM, IWM,
C                        F, JAC, RPAR, IPAR
C Call sequence output -- WM, IWM, IERPJ
C COMMON block variables accessed:
C     /DVOD01/  CCMXJ, DRC, H, RL1, TN, UROUND, ICF, JCUR, LOCJS,
C               MITER, MSBJ, N, NSLJ
C     /DVOD02/  NFE, NST, NJE, NLU
C
C Subroutines called by DVJAC: F, JAC, DACOPY, DCOPY, DGBFA, DGEFA,
C                              DSCAL
C Function routines called by DVJAC: DVNORM
C-----------------------------------------------------------------------
C DVJAC is called by DVNLSD to compute and process the matrix
C P = I - h*rl1*J , where J is an approximation to the Jacobian.
C Here J is computed by the user-supplied routine JAC if
C MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
C If MITER = 3, a diagonal approximation to J is used.
C If JSV = -1, J is computed from scratch in all cases.
C If JSV = 1 and MITER = 1, 2, 4, or 5, and if the saved value of J is
C considered acceptable, then P is constructed from the saved J.
C J is stored in wm and replaced by P.  If MITER .ne. 3, P is then
C subjected to LU decomposition in preparation for later solution
C of linear systems with P as coefficient matrix. This is done
C by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
C
C Communication with DVJAC is done with the following variables.  (For
C more details, please see the comments in the driver subroutine.)
C Y          = Vector containing predicted values on entry.
C YH         = The Nordsieck array, an LDYH by LMAX array, input.
C LDYH       = A constant .ge. N, the first dimension of YH, input.
C EWT        = An error weight vector of length N.
C SAVF       = Array containing f evaluated at predicted y, input.
C WM         = Real work space for matrices.  In the output, it containS
C              the inverse diagonal matrix if MITER = 3 and the LU
C              decomposition of P if MITER is 1, 2 , 4, or 5.
C              Storage of matrix elements starts at WM(3).
C              Storage of the saved Jacobian starts at WM(LOCJS).
C              WM also contains the following matrix-related data:
C              WM(1) = SQRT(UROUND), used in numerical Jacobian step.
C              WM(2) = H*RL1, saved for later use if MITER = 3.
C IWM        = Integer work space containing pivot information,
C              starting at IWM(31), if MITER is 1, 2, 4, or 5.
C              IWM also contains band parameters ML = IWM(1) and
C              MU = IWM(2) if MITER is 4 or 5.
C F          = Dummy name for the user supplied subroutine for f.
C JAC        = Dummy name for the user supplied Jacobian subroutine.
C RPAR, IPAR = Dummy names for user's real and integer work arrays.
C RL1        = 1/EL(2) (input).
C IERPJ      = Output error flag,  = 0 if no trouble, 1 if the P
C              matrix is found to be singular.
C JCUR       = Output flag to indicate whether the Jacobian matrix
C              (or approximation) is now current.
C              JCUR = 0 means J is not current.
C              JCUR = 1 means J is current.
C-----------------------------------------------------------------------
C
C Type declarations for labeled COMMON block DVOD01 --------------------
C
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for labeled COMMON block DVOD02 --------------------
C
      DOUBLE PRECISION HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C
C Type declarations for local variables --------------------------------
C
      DOUBLE PRECISION CON, DI, FAC, HRL1, ONE, PT1, R, R0, SRUR, THOU,
     1     YI, YJ, YJJ, ZERO
      INTEGER I, I1, I2, IER, II, J, J1, JJ, JOK, LENP, MBA, MBAND,
     1        MEB1, MEBAND, ML, ML3, MU, NP1
C
C Type declaration for function subroutines called ---------------------
C
      DOUBLE PRECISION DVNORM
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this subroutine.
C-----------------------------------------------------------------------
      SAVE ONE, PT1, THOU, ZERO
C-----------------------------------------------------------------------
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C
      DATA ONE /1.0D0/, THOU /1000.0D0/, ZERO /0.0D0/, PT1 /0.1D0/
C
      IERPJ = 0
      HRL1 = H*RL1
C See whether J should be evaluated (JOK = -1) or not (JOK = 1). -------
      JOK = JSV
      IF (JSV .EQ. 1) THEN
        IF (NST .EQ. 0 .OR. NST .GT. NSLJ+MSBJ) JOK = -1
        IF (ICF .EQ. 1 .AND. DRC .LT. CCMXJ) JOK = -1
        IF (ICF .EQ. 2) JOK = -1
      ENDIF
C End of setting JOK. --------------------------------------------------
C
      IF (JOK .EQ. -1 .AND. MITER .EQ. 1) THEN
C If JOK = -1 and MITER = 1, call JAC to evaluate Jacobian. ------------
      NJE = NJE + 1
      NSLJ = NST
      JCUR = 1
      LENP = N*N
      DO 110 I = 1,LENP
 110    WM(I+2) = ZERO
      CALL JAC (N, TN, Y, 0, 0, WM(3), N, RPAR, IPAR)
      IF (JSV .EQ. 1) CALL DCOPY (LENP, WM(3), 1, WM(LOCJS), 1)
      ENDIF
C
      IF (JOK .EQ. -1 .AND. MITER .EQ. 2) THEN
C If MITER = 2, make N calls to F to approximate the Jacobian. ---------
      NJE = NJE + 1
      NSLJ = NST
      JCUR = 1
      FAC = DVNORM (N, SAVF, EWT)
      R0 = THOU*ABS(H)*UROUND*REAL(N)*FAC
      IF (R0 .EQ. ZERO) R0 = ONE
      SRUR = WM(1)
      J1 = 2
      DO 230 J = 1,N
        YJ = Y(J)
        R = MAX(SRUR*ABS(YJ),R0/EWT(J))
        Y(J) = Y(J) + R
        FAC = ONE/R
        CALL F (N, TN, Y, FTEM, RPAR, IPAR)
        DO 220 I = 1,N
 220      WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
        Y(J) = YJ
        J1 = J1 + N
 230    CONTINUE
c      write(*,*) 'y(13) 6: ',y(13)
      NFE = NFE + N
      LENP = N*N
      IF (JSV .EQ. 1) CALL DCOPY (LENP, WM(3), 1, WM(LOCJS), 1)
      ENDIF
C
      IF (JOK .EQ. 1 .AND. (MITER .EQ. 1 .OR. MITER .EQ. 2)) THEN
      JCUR = 0
      LENP = N*N
      CALL DCOPY (LENP, WM(LOCJS), 1, WM(3), 1)
      ENDIF
C
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
C Multiply Jacobian by scalar, add identity, and do LU decomposition. --
      CON = -HRL1
      CALL DSCAL (LENP, CON, WM(3), 1)
      J = 3
      NP1 = N + 1
      DO 250 I = 1,N
        WM(J) = WM(J) + ONE
 250    J = J + NP1
      NLU = NLU + 1
      CALL DGEFA (WM(3), N, N, IWM(31), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
      ENDIF
C End of code block for MITER = 1 or 2. --------------------------------
C
      IF (MITER .EQ. 3) THEN
C If MITER = 3, construct a diagonal approximation to J and P. ---------
      NJE = NJE + 1
      JCUR = 1
      WM(2) = HRL1
      R = RL1*PT1
      DO 310 I = 1,N
 310    Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))
c      write(*,*) 'y(13) 7: ',y(13)
      CALL F (N, TN, Y, WM(3), RPAR, IPAR)
      NFE = NFE + 1
      DO 320 I = 1,N
        R0 = H*SAVF(I) - YH(I,2)
        DI = PT1*R0 - H*(WM(I+2) - SAVF(I))
        WM(I+2) = ONE
        IF (ABS(R0) .LT. UROUND/EWT(I)) GO TO 320
        IF (ABS(DI) .EQ. ZERO) GO TO 330
        WM(I+2) = PT1*R0/DI
 320    CONTINUE
      RETURN
 330  IERPJ = 1
      RETURN
      ENDIF
C End of code block for MITER = 3. -------------------------------------
C
C Set constants for MITER = 4 or 5. ------------------------------------
      ML = IWM(1)
      MU = IWM(2)
      ML3 = ML + 3
      MBAND = ML + MU + 1
      MEBAND = MBAND + ML
      LENP = MEBAND*N
C
      IF (JOK .EQ. -1 .AND. MITER .EQ. 4) THEN
C If JOK = -1 and MITER = 4, call JAC to evaluate Jacobian. ------------
      NJE = NJE + 1
      NSLJ = NST
      JCUR = 1
      DO 410 I = 1,LENP
 410    WM(I+2) = ZERO
      CALL JAC (N, TN, Y, ML, MU, WM(ML3), MEBAND, RPAR, IPAR)
      IF (JSV .EQ. 1)
     1   CALL DACOPY (MBAND, N, WM(ML3), MEBAND, WM(LOCJS), MBAND)
      ENDIF
C
      IF (JOK .EQ. -1 .AND. MITER .EQ. 5) THEN
C If MITER = 5, make ML+MU+1 calls to F to approximate the Jacobian. ---
      NJE = NJE + 1
      NSLJ = NST
      JCUR = 1
      MBA = MIN(MBAND,N)
      MEB1 = MEBAND - 1
      SRUR = WM(1)
      FAC = DVNORM (N, SAVF, EWT)
      R0 = THOU*ABS(H)*UROUND*REAL(N)*FAC
      IF (R0 .EQ. ZERO) R0 = ONE
      DO 560 J = 1,MBA
        DO 530 I = J,N,MBAND
          YI = Y(I)
          R = MAX(SRUR*ABS(YI),R0/EWT(I))
 530      Y(I) = Y(I) + R
        CALL F (N, TN, Y, FTEM, RPAR, IPAR)
        DO 550 JJ = J,N,MBAND
          Y(JJ) = YH(JJ,1)
          YJJ = Y(JJ)
          R = MAX(SRUR*ABS(YJJ),R0/EWT(JJ))
          FAC = ONE/R
          I1 = MAX(JJ-MU,1)
          I2 = MIN(JJ+ML,N)
          II = JJ*MEB1 - ML + 2
          DO 540 I = I1,I2
 540        WM(II+I) = (FTEM(I) - SAVF(I))*FAC
 550      CONTINUE
 560    CONTINUE
c      write(*,*) 'y(13) 8: ',y(13)
      NFE = NFE + MBA
      IF (JSV .EQ. 1)
     1   CALL DACOPY (MBAND, N, WM(ML3), MEBAND, WM(LOCJS), MBAND)
      ENDIF
C
      IF (JOK .EQ. 1) THEN
      JCUR = 0
      CALL DACOPY (MBAND, N, WM(LOCJS), MBAND, WM(ML3), MEBAND)
      ENDIF
C
C Multiply Jacobian by scalar, add identity, and do LU decomposition.
      CON = -HRL1
      CALL DSCAL (LENP, CON, WM(3), 1 )
      II = MBAND + 2
      DO 580 I = 1,N
        WM(II) = WM(II) + ONE
 580    II = II + MEBAND
      NLU = NLU + 1
      CALL DGBFA (WM(3), MEBAND, N, ML, MU, IWM(31), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C End of code block for MITER = 4 or 5. --------------------------------
C
C----------------------- End of Subroutine DVJAC -----------------------
      END
*DECK DACOPY
      SUBROUTINE DACOPY (NROW, NCOL, A, NROWA, B, NROWB)
      DOUBLE PRECISION A, B
      INTEGER NROW, NCOL, NROWA, NROWB
      DIMENSION A(NROWA,NCOL), B(NROWB,NCOL)
C-----------------------------------------------------------------------
C Call sequence input -- NROW, NCOL, A, NROWA, NROWB
C Call sequence output -- B
C COMMON block variables accessed -- None
C
C Subroutines called by DACOPY: DCOPY
C Function routines called by DACOPY: None
C-----------------------------------------------------------------------
C This routine copies one rectangular array, A, to another, B,
C where A and B may have different row dimensions, NROWA and NROWB.
C The data copied consists of NROW rows and NCOL columns.
C-----------------------------------------------------------------------
      INTEGER IC
C
      DO 20 IC = 1,NCOL
        CALL DCOPY (NROW, A(1,IC), 1, B(1,IC), 1)
 20     CONTINUE
C
      RETURN
C----------------------- End of Subroutine DACOPY ----------------------
      END
*DECK DVSOL
      SUBROUTINE DVSOL (WM, IWM, X, IERSL)
      DOUBLE PRECISION WM, X
      INTEGER IWM, IERSL
      DIMENSION WM(*), IWM(*), X(*)
C-----------------------------------------------------------------------
C Call sequence input -- WM, IWM, X
C Call sequence output -- X, IERSL
C COMMON block variables accessed:
C     /DVOD01/ -- H, RL1, MITER, N
C
C Subroutines called by DVSOL: DGESL, DGBSL
C Function routines called by DVSOL: None
C-----------------------------------------------------------------------
C This routine manages the solution of the linear system arising from
C a chord iteration.  It is called if MITER .ne. 0.
C If MITER is 1 or 2, it calls DGESL to accomplish this.
C If MITER = 3 it updates the coefficient H*RL1 in the diagonal
C matrix, and then computes the solution.
C If MITER is 4 or 5, it calls DGBSL.
C Communication with DVSOL uses the following variables:
C WM    = Real work space containing the inverse diagonal matrix if
C         MITER = 3 and the LU decomposition of the matrix otherwise.
C         Storage of matrix elements starts at WM(3).
C         WM also contains the following matrix-related data:
C         WM(1) = SQRT(UROUND) (not used here),
C         WM(2) = HRL1, the previous value of H*RL1, used if MITER = 3.
C IWM   = Integer work space containing pivot information, starting at
C         IWM(31), if MITER is 1, 2, 4, or 5.  IWM also contains band
C         parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C X     = The right-hand side vector on input, and the solution vector
C         on output, of length N.
C IERSL = Output flag.  IERSL = 0 if no trouble occurred.
C         IERSL = 1 if a singular matrix arose with MITER = 3.
C-----------------------------------------------------------------------
C
C Type declarations for labeled COMMON block DVOD01 --------------------
C
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for local variables --------------------------------
C
      INTEGER I, MEBAND, ML, MU
      DOUBLE PRECISION DI, HRL1, ONE, PHRL1, R, ZERO
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE ONE, ZERO
C
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
C
      DATA ONE /1.0D0/, ZERO /0.0D0/
C
      IERSL = 0
      GO TO (100, 100, 300, 400, 400), MITER
 100  CALL DGESL (WM(3), N, N, IWM(31), X, 0)
      RETURN
C
 300  PHRL1 = WM(2)
      HRL1 = H*RL1
      WM(2) = HRL1
      IF (HRL1 .EQ. PHRL1) GO TO 330
      R = HRL1/PHRL1
      DO 320 I = 1,N
        DI = ONE - R*(ONE - ONE/WM(I+2))
        IF (ABS(DI) .EQ. ZERO) GO TO 390
 320    WM(I+2) = ONE/DI
C
 330  DO 340 I = 1,N
 340    X(I) = WM(I+2)*X(I)
      RETURN
 390  IERSL = 1
      RETURN
C
 400  ML = IWM(1)
      MU = IWM(2)
      MEBAND = 2*ML + MU + 1
      CALL DGBSL (WM(3), MEBAND, N, ML, MU, IWM(31), X, 0)
      RETURN
C----------------------- End of Subroutine DVSOL -----------------------
      END
*DECK DVSRCO
      SUBROUTINE DVSRCO (RSAV, ISAV, JOB)
      DOUBLE PRECISION RSAV
      INTEGER ISAV, JOB
      DIMENSION RSAV(*), ISAV(*)
C-----------------------------------------------------------------------
C Call sequence input -- RSAV, ISAV, JOB
C Call sequence output -- RSAV, ISAV
C COMMON block variables accessed -- All of /DVOD01/ and /DVOD02/
C
C Subroutines/functions called by DVSRCO: None
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of the
C COMMON blocks DVOD01 and DVOD02, which are used internally by DVODE.
C
C RSAV = real array of length 49 or more.
C ISAV = integer array of length 41 or more.
C JOB  = flag indicating to save or restore the COMMON blocks:
C        JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV).
C        JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV).
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
      DOUBLE PRECISION RVOD1, RVOD2
      INTEGER IVOD1, IVOD2
      INTEGER I, LENIV1, LENIV2, LENRV1, LENRV2
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE LENRV1, LENIV1, LENRV2, LENIV2
C
      COMMON /DVOD01/ RVOD1(48), IVOD1(33)
      COMMON /DVOD02/ RVOD2(1), IVOD2(8)
      DATA LENRV1/48/, LENIV1/33/, LENRV2/1/, LENIV2/8/
C
      IF (JOB .EQ. 2) GO TO 100
      DO 10 I = 1,LENRV1
 10     RSAV(I) = RVOD1(I)
      DO 15 I = 1,LENRV2
 15     RSAV(LENRV1+I) = RVOD2(I)
C
      DO 20 I = 1,LENIV1
 20     ISAV(I) = IVOD1(I)
      DO 25 I = 1,LENIV2
 25     ISAV(LENIV1+I) = IVOD2(I)
C
      RETURN
C
 100  CONTINUE
      DO 110 I = 1,LENRV1
 110     RVOD1(I) = RSAV(I)
      DO 115 I = 1,LENRV2
 115     RVOD2(I) = RSAV(LENRV1+I)
C
      DO 120 I = 1,LENIV1
 120     IVOD1(I) = ISAV(I)
      DO 125 I = 1,LENIV2
 125     IVOD2(I) = ISAV(LENIV1+I)
C
      RETURN
C----------------------- End of Subroutine DVSRCO ----------------------
      END
*DECK DEWSET
      SUBROUTINE DEWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)
C***BEGIN PROLOGUE  DEWSET
C***SUBSIDIARY
C***PURPOSE  Set error weight vector.
C***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This subroutine sets the error weight vector EWT according to
C      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
C  with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
C  depending on the value of ITOL.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C***END PROLOGUE  DEWSET
C**End
      INTEGER N, ITOL
      INTEGER I
      DOUBLE PRECISION RTOL, ATOL, YCUR, EWT
      DIMENSION RTOL(*), ATOL(*), YCUR(N), EWT(N)
C
C***FIRST EXECUTABLE STATEMENT  DEWSET
      GO TO (10, 20, 30, 40), ITOL
 10   CONTINUE
      DO 15 I = 1,N
 15     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
      RETURN
 20   CONTINUE
      DO 25 I = 1,N
 25     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
      RETURN
 30   CONTINUE
      DO 35 I = 1,N
 35     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
      RETURN
 40   CONTINUE
      DO 45 I = 1,N
 45     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
      RETURN
C----------------------- END OF SUBROUTINE DEWSET ----------------------
      END
*DECK DVNORM
      DOUBLE PRECISION FUNCTION DVNORM (N, V, W)
C***BEGIN PROLOGUE  DVNORM
C***SUBSIDIARY
C***PURPOSE  Weighted root-mean-square vector norm.
C***TYPE      DOUBLE PRECISION (SVNORM-S, DVNORM-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This function routine computes the weighted root-mean-square norm
C  of the vector of length N contained in the array V, with weights
C  contained in the array W of length N:
C    DVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 )
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C***END PROLOGUE  DVNORM
C**End
      INTEGER N,   I
      DOUBLE PRECISION V, W,   SUM
      DIMENSION V(N), W(N)
C
C***FIRST EXECUTABLE STATEMENT  DVNORM
      SUM = 0.0D0
      DO 10 I = 1,N
 10     SUM = SUM + (V(I)*W(I))**2
      DVNORM = SQRT(SUM/N)
      RETURN
C----------------------- END OF FUNCTION DVNORM ------------------------
      END
*DECK XERRWD
      SUBROUTINE XERRWD (MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2)
C***BEGIN PROLOGUE  XERRWD
C***SUBSIDIARY
C***PURPOSE  Write error message with values.
C***CATEGORY  R3C
C***TYPE      DOUBLE PRECISION (XERRWV-S, XERRWD-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  Subroutines XERRWD, XSETF, XSETUN, and the function routine IXSAV,
C  as given here, constitute a simplified version of the SLATEC error
C  handling package.
C
C  All arguments are input arguments.
C
C  MSG    = The message (character array).
C  NMES   = The length of MSG (number of characters).
C  NERR   = The error number (not used).
C  LEVEL  = The error level..
C           0 or 1 means recoverable (control returns to caller).
C           2 means fatal (run is aborted--see note below).
C  NI     = Number of integers (0, 1, or 2) to be printed with message.
C  I1,I2  = Integers to be printed, depending on NI.
C  NR     = Number of reals (0, 1, or 2) to be printed with message.
C  R1,R2  = Reals to be printed, depending on NR.
C
C  Note..  this routine is machine-dependent and specialized for use
C  in limited context, in the following ways..
C  1. The argument MSG is assumed to be of type CHARACTER, and
C     the message is printed with a format of (1X,A).
C  2. The message is assumed to take only one line.
C     Multi-line messages are generated by repeated calls.
C  3. If LEVEL = 2, control passes to the statement   STOP
C     to abort the run.  This statement may be machine-dependent.
C  4. R1 and R2 are assumed to be in double precision and are printed
C     in D21.13 format.
C
C***ROUTINES CALLED  IXSAV
C***REVISION HISTORY  (YYMMDD)
C   920831  DATE WRITTEN
C   921118  Replaced MFLGSV/LUNSAV by IXSAV. (ACH)
C   930329  Modified prologue to SLATEC format. (FNF)
C   930407  Changed MSG from CHARACTER*1 array to variable. (FNF)
C   930922  Minor cosmetic change. (FNF)
C***END PROLOGUE  XERRWD
C
C*Internal Notes:
C
C For a different default logical unit number, IXSAV (or a subsidiary
C routine that it calls) will need to be modified.
C For a different run-abort command, change the statement following
C statement 100 at the end.
C-----------------------------------------------------------------------
C Subroutines called by XERRWD.. None
C Function routine called by XERRWD.. IXSAV
C-----------------------------------------------------------------------
C**End
C
C  Declare arguments.
C
      DOUBLE PRECISION R1, R2
      INTEGER NMES, NERR, LEVEL, NI, I1, I2, NR
      CHARACTER*(*) MSG
C
C  Declare local variables.
C
      INTEGER LUNIT, IXSAV, MESFLG
C
C  Get logical unit number and message print flag.
C
C***FIRST EXECUTABLE STATEMENT  XERRWD
      LUNIT = IXSAV (1, 0, .FALSE.)
      MESFLG = IXSAV (2, 0, .FALSE.)
      IF (MESFLG .EQ. 0) GO TO 100
C
C  Write the message.
C
      WRITE (LUNIT,10)  MSG
 10   FORMAT(1X,A)
      IF (NI .EQ. 1) WRITE (LUNIT, 20) I1
 20   FORMAT(6X,'In above message,  I1 =',I10)
      IF (NI .EQ. 2) WRITE (LUNIT, 30) I1,I2
 30   FORMAT(6X,'In above message,  I1 =',I10,3X,'I2 =',I10)
      IF (NR .EQ. 1) WRITE (LUNIT, 40) R1
 40   FORMAT(6X,'In above message,  R1 =',D21.13)
      IF (NR .EQ. 2) WRITE (LUNIT, 50) R1,R2
 50   FORMAT(6X,'In above,  R1 =',D21.13,3X,'R2 =',D21.13)
C
C  Abort the run if LEVEL = 2.
C
 100  IF (LEVEL .NE. 2) RETURN
      STOP
C----------------------- End of Subroutine XERRWD ----------------------
      END
*DECK XSETF
      SUBROUTINE XSETF (MFLAG)
C***BEGIN PROLOGUE  XSETF
C***PURPOSE  Reset the error print control flag.
C***CATEGORY  R3A
C***TYPE      ALL (XSETF-A)
C***KEYWORDS  ERROR CONTROL
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C   XSETF sets the error print control flag to MFLAG:
C      MFLAG=1 means print all messages (the default).
C      MFLAG=0 means no printing.
C
C***SEE ALSO  XERRWD, XERRWV
C***REFERENCES  (NONE)
C***ROUTINES CALLED  IXSAV
C***REVISION HISTORY  (YYMMDD)
C   921118  DATE WRITTEN
C   930329  Added SLATEC format prologue. (FNF)
C   930407  Corrected SEE ALSO section. (FNF)
C   930922  Made user-callable, and other cosmetic changes. (FNF)
C***END PROLOGUE  XSETF
C
C Subroutines called by XSETF.. None
C Function routine called by XSETF.. IXSAV
C-----------------------------------------------------------------------
C**End
      INTEGER MFLAG, JUNK, IXSAV
C
C***FIRST EXECUTABLE STATEMENT  XSETF
      IF (MFLAG .EQ. 0 .OR. MFLAG .EQ. 1) JUNK = IXSAV (2,MFLAG,.TRUE.)
      RETURN
C----------------------- End of Subroutine XSETF -----------------------
      END
*DECK XSETUN
      SUBROUTINE XSETUN (LUN)
C***BEGIN PROLOGUE  XSETUN
C***PURPOSE  Reset the logical unit number for error messages.
C***CATEGORY  R3B
C***TYPE      ALL (XSETUN-A)
C***KEYWORDS  ERROR CONTROL
C***DESCRIPTION
C
C   XSETUN sets the logical unit number for error messages to LUN.
C
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***SEE ALSO  XERRWD, XERRWV
C***REFERENCES  (NONE)
C***ROUTINES CALLED  IXSAV
C***REVISION HISTORY  (YYMMDD)
C   921118  DATE WRITTEN
C   930329  Added SLATEC format prologue. (FNF)
C   930407  Corrected SEE ALSO section. (FNF)
C   930922  Made user-callable, and other cosmetic changes. (FNF)
C***END PROLOGUE  XSETUN
C
C Subroutines called by XSETUN.. None
C Function routine called by XSETUN.. IXSAV
C-----------------------------------------------------------------------
C**End
      INTEGER LUN, JUNK, IXSAV
C
C***FIRST EXECUTABLE STATEMENT  XSETUN
      IF (LUN .GT. 0) JUNK = IXSAV (1,LUN,.TRUE.)
      RETURN
C----------------------- End of Subroutine XSETUN ----------------------
      END
*DECK IXSAV
      INTEGER FUNCTION IXSAV (IPAR, IVALUE, ISET)
C***BEGIN PROLOGUE  IXSAV
C***SUBSIDIARY
C***PURPOSE  Save and recall error message control parameters.
C***CATEGORY  R3C
C***TYPE      ALL (IXSAV-A)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  IXSAV saves and recalls one of two error message parameters:
C    LUNIT, the logical unit number to which messages are printed, and
C    MESFLG, the message print flag.
C  This is a modification of the SLATEC library routine J4SAVE.
C
C  Saved local variables..
C   LUNIT  = Logical unit number for messages.  The default is obtained
C            by a call to IUMACH (may be machine-dependent).
C   MESFLG = Print control flag..
C            1 means print all messages (the default).
C            0 means no printing.
C
C  On input..
C    IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
C    IVALUE = The value to be set for the parameter, if ISET = .TRUE.
C    ISET   = Logical flag to indicate whether to read or write.
C             If ISET = .TRUE., the parameter will be given
C             the value IVALUE.  If ISET = .FALSE., the parameter
C             will be unchanged, and IVALUE is a dummy argument.
C
C  On return..
C    IXSAV = The (old) value of the parameter.
C
C***SEE ALSO  XERRWD, XERRWV
C***ROUTINES CALLED  IUMACH
C***REVISION HISTORY  (YYMMDD)
C   921118  DATE WRITTEN
C   930329  Modified prologue to SLATEC format. (FNF)
C   930915  Added IUMACH call to get default output unit.  (ACH)
C   930922  Minor cosmetic changes. (FNF)
C   010425  Type declaration for IUMACH added. (ACH)
C***END PROLOGUE  IXSAV
C
C Subroutines called by IXSAV.. None
C Function routine called by IXSAV.. IUMACH
C-----------------------------------------------------------------------
C**End
      LOGICAL ISET
      INTEGER IPAR, IVALUE
C-----------------------------------------------------------------------
      INTEGER IUMACH, LUNIT, MESFLG
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this routine.
C-----------------------------------------------------------------------
      SAVE LUNIT, MESFLG
      DATA LUNIT/-1/, MESFLG/1/
C
C***FIRST EXECUTABLE STATEMENT  IXSAV
      IF (IPAR .EQ. 1) THEN
        IF (LUNIT .EQ. -1) LUNIT = IUMACH()
        IXSAV = LUNIT
        IF (ISET) LUNIT = IVALUE
        ENDIF
C
      IF (IPAR .EQ. 2) THEN
        IXSAV = MESFLG
        IF (ISET) MESFLG = IVALUE
        ENDIF
C
      RETURN
C----------------------- End of Function IXSAV -------------------------
      END
*DECK IUMACH
      INTEGER FUNCTION IUMACH()
C***BEGIN PROLOGUE  IUMACH
C***PURPOSE  Provide standard output unit number.
C***CATEGORY  R1
C***TYPE      INTEGER (IUMACH-I)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C *Usage:
C        INTEGER  LOUT, IUMACH
C        LOUT = IUMACH()
C
C *Function Return Values:
C     LOUT : the standard logical unit for Fortran output.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   930915  DATE WRITTEN
C   930922  Made user-callable, and other cosmetic changes. (FNF)
C***END PROLOGUE  IUMACH
C
C*Internal Notes:
C  The built-in value of 6 is standard on a wide range of Fortran
C  systems.  This may be machine-dependent.
C**End
C***FIRST EXECUTABLE STATEMENT  IUMACH
      IUMACH = 6
C
      RETURN
C----------------------- End of Function IUMACH ------------------------
      END
*DECK DUMACH
      DOUBLE PRECISION FUNCTION DUMACH ()
C***BEGIN PROLOGUE  DUMACH
C***PURPOSE  Compute the unit roundoff of the machine.
C***CATEGORY  R1
C***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C *Usage:
C        DOUBLE PRECISION  A, DUMACH
C        A = DUMACH()
C
C *Function Return Values:
C     A : the unit roundoff of the machine.
C
C *Description:
C     The unit roundoff is defined as the smallest positive machine
C     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
C     in a machine-independent manner.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DUMSUM
C***REVISION HISTORY  (YYYYMMDD)
C   19930216  DATE WRITTEN
C   19930818  Added SLATEC-format prologue.  (FNF)
C   20030707  Added DUMSUM to force normal storage of COMP.  (ACH)
C***END PROLOGUE  DUMACH
C
      DOUBLE PRECISION U, COMP
C***FIRST EXECUTABLE STATEMENT  DUMACH
      U = 1.0D0
 10   U = U*0.5D0
      CALL DUMSUM(1.0D0, U, COMP)
      IF (COMP .NE. 1.0D0) GO TO 10
      DUMACH = U*2.0D0
      RETURN
C----------------------- End of Function DUMACH ------------------------
      END
      SUBROUTINE DUMSUM(A,B,C)
C     Routine to force normal storing of A + B, for DUMACH.
      DOUBLE PRECISION A, B, C
      C = A + B
      RETURN
      END
C
C
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     constant times a vector plus a vector.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (DA.EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,4)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N.LT.4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
          DY(I) = DY(I) + DA*DX(I)
          DY(I+1) = DY(I+1) + DA*DX(I+1)
          DY(I+2) = DY(I+2) + DA*DX(I+2)
          DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
      END
C
C
      SUBROUTINE DSCAL(N,DA,DX,INCX)
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  Purpose
*  =======
**
*     scales a vector by a constant.
*     uses unrolled loops for increment equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
*
*        code for increment equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N.LT.5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          DX(I) = DA*DX(I)
          DX(I+1) = DA*DX(I+1)
          DX(I+2) = DA*DX(I+2)
          DX(I+3) = DA*DX(I+3)
          DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END
C
C
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     copies a vector, x, to a vector, y.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DY(IY) = DX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,7)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DY(I) = DX(I)
   30 CONTINUE
      IF (N.LT.7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
          DY(I) = DX(I)
          DY(I+1) = DX(I+1)
          DY(I+2) = DX(I+2)
          DY(I+3) = DX(I+3)
          DY(I+4) = DX(I+4)
          DY(I+5) = DX(I+5)
          DY(I+6) = DX(I+6)
   50 CONTINUE
      RETURN
      END
c
c     end dvode
c
c     begin linpack routines used by dvode
c
      SUBROUTINE DGEFA (A, LDA, N, IPVT, INFO)
C
C     This is an all-FORTRAN equivalent to the LINPACK routine of the
C     same name.  The LINPACK version should be used if it is available.
C
C   Part of DDRIV1,2,3 package, double precision version,
C   for use with a FORTRAN 77 compiler.
C   D. K. Kahaner, National Bureau of Standards, and
C   C. D. Sutherland, Los Alamos National Laboratory, January 2, 1983
C
      DOUBLE PRECISION A(LDA,N), T
      INTEGER IPVT(N)
C
      INFO = 0
      IF(N.EQ.1) THEN
        IF(A(1,1).EQ.0.D0) INFO = 1
        RETURN
      END IF
      DO 80 K = 1,N-1
        L = K
        DO 20 I = K+1,N
 20       IF(ABS(A(I,K)).GT.ABS(A(L,K))) L = I
        IPVT(K) = L
        IF(A(L,K).EQ.0.D0) THEN
          INFO = K
        ELSE
          IF(L.NE.K) THEN
            T = A(L,K)
            A(L,K) = A(K,K)
            A(K,K) = T
          END IF
          T = -1.D0/A(K,K)
          DO 40 I = K+1,N
 40         A(I,K) = T*A(I,K)
          DO 60 J = K+1,N
            T = A(L,J)
            IF (L.NE.K) THEN
              A(L,J) = A(K,J)
              A(K,J) = T
            END IF
            IF(T.NE.0.D0) THEN
              DO 50 I = K+1,N
 50             A(I,J) = A(I,J) + T*A(I,K)
            END IF
 60       CONTINUE
        END IF
 80     CONTINUE
      IPVT(N) = N
      IF (A(N,N).EQ.0.D0) INFO = N
      END
      SUBROUTINE DGESL (A, LDA, N, IPVT, B, JOB)
C
C     This is an all-FORTRAN equivalent to the LINPACK routine of the
C     same name.  The LINPACK version should be used if it is available.
C
C   Part of DDRIV1,2,3 package, double precision version,
C   for use with a FORTRAN 77 compiler.
C   D. K. Kahaner, National Bureau of Standards, and
C   C. D. Sutherland, Los Alamos National Laboratory, January 2, 1983
C
      DOUBLE PRECISION A(LDA,N), B(N), T
      INTEGER IPVT(N)
C
      IF(N.EQ.1) THEN
        B(1) = B(1)/A(1,1)
        RETURN
      END IF
      IF (JOB.EQ.0) THEN
        DO 40 K = 1,N-1
          L = IPVT(K)
          T = B(L)
          IF (L.NE.K) THEN
            B(L) = B(K)
            B(K) = T
          END IF
          IF (T.NE.0.D0) THEN
            DO 30 J = K+1,N
 30           B(J) = B(J) + T*A(J,K)
          END IF
 40       CONTINUE
        DO 60 K = N,1,-1
          B(K) = B(K)/A(K,K)
          IF (K.GT.1) THEN
            T = -B(K)
            IF(T.NE.0.D0) THEN
              DO 50 J = 1,K-1
 50             B(J) = B(J) + T*A(J,K)
            END IF
          END IF
 60       CONTINUE
      ELSE
C
C       JOB = NONZERO, SOLVE TRANS(A)*X = B
C       FIRST SOLVE TRANS(U)*Y = B
C
        DO 90 K = 1,N
          T = 0.D0
          DO 80 L = 1,K-1
 80         T = T + A(L,K)*B(L)
 90       B(K) = (B(K) - T)/A(K,K)
C
C       NOW SOLVE TRANS(L)*X = Y
C
        DO  110 K = N-1,1,-1
          T = 0.D0
          DO 100 L = K+1,N
 100        T = T + A(L,K)*B(L)
          B(K) = B(K) + T
          L = IPVT(K)
          IF (L.NE.K) THEN
            T = B(L)
            B(L) = B(K)
            B(K) = T
          END IF
 110      CONTINUE
      END IF
      END
      SUBROUTINE DGBFA (A, LDA, N, ML, MU, IPVT, INFO)
C
C     This is an all-FORTRAN equivalent to the LINPACK routine of the
C     same name.  The LINPACK version should be used if it is available.
C
C   Part of DDRIV1,2,3 package, double precision version,
C   for use with a FORTRAN 77 compiler.
C   D. K. Kahaner, National Bureau of Standards, and
C   C. D. Sutherland, Los Alamos National Laboratory, January 2, 1983
C
      DOUBLE PRECISION A(LDA,N), T
      INTEGER IPVT(N)
C
      INFO = 0
      IF (N.EQ.1) THEN
        IF (A(1,1).EQ.0.D0) INFO = 1
        RETURN
      END IF
      IF (ML.EQ.0) THEN
        DO 110 J = 1,N
          IPVT(J) = J
          IF (A(MU+1,J).EQ.0.D0) INFO = J
 110      CONTINUE
        RETURN
      END IF
      MW = ML + MU + 1
      J0 = MU + 2
      J1 = MIN(N, MW) - 1
      DO 130 JZ = J0,J1
        DO 120 I = MW+1-JZ,ML
 120      A(I,JZ) = 0.D0
 130    CONTINUE
      JZ = J1
      JU = 0
      DO 210 K = 1,N-1
        JZ = JZ + 1
        IF (JZ.LE.N) THEN
          DO 140 I = 1,ML
 140        A(I,JZ) = 0.D0
        END IF
        L = MW
        J1 = MW + 1
        J2 = MIN(MW+N-K, MW+ML)
        DO 150 J = J1,J2
          IF (ABS(A(J,K)).GT.ABS(A(L,K))) L = J
 150      CONTINUE
        IPVT(K) = L + K - MW
        IF (A(L,K).EQ.0.D0) THEN
          INFO = K
        ELSE
          IF (L.NE.MW) THEN
            T = A(L,K)
            A(L,K) = A(MW,K)
            A(MW,K) = T
          END IF
          T = -1.D0/A(MW,K)
          DO 160 J = J1,J2
 160        A(J,K) = T*A(J,K)
          JU = MIN(MAX(JU, MU+IPVT(K)), N)
          MM = MW
          DO 190 J = K+1,JU
            L = L - 1
            MM = MM - 1
            T = A(L,J)
            IF (L.NE.MM) THEN
              A(L,J) = A(MM,J)
              A(MM,J) = T
            END IF
            IF (T.NE.0.D0) THEN
              DO 180 I = J1,J2
                I3 = I + MM - MW
 180            A(I3,J) = A(I3,J) + T*A(I,K)
            END IF
 190        CONTINUE
        END IF
 210    CONTINUE
      IPVT(N) = N
      IF (A(MW,N).EQ.0.D0) INFO = N
      END
      SUBROUTINE DGBSL (A, LDA, N, ML, MU, IPVT, B, JOB)
C
C     This is an all-FORTRAN equivalent to the LINPACK routine of the
C     same name.  The LINPACK version should be used if it is available.
C
C   Part of DDRIV1,2,3 package, double precision version,
C   for use with a FORTRAN 77 compiler.
C   D. K. Kahaner, National Bureau of Standards, and
C   C. D. Sutherland, Los Alamos National Laboratory, January 2, 1983
C
      DOUBLE PRECISION A(LDA,N), B(N), T
      INTEGER IPVT(N)
C
      IF (N.EQ.1) THEN
        B(1) = B(1)/A(1,1)
        RETURN
      END IF
C
      MW = ML + MU + 1
      IF (JOB.EQ.0) THEN
        IF (ML.NE.0) THEN
          DO 40 K = 1,N-1
            L = IPVT(K)
            T = B(L)
            IF (L.NE.K) THEN
              B(L) = B(K)
              B(K) = T
            END IF
            IF (T.NE.0.D0) THEN
              J1 = MW + 1
              J2 = MIN(MW+N-K, MW+ML)
              DO 30 J = J1,J2
 30             B(K+J-MW) = B(K+J-MW) + T*A(J,K)
            END IF
 40         CONTINUE
        END IF
        DO 70 K = N,1,-1
          B(K) = B(K)/A(MW,K)
          IF (K.NE.1) THEN
            IF (MW.NE.1) THEN
              T = -B(K)
              IF (T.NE.0.D0) THEN
                J1 = MAX(1, MW+1-K)
                J2 = MW - 1
                DO 60 J = J1,J2
 60               B(K+J-MW) = B(K+J-MW) + T*A(J,K)
              END IF
            END IF
          END IF
 70       CONTINUE
      ELSE
C
C       JOB = NONZERO, SOLVE TRANS(A)*X = B
C       FIRST SOLVE TRANS(U)*Y = B
C
        DO 100 K = 1,N
          LM = MIN(K, MW) - 1
          LA = MW - LM
          LB = K - LM
          T = 0.D0
          DO 90 L = 1,LM
 90         T = T + A(L+LA-1,K)*B(L+LB-1)
 100      B(K) = (B(K) - T)/A(MW,K)
C
C       NOW SOLVE TRANS(L)*X = Y
C
        IF (ML.NE.0) THEN
          DO  120 KB = 1,N-1
            K = N - KB
            LM = MIN(KB, ML)
            T = 0.D0
            DO 110 L = 1,LM
 110          T = T + A(L+MW,K)*B(L+K)
            B(K) = B(K) + T
            L = IPVT(K)
            IF (L.NE.K) THEN
              T = B(L)
              B(L) = B(K)
              B(K) = T
            END IF
 120        CONTINUE
        END IF
      END IF
      END
c
c     end linpack routines used by dvode
c


c ----------------------------------------------------------------------
      subroutine cubic(a, b, c, d, x, nroot)
c ----------------------------------------------------------------------
c Adapted from http://terpconnect.umd.edu/~nsw/ench250/cubiceq2.for
c ----------------------------------------------------------------------
c Solve a cubic equation where a, b, c, and d are real.
c   a*x**3 + b*x**2 + c*x + d = 0
c
c Public Variables
c   a, b, c, d  ... coefficients (input)
c   x(i)        ... three (generally) complex solutions (output)
c   nroot       ... number of roots (output)
c Local Variables:
c   y1, y2, y3  ... three transformed solutions
c
c Formula used are given in Tuma, "Engineering Mathematics Handbook", p7
c   (McGraw Hill, 1978).
c   Step 0: If a is 0. use the quadratic formula to avoid dividing by 0.
c   Step 1: Calculate p and q
c           p = ( 3*c/a - (b/a)**2 ) / 3
c           q = ( 2*(b/a)**3 - 9*b*c/a/a + 27*d/a ) / 27
c   Step 2: Calculate discriminant D
c           D = (p/3)**3 + (q/2)**2
c   Step 3: Depending on the sign of D, we follow different strategy.
c           If D<0, three distinct real roots.
c           If D=0, three real roots of which at least two are equal.
c           If D>0, one real and two complex roots.
c   Step 3a: For D>0 and D=0,
c           Calculate u and v
c           u = cubic_root(-q/2 + sqrt(D))
c           v = cubic_root(-q/2 - sqrt(D))
c           Find the three transformed roots
c           y1 = u + v
c           y2 = -(u+v)/2 + i (u-v)*sqrt(3)/2
c           y3 = -(u+v)/2 - i (u-v)*sqrt(3)/2
c   Step 3b Alternately, for D<0, a trigonometric formulation is more convenient
c           y1 =  2 * sqrt(|p|/3) * cos(phi/3)
c           y2 = -2 * sqrt(|p|/3) * cos((phi+pi)/3)
c           y3 = -2 * sqrt(|p|/3) * cos((phi-pi)/3)
c           where phi = acos(-q/2/sqrt(|p|**3/27))
c                 pi  = 3.141592654...
c   Step 4  Finally, find the three roots
c           x = y - b/a/3
c ----------------------------------------------------------------------

c Declare variables
      implicit none
      integer nroot
      doubleprecision a,b,c,d,x(3),pi,p,q,DD,phi,temp1,temp2,y1,y2,y3
      doubleprecision u,v,y2r,y2i
      data pi/3.14159265358979324d0/

c Step 0: If a is 0 use the quadratic formula. -------------------------
      if(a.eq.0d0)then
        write(*,*) 'ERROR: a=0 not implemented in subroutine cubic.'
        stop 33001
c        call quad(b, c, d, x, nroot)
        return
      endif

c Cubic equation with 3 roots
      nroot = 3

c Step 1: Calculate p and q --------------------------------------------
      p  = c/a - b*b/(3d0*a*a)
      q  = (2d0*(b/a)**3 - 9d0*b*c/(a*a) + 27d0*d/a) / 27d0

c Step 2: Calculate DD (discriminant) ----------------------------------
      DD = p*p*p/27d0 + 0.25d0*q*q

c Step 3: Branch to different algorithms based on DD -------------------
      if(DD.lt.0d0)then
c       Step 3b:
c       3 real unequal roots -- use the trigonometric formulation
        phi = acos(-0.5d0*q/sqrt(abs(p*p*p)/27d0))
        temp1 = 2d0*sqrt(abs(p)/3d0)
        y1 =  temp1*cos(phi/3d0)
        y2 = -temp1*cos((phi+pi)/3d0)
        y3 = -temp1*cos((phi-pi)/3d0)
      else
c       Step 3a:
c       1 real root & 2 conjugate complex roots OR 3 real roots (some are equal)
        temp1 = -0.5d0*q + sqrt(DD)
        temp2 = -0.5d0*q - sqrt(DD)
        u = abs(temp1)**(1d0/3d0)
        v = abs(temp2)**(1d0/3d0)
        if(temp1.lt.0d0) u=-u
        if(temp2.lt.0d0) v=-v
        y1  = u + v
        y2r = -(u+v)*0.5d0
        y2i =  (u-v)*sqrt(3d0)*0.5d0
      endif

c Step 4: Final transformation -----------------------------------------
      temp1 = b/(a*3d0)
      y1 = y1-temp1
      y2 = y2-temp1
      y3 = y3-temp1
      y2r=y2r-temp1

c Assign answers -------------------------------------------------------
      if(DD.lt.0d0)then
        x(1) = y1
        x(2) = y2
        x(3) = y3
      elseif(DD.eq.0d0)then
        x(1) = y1
        x(2) = y2r
        x(3) = y2r
      else
c        write(*,*) 'ERROR: complex roots in subroutine cubic'
c        stop 33002
        x(1) = y1
        x(2) = 0d0
        x(3) = 0d0
c        x(1) = cmplx( y1,  0.)
c        x(2) = cmplx(y2r, y2i)
c        x(3) = cmplx(y2r,-y2i)
      endif
c
      return
      end
c
c
c
      SUBROUTINE cbcrt (a,z)
!-----------------------------------------------------------------------

!        CBCRT COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
!              A(1) + A(2)*Z + A(3)*Z**2 + A(4)*Z**3
!     AND STORES THE RESULTS IN Z. IT IS ASSUMED THAT A(4) IS NONZERO.

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS
!        NAVAL SURFACE WEAPONS CENTER
!        DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------
!     Converted to be compatible with ELF90 by
!        Alan Miller
!        amiller @ bigpond.net.au
!     WWW-page: http://users.bigpond.net.au/amiller
!     Latest revision - 27 February 1997
!-----------------------------------------------------------------------
!     Adapted for use in DISKEVOL by Ruud Visser on 8 June 2009.
!     We only need real non-negative solutions, so complex solutions
!     are set to -1.
!-----------------------------------------------------------------------
      IMPLICIT NONE
      doubleprecision a(*),z(*),aq(3),arg,c,cf,d,eps,p,p1,q,q1,r,ra,rb
      doubleprecision rq,rt,rt3,r1,s,sf,sq,sum,t,tol,t1,w,w1,w2,x,x1,x2
      doubleprecision x3,y,y1,y2,y3,cbrt
      parameter(rt3=1.7320508075689D0,eps=1d-12)
!-----------------------------------------------------------------------
      IF (a(1).eq.0d0) GO TO 100
      p = a(3)/(3.0D0*a(4))
      q = a(2)/a(4)
      r = a(1)/a(4)
      tol = 4.0D0*eps

      c = 0d0
      t = a(2) - p*a(3)
      IF (ABS(t).gt.tol*ABS(a(2))) c = t/a(4)

      t = 2.0D0*p*p - q
      IF (ABS(t).le.tol*ABS(q)) t = 0d0
      d = r + p*t
      IF (ABS(d).le.tol*ABS(r)) GO TO 110

!           SET  SQ = (A(4)/S)**2 * (C**3/27 + D**2/4)

      s = MAX(ABS(a(1)), ABS(a(2)), ABS(a(3)))
      p1 = a(3)/(3.0D0*s)
      q1 = a(2)/s
      r1 = a(1)/s

      t1 = q - 2.25D0*p*p
      IF (ABS(t1).le.tol*ABS(q)) t1 = 0d0
      w = 0.25D0*r1*r1
      w1 = 0.5D0*p1*r1*t
      w2 = q1*q1*t1/27.0D0
      IF (w1.lt.0d0) GO TO 10
      w = w + w1
      sq = w + w2
      GO TO 12
   10 IF (w2.lt.0d0) GO TO 11
      w = w + w2
      sq = w + w1
      GO TO 12
   11 sq = w + (w1 + w2)
   12 IF (ABS(sq).le.tol*w) sq = 0d0
      rq = ABS(s/a(4))*SQRT(ABS(sq))
      IF (sq.ge.0d0) GO TO 40

!                   ALL ROOTS ARE REAL

      arg = ATAN2(rq, -0.5D0*d)
      cf = COS(arg/3.0D0)
      sf = SIN(arg/3.0D0)
      rt = SQRT(-c/3.0D0)
      y1 = 2.0D0*rt*cf
      y2 = -rt*(cf + rt3*sf)
      y3 = -(d/y1)/y2

      x1 = y1 - p
      x2 = y2 - p
      x3 = y3 - p
      IF (ABS(x1).le.ABS(x2)) GO TO 20
      t = x1
      x1 = x2
      x2 = t
   20 IF (ABS(x2).le.ABS(x3)) GO TO 30
      t = x2
      x2 = x3
      x3 = t
      IF (ABS(x1).le.ABS(x2)) GO TO 30
      t = x1
      x1 = x2
      x2 = t

   30 w = x3
      IF (ABS(x2).lt.0.1D0*ABS(x3)) GO TO 70
      IF (ABS(x1).lt.0.1D0*ABS(x2)) x1 = - (r/x3)/x2
      z(1) = x1
      z(2) = x2
      z(3) = x3
      GO TO 200

!                  REAL AND COMPLEX ROOTS

   40 ra = cbrt(-0.5D0*d - SIGN(rq,d))
      rb = -c/(3.0D0*ra)
      t = ra + rb
      w = -p
      x = -p
      IF (ABS(t).le.tol*ABS(ra)) GO TO 41
      w = t - p
      x = -0.5D0*t - p
      IF (ABS(x).le.tol*ABS(p)) x = 0d0
   41 t = ABS(ra - rb)
      y = 0.5D0*rt3*t

      IF (t.le.tol*ABS(ra)) GO TO 60
      IF (ABS(x).lt.ABS(y)) GO TO 50
      s = ABS(x)
      t = y/x
      GO TO 51
   50 s = ABS(y)
      t = x/y
   51 IF (s.lt.0.1D0*ABS(w)) GO TO 70
      w1 = w/s
      sum = 1.0D0 + t*t
      IF (w1*w1.lt.0.01D0*sum) w = - ((r/sum)/s)/s
      z(1) = w
      z(2) = -1d0
      z(3) = -1d0
      GO TO 200

!               AT LEAST TWO ROOTS ARE EQUAL

   60 IF (ABS(x).lt.ABS(w)) GO TO 61
      IF (ABS(w).lt.0.1D0*ABS(x)) w = - (r/x)/x
      z(1) = w
      z(2) = x
      z(3) = z(2)
      GO TO 200
   61 IF (ABS(x).lt.0.1D0*ABS(w)) GO TO 70
      z(1) = x
      z(2) = z(1)
      z(3) = w
      GO TO 200

!     HERE W IS MUCH LARGER IN MAGNITUDE THAN THE OTHER ROOTS.
!     AS A RESULT, THE OTHER ROOTS MAY BE EXCEEDINGLY INACCURATE
!     BECAUSE OF ROUNDOFF ERROR.  TO DEAL WITH THIS, A QUADRATIC
!     IS FORMED WHOSE ROOTS ARE THE SAME AS THE SMALLER ROOTS OF
!     THE CUBIC.  THIS QUADRATIC IS THEN SOLVED.

!     THIS CODE WAS WRITTEN BY WILLIAM L. DAVIS (NSWC).

   70 aq(1) = a(1)
      aq(2) = a(2) + a(1)/w
      aq(3) = -a(4)*w
      CALL qdcrt(aq, z)
      z(3) = w
      GO TO 200

!-----------------------------------------------------------------------

!                  CASE WHEN A(1) = 0

  100 a(1) = a(2)
      a(2) = a(3)
      a(3) = a(4)
      CALL qdcrt(a, z)
      z(3) = 0d0
      GO TO 200

!                   CASE WHEN D = 0

  110 z(1) = -p
      w = SQRT(ABS(c))
      IF (c.lt.0d0) GO TO 120
      z(2) = -1d0
      z(3) = -1d0
      GO TO 200

  120 IF (p.ne.0d0) GO TO 130
      z(2) = w
      z(3) = -w
      GO TO 200

  130 x = -(p + SIGN(w,p))
      z(3) = x
      t = 3.0D0*a(1)/(a(3)*x)
      IF (ABS(p).gt.ABS(t)) GO TO 131
      z(2) = t
      GO TO 200
  131 z(2) = z(1)
      z(1) = t
      GO TO 200

!     Finally, make sure the real non-negative root is in z(1).

  200 if(z(1).lt.0d0) then
          if(z(2).lt.0d0) then
              if(z(3).ge.0d0) then
                  t = z(1)
                  z(1) = z(3)
                  z(3) = t
              endif
          else
              if(z(3).lt.0d0) then
                  t = z(1)
                  z(1) = z(2)
                  z(2) = t
              endif
          endif
      endif
      RETURN

      END
c
c
c
      SUBROUTINE qdcrt (a, z)
!-----------------------------------------------------------------------

!        QDCRT COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
!              A(1) + A(2)*Z + A(3)*Z**2
!     AND STORES THE RESULTS IN Z.  IT IS ASSUMED THAT A(3) IS NONZERO.
!-----------------------------------------------------------------------
!     Converted to be compatible with ELF90 by
!        Alan Miller
!        amiller @ bigpond.net.au
!     WWW-page: http://users.bigpond.net.au/amiller
!     Latest revision - 27 February 1997
!-----------------------------------------------------------------------
!     Adapted for use in DISKEVOL by Ruud Visser on 8 June 2009
!     We only need real non-negative solutions, so complex solutions
!     are set to -1.
!-----------------------------------------------------------------------

      IMPLICIT NONE
      doubleprecision a(*),z(*),d,eps,r,w,x,y
      parameter(eps=1d-12)
!-------------------
      IF (a(1).eq.0d0) GO TO 40
      d = a(2)*a(2) - 4.0D0*a(1)*a(3)
      IF (ABS(d).le.2.0D0*eps*a(2)*a(2)) GO TO 20
      r = SQRT(ABS(d))
      IF (d.lt.0d0) GO TO 30

!                 DISTINCT REAL ROOTS

      IF (a(2).ne.0d0) GO TO 10
      x = ABS(0.5D0*r/a(3))
      z(1) = x
      z(2) = -x
      RETURN

   10 w = -(a(2) + SIGN(r,a(2)))
      z(1) = 2.0D0*a(1)/w
      z(2) = 0.5D0*w/a(3)
      RETURN

!                  EQUAL REAL ROOTS

   20 z(1) = -0.5D0*a(2)/a(3)
      z(2) = z(1)
      RETURN

!                   COMPLEX ROOTS

   30 z(1) = -1d0
      z(2) = -1d0
      RETURN

!                 CASE WHEN A(1) = 0

   40 z(1) = 0d0
      z(2) = -a(2)/a(3)
      RETURN

      END
c
c
c
      FUNCTION cbrt(x)
!-----------------------------------------------------------------------
!                   CUBE ROOT OF A REAL NUMBER
!-----------------------------------------------------------------------
!     Adapted from Alan Miller for use in DISKEVOL by Ruud Visser on
!     8 June 2009
!-----------------------------------------------------------------------
      IMPLICIT NONE
      doubleprecision x,r,cbrt
      
      IF (x.lt.0d0) THEN
          r = LOG(-x) / 3d0
          cbrt = -EXP(r)
      ELSE IF (x.gt.0d0) THEN
          r = LOG(x) / 3d0
          cbrt = EXP(r)
      ELSE
          cbrt = 0d0
      END IF

      RETURN
      END
c
