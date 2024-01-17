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

c        index of centrifugal radius
          icen = 0
          do ir=1,nr
              if(r(ir).lt.rcentr) icen=ir
          enddo

c        index of centrifugal radius *2
          icen2 = 0
          do ir=1,nr
              if(r(ir).lt.2d0*rcentr) icen2=ir
          enddo

          psi(1) = 0.d0
          mdotcap = 0d0
          mdotdisk = 0d0
          mdotflow = 0d0

!           if(idistr.eq.0) then
! c
! c             (Default) Distribute the matter over the disk in the Ulrich way
! c
!               do ir=1,icen-1
!                   mmu0      = sqrt(1-r(ir)/rcentr)
!                   rhocmueq  = mdot*(r(ir)/rcentr)/
!      %                 ((4*pi*sqrt(GG*mcentr*r(ir)**3))*(2*mmu0**2))
!                   vz        = sqrt(GG*mcentr/r(ir))*mmu0
!                   psi(ir)   = 2*rhocmueq*vz
!                   subkep(ir) = sqrt(1d0-mmu0**2)
!               enddo
!           else if(idistr.eq.1) then
! c
! c             Distribute the matter over the disk in such a way that
! c             there is no friction of the infalling matter with the disk
! c             (see Hueso & Guillot 2005)
! c
!               do ir=1,icen-1
!                   psi(ir)   = (mdot/(8*pi*r(ir)*rcentr)) *
!      %                        ((r(ir)/rcentr)*
!      %                        (1.d0-sqrt(r(ir)/rcentr)))**(-0.5)
!                   subkep(ir) = 1.d0
!               enddo
!           else if(idistr.eq.2) then
! c
! c     Distribute the matter over the disk in the Cassen & Moosman (1981)
! c     way as described in Terebey, Shu & Cassen (1984).
! c
! c              do ir=1,icen-1
!               do ir=1,icen
!                   zeta = rcentr / r(ir)
!                   mmu0 = sqrt(1d0-1d0/zeta)
!                   legendre = 0.5d0 * (3d0*mmu0**2-1d0)
!                   vx = -sqrt(GG*mcentr/r(ir))
!                   vz = -mmu0 * vx
!                   rhocmueq = -mdot / (4d0 * pi * r(ir)**2 * vx *
!      $                                (1d0 + 2d0*zeta*legendre))
!                   psi(ir) = 2d0 * rhocmueq * vz
!                   subkep(ir) = sqrt(1d0-mmu0**2)
!               enddo
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
    !  I should use vertden for now
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
! c
! c     Check if material can indeed reach each point on the surface.
! c
!              do ir=nr,1,-1
!                  if(accden(ir).eq.0d0) goto 100
!                  if(ir.ge.irzero) then
!                      acc_surface(ir) = 0d0
!                      if(r(ir).lt.rcentr) then
!                          zeta = rcentr / r(ir)
!                          disk_mu0(ir) = sqrt(1d0-1d0/zeta)
!                          call get_vertden(r(ir),0d0,r(ir),0.5d0*pi,
!      $                                    time,disk_mu0(ir),0)
!                          accden(ir) = mden * mpeffinv
!                      else
!                          accden(ir) = 0d0
!                      endif
!                      goto 100
!                  endif
!   232            continue
!                  docosth0 = 1
!                  costh0 = -1d0
!                  jr = ir
!                  ztemp = acc_surface(ir)
!                  ttemp = time
!   230            continue
!                  if(jr.ge.irzero) goto 100
!                  rsph = sqrt(r(jr)**2+ztemp**2)
!                  theta = atan2(r(jr),ztemp)
!                  costh = ztemp / rsph
!                  sinth = r(jr) / rsph
!                  call get_vertden(r(jr),ztemp,rsph,theta,time,costh0,0)
!                  if(docosth0.eq.1) then
!                      costh0 = global_mu0
!                      docosth0 = 0
!                  endif
!                  vx = ux
!                  vz = uz
! c
! c     If the velocities are zero, there is no accretion, so skip ahead.
! c
!                  if(jr.eq.ir.and.vx.eq.0d0.and.vz.eq.0d0) goto 231
! c
! c     Now take small steps back along the infall trajectory to see if it
! c     intersects the disk anywhere. If so, this trajectory cannot reach
! c     the point where we started, so we have to raise the disk surface
! c     there. If the trajectory originated from inside the outflow
! c     region, set the accretion density to zero and skip ahead.
! c
!                  if(jr.lt.nr) then
!                      dt = abs((r(jr+1)-r(jr)) / vx)
!                      ztemp = ztemp - vz * dt
!                      ttemp = ttemp - dt
!                      if(ttemp.le.0d0) goto 100
!                      if(number_invalid(vx).ne.0.or.
!      $                                     number_invalid(vz).ne.0) then
!                          write(*,*) time/year,r(jr)/au,ztemp/au,
!      $                              vx,vz,dt
!                          write(*,*) ir,disk_surface(ir)/au,
!      $rsph/(cloud_asound*ttemp),(cloud_omega*ttemp)**2
!                          stop 229
!                      endif
!                      outangle = 0.5519d0*(ttemp/tcloudm)**0.26d0
!                      if(outangle.gt.pi/2) outangle = pi/2
!                      aout = 0.191d0 / (ttemp/tcloudm)**3
!                      if(aout.lt.(0.191d0)) aout = 0.191d0
!                      if(ztemp.gt.aout*au*(r(jr+1)/au)**1.5d0.and.
!      $                                       time.gt.0.2d0*tcloudm) then
!                          accdenf(ir) = accden(ir)
!                          accden(ir) = 0d0
!                          goto 100
!                      endif
!                      if(ztemp.lt.disk_surface(jr+1)) then
!                          accden(ir) = 0d0
!                          if(acc_surface(ir).gt.0d0) then
!                              acc_surface(ir) = acc_surface(ir) +
!      $                                        (disk_surface(jr+1)-ztemp)
!                              goto 232
!                          endif
!                          goto 100
!                      endif
! c
! c     Take another step if we haven't reached the edge of the grid.
! c
!                      if(jr.lt.nr-1) then
!                          jr = jr + 1
!                          goto 230
!                      endif
!                  endif
!                  goto 100
!   231            continue
!                  accden(ir) = 0d0
!   100            continue
!              enddo
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
!           if(idistr.ne.3) then
!              do ir=max(icen,1),nr
!                 psi(ir) = 0.d0
!              enddo
!           endif
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
