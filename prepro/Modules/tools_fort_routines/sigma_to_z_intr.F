#define CUBIC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigma_to_z_intr (Lm,Mm,N, nz, z_r, z_w, rmask, var,
     &                    z_lev, var_zlv, imin,jmin,kmin, FillValue)
!
! Interpolate field "var" defined in sigma-space to 3-D z_lev. 
!


      implicit none

      integer Lm,Mm,N, nz,  imin,imax,jmin,jmax, kmin, i,j,k,m

      integer km(0:Lm+1)

      real*8 var(imin:Lm+1,jmin:Mm+1,kmin:N), 
     &       z_r(0:Lm+1,0:Mm+1,N),   rmask(0:Lm+1,0:Mm+1),
     &       z_w(0:Lm+1,0:Mm+1,0:N), z_lev(imin:Lm+1,jmin:Mm+1,nz),
     &       FillValue, var_zlv(imin:Lm+1,jmin:Mm+1,nz),
     &       zz(0:Lm+1,0:N+1), dpth
#ifdef CUBIC
     &     , dz(0:Lm+1,kmin-1:N), FC(0:Lm+1,kmin-1:N), p,q,cff
#endif
      integer numthreads, trd, chunk_size, margin, jstr,jend
C$    integer omp_get_num_threads, omp_get_thread_num


      imax=Lm+1
      jmax=Mm+1

      numthreads=1
C$    numthreads=omp_get_num_threads()
      trd=0
C$    trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=jmin !max(  trd   *chunk_size  -margin, jmin )
      jend=jmax !min( (trd+1)*chunk_size-1-margin, jmax )


Cf2py intent(in) Lm,Mm,N, nz, z_r, z_w, rmask, var, z_lev, imin,jmin,kmin, FillValue
Cf2py intent(out) var_zlv


#ifdef VERBOSE
      write(*,'(2x,2A,3I3,2x,A,3I4)') 'sigma_to_z_intr: imin,jmin,',
     &   'kmin =', imin,jmin,kmin, 'trd,jstr,jend =', trd,jstr,jend
      write(*,*) Lm,Mm,N
#endif


      do j=jstr,jend
        if (kmin.eq.1) then
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=z_r(i,j,k)
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=z_w(i,j,0)
              zz(i,N+1)=z_w(i,j,N)
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.5D0*(z_r(i,j,k)+z_r(i-1,j,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.5D0*(z_w(i-1,j,0)+z_w(i,j,0))
              zz(i,N+1)=0.5D0*(z_w(i-1,j,N)+z_w(i,j,N))
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.5*(z_r(i,j,k)+z_r(i,j-1,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.5D0*(z_w(i,j,0)+z_w(i,j-1,0))
              zz(i,N+1)=0.5D0*(z_w(i,j,N)+z_w(i,j-1,N))
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.25D0*( z_r(i,j,k)+z_r(i-1,j,k)
     &                      +z_r(i,j-1,k)+z_r(i-1,j-1,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.25D0*( z_w(i,j,0)+z_w(i-1,j,0)
     &                     +z_w(i,j-1,0)+z_w(i-1,j-1,0))

              zz(i,N+1)=0.25D0*( z_w(i,j,N)+z_w(i-1,j,N)
     &                      +z_w(i,j-1,N)+z_w(i-1,j-1,N))
             enddo
          endif
        else
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=z_w(i,j,k)
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.5D0*(z_w(i,j,k)+z_w(i-1,j,k))
              enddo
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.5*(z_w(i,j,k)+z_w(i,j-1,k))
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.25D0*( z_w(i,j,k)+z_w(i-1,j,k)
     &                      +z_w(i,j-1,k)+z_w(i-1,j-1,k))
              enddo
            enddo
          endif
        endif
#ifdef CUBIC
        do k=kmin,N-1
          do i=imin,imax
            dz(i,k)=zz(i,k+1)-zz(i,k)
            FC(i,k)=var(i,j,k+1)-var(i,j,k)
          enddo
        enddo
        do i=imin,imax
          dz(i,kmin-1)=dz(i,kmin)
          FC(i,kmin-1)=FC(i,kmin)

          dz(i,N)=dz(i,N-1)
          FC(i,N)=FC(i,N-1)
        enddo
        do k=N,kmin,-1   !--> irreversible
          do i=imin,imax
            cff=FC(i,k)*FC(i,k-1)
            if (cff.gt.0.D0) then
              FC(i,k)=cff*(dz(i,k)+dz(i,k-1))/( (FC(i,k)+FC(i,k-1))
     &                                          *dz(i,k)*dz(i,k-1) )
            else
              FC(i,k)=0.D0
            endif
          enddo
        enddo
#endif
        do m=1,nz
       

          if (kmin.eq.0) then                      !
            do i=imin,imax                         !
              dpth=zz(i,N)-zz(i,0) 
              if (rmask(i,j).lt.0.5) then
                km(i)=-3     !--> masked out
              elseif (dpth*(z_lev(i,j,m)-zz(i,N)).gt.0.) then
                km(i)=N+2    !<-- above surface
              elseif (dpth*(zz(i,0)-z_lev(i,j,m)).gt.0.) then
                km(i)=-2     !<-- below bottom
              else
                km(i)=-1     !--> to search
              endif
            enddo
          else
            do i=imin,imax
              dpth=zz(i,N+1)-zz(i,0)     
              if (rmask(i,j).lt.0.5) then
                km(i)=-3     !--> masked out
              elseif (dpth*(z_lev(i,j,m)-zz(i,N+1)).gt.0.) then
                km(i)=N+2    !<-- above surface

              elseif (dpth*(z_lev(i,j,m)-zz(i,N)).gt.0.) then
                km(i)=N      !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(i,0)-z_lev(i,j,m)).gt.0.) then
                km(i)=-2     !<-- below bottom
              elseif (dpth*(zz(i,1)-z_lev(i,j,m)).gt.0.) then
                km(i)=0      !<-- above bottom, but below z_r(1)
              else
                km(i)=-1     !--> to search
              endif
            enddo
          endif
          do k=N-1,kmin,-1
            do i=imin,imax
              if (km(i).eq.-1) then
                if((zz(i,k+1)-z_lev(i,j,m))*(z_lev(i,j,m)-zz(i,k))
     & .ge. 0.) km(i)=k
              endif
            enddo
          enddo

          do i=imin,imax
            if (km(i).eq.-3) then
              var_zlv(i,j,m)=0.             !<-- masked out
            elseif (km(i).eq.-2) then
              var_zlv(i,j,m)=FillValue      !<-- below bottom
            elseif (km(i).eq.N+2) then
              var_zlv(i,j,m)=-FillValue     !<-- above surface
            elseif (km(i).eq.N) then
              var_zlv(i,j,m)=var(i,j,N)     !-> R-point, above z_r(N)
#ifdef CUBIC
     &              +FC(i,N)*(z_lev(i,j,m)-zz(i,N))
#else
     &              +(z_lev(i,j,m)-zz(i,N))*(var(i,j,N)-var(i,j,N-1))
     &                                   /(zz(i,N)-zz(i,N-1))
#endif
            elseif (km(i).eq.kmin-1) then   !-> R-point below z_r(1),
              var_zlv(i,j,m)=var(i,j,kmin)  !     but above bottom  
#ifdef CUBIC
     &           -FC(i,kmin)*(zz(i,kmin)-z_lev(i,j,m))
#else
     &          -(zz(i,kmin)-z_lev(i,j,m))*(var(i,j,kmin+1)-var(i,j,kmin))
     &                                  /(zz(i,kmin+1)-zz(i,kmin))
#endif
            else
              k=km(i)
              !write(*,*) k,km
#ifdef CUBIC
              cff=1.D0/(zz(i,k+1)-zz(i,k))
              p=z_lev(i,j,m)-zz(i,k)
              q=zz(i,k+1)-z_lev(i,j,m)

              var_zlv(i,j,m)=cff*( q*var(i,j,k) + p*var(i,j,k+1)
     &             -cff*p*q*( cff*(q-p)*(var(i,j,k+1)-var(i,j,k))
     &                                  +p*FC(i,k+1) -q*FC(i,k) )
     &                                                          )
#else
              var_zlv(i,j,m)=( var(i,j,k)*(zz(i,k+1)-z_lev(i,j,m))
     &                        +var(i,j,k+1)*(z_lev(i,j,m)-zz(i,k))
     &                               )/(zz(i,k+1)-zz(i,k))
#endif


            !write(*,*) 'bof',i,j,k,zz(i,k), zz(i,k+1), z_lev(i,j,m), m



#ifdef VERBOSE
              if ((zz(i,k+1)-z_lev(i,j,m))*(z_lev(i,j,m)-zz(i,k))
     & .lt.0.D0) then

                write(*,*) zz(i,k), zz(i,k+1), z_lev(i,j,m), m
                write(*,'(/1x,2A/)') '### ERROR: sigma_to_z_intr ',
     &                               ':: Search algoritm failure.'
                stop 
              endif
#endif
            endif
          enddo
        enddo ! <-- m
      enddo !<-- j

      return
      end
 
