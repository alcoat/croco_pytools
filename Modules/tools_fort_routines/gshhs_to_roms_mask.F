! Reads GSHHS shoreline data set and generates "raw" ROMS land mask
! placing it into a separate netCDF file.
! Usage:
!                gshhs_to_roms_mask roms_grid_file.nc
! followed by
!                copymask mask.nc roms_grid_file.nc
!
! followed by some inspection, perhaps editing by "editmask" matlab
! tool, and/or
!                single_connect i0 j0 roms_grid_file.nc
!
! where "copymask" and "single_connect" are in this package.

! input: (1) "lon_rho","lat_rho" and their dimensions read from
!                                           roms_grid_file.nc, and

!        (2) GSHHS shoreline data, usually named "gshhs_f.b" from
!
!            http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/latest/
!
!            which can be found inside  gshhg-bin-VERSION.zip archive
!            Because shoreline data comes in a very specific format
!            [binary, big-endian, integer-only, written sequentially
!            by a C-program in records with unpredictable length],
!            specialized programs are needed to read it.  Both C- and
!            Fortran versions, "read_gshhs.c" and "read_gshhs_data.F",
!            provided with this package can be used interchangeably
!            with equivalent outcome.  Within this file the selection
!            between then is done by CPP-switch USE_C_READER (defined
!            for C-, undefined for Fortran). Because this program
!            uses module "mod_gshhs" defined in "read_gshhs_data.F",
!            the latter needs to be always included into compilation
!            line in even in the case when C-version of reader and
!            decoder is used. Thus, to avoid naming conflict, the
!            reading and decoding functions in C- and Fortran versions
!            have different names, resulting in appearance of a few
!            extra #ifdefs within this file.  See also "gshhs_test.F".
!
!            IMPORTANT: File "gshhs_f.b" must be placed somewhere on
!            the local computer, after which edit "read_gshhs_data.F"
!            and find lines which look like
!
!                     #ifdef CRUDE
!                        fname='/path/to/dir/gshhs_i.b'
!                     #else
!                        fname='/path/to/dir/gshhs_f.b'
!                     #endif
!
!            point to the correct location. Similarly edit C-program
!            "read_gshhs.c" and make sure that lines
!
!                     char fname []="/path/to/dir/gshhs_f.b";
!
!            are set correctly as well.
!
!           [CPP-switch CRUDE above is set by Makefile to create an
!            alternative executable "gshhs_to_roms_mask_crude" which
!            is just a faster version due to reduced dataset -- with
!            loss of quality of the outcome accordingly. Nevertheless
!            it may be useful if grid building procedure involves
!            iterations (such as adjusting boundaries of the model
!            domain to minimize land masking). In this case using
!            crude version during intermediate stages may speed up
!            the process significantly.

! output: file "mask.nc" created by this program; it contains an
!         integer variable "land_rho" with dimensions matching the
!         above, and land_rho=0 for water points, nonzero-value(s)
!         (typically +2) for land, but may produce different values
!         for different land masses and object, depending on tuning
!         of this program below.

! This program uses explicit OpenMP threads with dynamic scheduling.

! Created and maintained by Alexander Shchepetkin, old_galaxy@yahoo.com


! WARNING: PRE_SELECT option designed for acceleration by reconnecting
! contours leaving model domain via its perimeter (hence reducing the
! number of points involves) is not fully developed and should not
! be activated until further notice.

      subroutine gshhs_to_roms_mask(ic,jc,xr,yr,fname,mask)

#define MERIDIAN
#define FAST_MODE
#define USE_C_READER
c-#define PRE_SELECT
#define TIMING
c--#define VERBOSE

      implicit none
      character(len=100) fname
      real(kind=8), dimension(ic+2,jc+2) :: xr,yr
      integer(kind=4), dimension(ic+2,jc+2) :: xc,yc
      integer(kind=2), dimension(ic+2,jc+2) :: mask

      integer :: header_size
      integer(kind=4) nsegm, npts_total
      integer(kind=4), dimension(:), allocatable :: isegm
      integer(kind=4), dimension(:,:), allocatable :: segmhdr,segmpt
#ifdef PRE_SELECT
      integer nrim
      integer(kind=4), dimension(:), allocatable :: xrim,yrim
#endif
      integer nargs, iargc, ncgrd, ncid, varid, ierr, nx,ny, part(4),
     &        west_edge, east_edge, south_edge, north_edge, r2dgrd(2),
     &                                  ic,jc, k, lgrd, ltrg
#ifdef USE_C_READER
     &                                                  , read_gshhs
#endif
      include "netcdf.inc"
#ifdef TIMING
      integer(kind=4) iclk(2), nclk, clk_rate, clk_max
      integer(kind=8) inc_clk, net_clk
      nclk=1 ; inc_clk=0 ; net_clk=0
      call system_clock(iclk(nclk), clk_rate, clk_max)
#endif
      header_size=11

Cf2py intent(in) ic,jc,xr,yr,fname
Cf2py intent(out) mask

! Read grid file first: since it will not be modified, open it with
! no-write status; because the intent is to convers double-precision
! number into integer, allocate integer arrays first, then double,
! and deallocate latter as soon as conversion is complete; save
! partition attribute (if exists): this may be useful for processing
! very large grids by parts and joining the output files.

      nx=ic+2 ; ny=jc+2
      part(1)=-1

      west_edge= 361000000 ; east_edge =-361000000
      south_edge=361000000 ; north_edge=-361000000
      do jc=1,ny
        do ic=1,nx
          xc(ic,jc)=nint(xr(ic,jc)*1.0D+6)   ! convert to
          yc(ic,jc)=nint(yr(ic,jc)*1.0D+6)   ! integer numbers

          west_edge=min(west_edge,xc(ic,jc))
          east_edge=max(east_edge,xc(ic,jc))
          south_edge=min(south_edge,yc(ic,jc))
          north_edge=max(north_edge,yc(ic,jc))
        enddo
      enddo

      write(*,*) 'Geographical limits of ROMS grid:'
      write(*,*) '  west =', west_edge*1.D-6,
     &           '  east =', east_edge*1.D-6
      write(*,*) ' south =', south_edge*1.D-6,
     &           ' north =', north_edge*1.D-6


#ifdef PRE_SELECT
      nrim=2*(nx+ny) ; allocate(xrim(0:nrim),yrim(0:nrim))
      k=0
      do ic=1,nx
        k=k+1
        xrim(k)=nint( (2.D0*xr(ic,1)-xr(ic,2))*1.0D+6 )
        yrim(k)=nint( (2.D0*yr(ic,1)-yr(ic,2))*1.0D+6 )
      enddo
      do jc=1,ny
        k=k+1
        xrim(k)=nint( (2.D0*xr(nx,jc)-xr(nx-1,jc))*1.0D+6 )
        yrim(k)=nint( (2.D0*yr(nx,jc)-yr(nx-1,jc))*1.0D+6 )
      enddo
      do ic=nx,1,-1
        k=k+1
        xrim(k)=nint( (2.D0*xr(ic,ny)-xr(ic,ny-1))*1.0D+6 )
        yrim(k)=nint( (2.D0*yr(ic,ny)-yr(ic,ny-1))*1.0D+6 )
      enddo
      do jc=ny,1,-1
        k=k+1
        xrim(k)=nint( (2.D0*xr(1,jc)-xr(2,jc))*1.0D+6 )
        yrim(k)=nint( (2.D0*yr(1,jc)-yr(2,jc))*1.0D+6 )
      enddo
      xrim(0)=xrim(k) ; yrim(0)=yrim(k)
#endif
!      deallocate(xr,yr)
!        else
!          write(*,'(/1x,4A/12x,A)')  '### ERROR: Cannot open netCDF ',
!     &         'file ''', grid_file(1:lgrd), '''.', nf_strerror(ierr)
!          stop
!        endif
!      else
!        write(*,'(/1x,A/12x,2A/)') 'Usage:', 'gshhs_to_roms_mask ',
!     &                                        'roms_grid_file.nc'
!        stop
!      endif


! Read GSHHS shoreline dataset...

#ifdef USE_C_READER
      nsegm=200000 ; npts_total=11000000
      allocate(isegm(nsegm), segmhdr(header_size,nsegm),
     &                            segmpt(2,npts_total))
      nsegm=0
      write(*,*) 'fname:',trim(fname),len(trim(fname))
      ierr=read_gshhs(nsegm, isegm, segmhdr,segmpt,trim(fname))
      if (ierr /= 0) stop
#else
      call read_gshhs_data(fname)
#endif


#ifdef TIMING
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate, clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(1x,2A,F5.2,1x,A)')   'reading of GSHHS dataset ',
     &        'complete in', dble(inc_clk)/dble(clk_rate), 'sec'
#endif

      call check_sanity(nsegm, npts_total, isegm, segmhdr, segmpt)

#ifdef PRE_SELECT
      call pre_select(nsegm, npts_total, isegm, segmhdr, segmpt,
     &                                             nrim,xrim,yrim)
#endif
!
! Generate land mask....
!
!      allocate(mask(nx,ny))
      do jc=1,ny
        do ic=1,nx
          mask(ic,jc)=0 !<-- initialize
        enddo
      enddo
!
      write(*,'(/1x,A)') 'start building land mask...'
C$OMP PARALLEL SHARED(nsegm, npts_total, isegm, segmhdr, segmpt,
C$OMP&                                        nx,ny, xc,yc, mask)
      call gshhs_to_mask_thread(nsegm, npts_total, isegm, segmhdr,
     &                                segmpt, nx,ny, xc, yc,mask)
C$OMP END PARALLEL


      end


      subroutine check_sanity(nsegm,npts_total, isegm,segmhdr,segmpt)
      implicit none
      integer nsegm, npts_total
      integer isegm(nsegm), segmhdr(11,nsegm), segmpt(2,npts_total)
      integer i,j,m, istr,iend, ierr,jj
c*    integer west,east

! Sanity check: isegm(m) is starting index "m"-th segment,
! segmhdr(2,m) is number of points in it;

      ierr=0
      do m=1,nsegm-1
        j=isegm(m+1)-isegm(m)-segmhdr(2,m)  !<-- should always be 0
        if (j /= 0) then
          ierr=ierr+1
          write(*,*) 'check_sanity :: ambiguous length of polygon',
     &                     segmhdr(1,m), isegm(m), segmhdr(2,m), j
        endif
      enddo
      if (ierr == 0) write(*,'(/1x,2A)') 'check_sanity :: passed ',
     &         'segment start, end, and length consistency check.'
      if (ierr /= 0) stop

      do m=1,nsegm
        istr=isegm(m) ; iend=istr+segmhdr(2,m)-1
        if ( segmpt(1,istr) /= segmpt(1,iend) .or.
     &       segmpt(2,istr) /= segmpt(2,iend) ) then
          if ( segmpt(1,istr) ==  180000000    .and.
     &         segmpt(1,iend) == -180000000      .and.
     &         segmpt(2,istr) == segmpt(2,iend) ) then
            write(*,'(1x,2A,I8)') 'check_sanity :: found ',
     &      '360-degree periodicity for polygon', segmhdr(1,m)
          else
            ierr=ierr+1
            write(*,*) 'check_sanity :: found non-closed polygon',
     &                                  segmhdr(1,m), istr,iend
            write(*,*)   segmpt(1,istr),segmpt(1,iend),
     &                   segmpt(2,istr),segmpt(2,iend)
          endif
        endif
c*      write(*,'(7I10)') m,istr,iend, segmpt(1,istr),segmpt(1,iend)
c*   &                               , segmpt(2,istr),segmpt(2,iend)
      enddo
      if (ierr == 0) write(*,*) 'check_sanity :: passed segment ',
     &                                 'start-end closure check.'
      if (ierr /= 0) stop

c*      do m=1,nsegm
c*        west=segmhdr(4,m) ; east=segmhdr(5,m)
c*        istr=isegm(m) ; iend=istr+segmhdr(2,m)-1
c*        do i=istr,iend
c*          if (segmpt(1,i) > east) then
c*            segmpt(1,i)=segmpt(1,i)-360000000
c*          endif
c*          if (segmpt(1,i) < west) then
c*            segmpt(1,i)=segmpt(1,i)+360000000
c*          endif
c*        enddo
c*      enddo

      do m=1,nsegm
        istr=isegm(m) ; iend=istr+segmhdr(2,m)-1
        ierr=0
        jj=0
!  1     continue
        do while(.true.) ! --> to avoid the goto
        jj=jj+1
        do i=istr+1,iend
          if ( segmpt(1,i)-segmpt(1,i-1) > 180000000 .or.
     &         segmpt(1,i-1)-segmpt(1,i) > 180000000 ) then
            ierr=+1

#ifdef VERBOSE
            write(*,'(17x,A,I11,1x,A,I8,1x,A,I12)')
     &              'increment', segmpt(1,i)-segmpt(1,i-1),
     &              'polygon',  segmhdr(1,m), 'at i =', i
#endif
          endif
        enddo
        
        if (ierr == +1) then
          write(*,*) m,i,ierr,jj
#ifdef VERBOSE
          write(*,'(1x,A,I8,1x,A)') 'check_sanity :: polygon',
     &    segmhdr(1,m), 'is discontinuous. Attempting to repair.'
#endif
          do i=istr+1,iend
            if (segmpt(1,i)-segmpt(1,i-1) > 180000000) then
                segmpt(1,i)=segmpt(1,i)-360000000
            elseif (segmpt(1,i-1)-segmpt(1,i) > 180000000) then
                    segmpt(1,i)=segmpt(1,i)+360000000
             endif
          enddo
          ierr=-1
!          goto 1 
        elseif (ierr == -1) then
          write(*,*) m,i,ierr,jj
          write(*,'(1x,2A,I8)') 'check_sanity :: repaired ',
     &               'continuity for polygon', segmhdr(1,m) 
         exit
        else
         exit
        endif
   
        enddo !!---> while loop
      enddo
      write(*,'(1x,A/)') 'check_sanity :: inspection complete.'
      end


#ifdef PRE_SELECT
      subroutine pre_select(nsegm,npts_total, isegm,segmhdr,segmpt,
     &                                            nrim, xrim,yrim)
      implicit none
      integer nsegm, npts_total, isegm(nsegm), segmhdr(11,nsegm),
     &    segmpt(2,npts_total), nrim, xrim(0:nrim), yrim(0:nrim)
!>
      integer xs,ys, west_edge, east_edge, south_edge, north_edge,
     &        shift, west, east, south, north, wstp,wstm, estp,estm,
     &   i,istr,iend, j,k,m, inc, parity, xntry,yntry,xexit,yexit,
     &   npts_pre_select, npts_reduced,  count, ncrs,jcrs(256)

      integer x_saved,y_saved
      logical found_overlap, start_outside, first_entry,
     &        saved_entry,   this_inside,   prev_inside
      real(kind=8) crss
      integer nbffr, max_segm_size
      integer, dimension(:,:), allocatable :: bffr

      max_segm_size=segmhdr(2,1)
      do m=2,nsegm
        max_segm_size=max(max_segm_size,segmhdr(2,m))
      enddo
      max_segm_size=max_segm_size+4*nrim
      allocate(bffr(2,max_segm_size))

      write(*,*) 'pre_select :: finding limits of the domain'

      west_edge= 361000000 ; east_edge =-361000000
      south_edge=361000000 ; north_edge=-361000000
      do i=1,nrim
        west_edge =min(west_edge ,xrim(i))
        east_edge =max(east_edge ,xrim(i))
        south_edge=min(south_edge,yrim(i))
        north_edge=max(north_edge,yrim(i))
      enddo

      write(*,*) 'pre_select :: start processing segments'

      npts_pre_select=0 ; npts_reduced=0
      do m=1,nsegm
        nbffr=0 ; found_overlap=.false.
        south=segmhdr(6,m) ; north=segmhdr(7,m)
        if (south <= north_edge .and. north >= south_edge) then
          west=segmhdr(4,m)   ; east=segmhdr(5,m)
          wstp=west+360000000 ; estp=east+360000000
          wstm=west-360000000 ; estm=east-360000000
          if (west <= east_edge .and. east >= west_edge) then
            found_overlap=.true. ; shift=0
          elseif(wstp <= east_edge .and. estp >= west_edge) then
            found_overlap=.true. ; shift=360000000
          elseif(wstm <= east_edge .and. estm >= west_edge) then
            found_overlap=.true. ; shift=-360000000
          endif
        endif
        if (found_overlap) then
          npts_pre_select=npts_pre_select+segmhdr(2,m)

          istr=isegm(m) ; iend=istr+segmhdr(2,m)-1
          prev_inside=.false.
          do i=istr,iend
            segmpt(1,i)=segmpt(1,i)+shift
            xs=segmpt(1,i) ; ys=segmpt(2,i)
!--->
            ncrs=0 ; count=0
            do j=1,nrim
              if ( yrim(j-1) <= ys .and. ys <= yrim(j) .or.
     &             yrim(j-1) >= ys .and. ys >= yrim(j) ) then
                ncrs=ncrs+1 ; jcrs(ncrs)=j
              endif
            enddo
            do k=1,ncrs
              j=jcrs(k)
              if (yrim(j) /= yrim(j-1)) then
                crss=( dble(xrim(j-1))*(dble(yrim(j))-dble(ys))
     &                +dble(xrim(j))*(dble(ys)-dble(yrim(j-1)))
     &                          )/(dble(yrim(j))-dble(yrim(j-1)))
                if (crss > dble(xs)) then
                  if (ys == yrim(j) .or. ys == yrim(j-1)) then
                    inc=1
                  else
                    inc=2
                  endif
                  if (yrim(j) > yrim(j-1)) then
                    count=count+inc
                  elseif (yrim(j) < yrim(j-1)) then
                    count=count-inc
                  endif
                endif
              endif
            enddo
            if (count == 2) then
              this_inside=.true.
            elseif (count == 0) then
              this_inside=.false.
            else
              write(*,*) '### ERROR: pre_select :: 1' ; stop
            endif
            if (i == istr) then
              parity=0
              if (this_inside) then
                start_outside=.false.
                write(*,*) 'starting inside'
              else
                start_outside=.true.
                write(*,*) 'starting outside'
              endif
              first_entry=.true.
              saved_entry=.false. !<-- trigger
            endif
!-->
            if (this_inside) then
              if (.not.prev_inside .and. i > istr) then
                parity=parity+1 ; xntry=xs ; yntry=ys
             write(*,*) 'entering rim', segmhdr(1,m), i
                if (first_entry .and. start_outside) then
                  first_entry=.false. ; x_saved=xs
                  saved_entry=.true.  ; y_saved=ys
                else
                  call close_contour(xexit,yexit, xntry,yntry,
     &                            nrim, xrim,yrim, nbffr,bffr)
                endif
              endif
              nbffr=nbffr+1
              bffr(1,nbffr)=xs ; bffr(2,nbffr)=ys
            elseif (prev_inside) then
              parity=parity-1 ; xexit=xs ; yexit=ys
           write(*,*) ' leaving rim', segmhdr(1,m), i
            endif
            prev_inside=this_inside
          enddo  !<-- i=istr,iend
          if (saved_entry) then
            write(*,*) 'closing saved entry'
            call close_contour(xexit,yexit, x_saved,y_saved,
     &                          nrim, xrim,yrim, nbffr,bffr)
          endif
          if (parity /= 0) then
            write(*,*) '### ERROR: parity of crossings ???'
          endif
        endif   !<-- found_overlap
        npts_reduced=npts_reduced+nbffr

        write(*,*) 'segment', segmhdr(1,m), segmhdr(2,m) , nbffr

        if (nbffr < segmhdr(2,m)) then
          segmhdr(2,m)=nbffr
          if (nbffr > 0) then
            istr=isegm(m) ; iend=istr+segmhdr(2,m)-1
            do i=istr,iend
              segmpt(1,i)=bffr(1,i-istr+1)
              segmpt(2,i)=bffr(2,i-istr+1)
            enddo
          endif
        elseif (nbffr > segmhdr(2,m)) then
          write(*,*) 'increased size for segment', segmhdr(1,m)
        endif
      enddo  !<-- m

      write(*,*) 'total number of points in dataset', npts_total
      write(*,*) '              pre-selected points', npts_pre_select
      write(*,*) '            after contour surgery', npts_reduced
      deallocate(bffr)
c       stop
      end

      subroutine close_contour(x1,y1, x2,y2, nrim, xrim,yrim,
     &                                           nbffr, bffr)
      implicit none
      integer x1,y1, x2,y2, nrim, xrim(0:nrim), yrim(0:nrim),
     &                                      nbffr, bffr(2,*)
      integer i, i1,i2
      call find_closest(x1,y1,i1, nrim, xrim,yrim)
      call find_closest(x2,y2,i2, nrim, xrim,yrim)
      write(*,*) 'closing contour, i1,i2 =',i1,i2
      if (i1 <= i2) then
        do i=i1,i2,+1
          nbffr=nbffr+1
          bffr(1,nbffr)=xrim(i) ; bffr(2,nbffr)=yrim(i)
        enddo
      else
        do i=i1,i2,-1
          nbffr=nbffr+1
          bffr(1,nbffr)=xrim(i) ; bffr(2,nbffr)=yrim(i)
        enddo
!/*
!        do i=i1,nrim
!          nbffr=nbffr+1
!          bffr(1,nbffr)=xrim(i) ; bffr(2,nbffr)=yrim(i)
!        enddo
!        do i=1,i1
!          nbffr=nbffr+1
!          bffr(1,nbffr)=xrim(i) ; bffr(2,nbffr)=yrim(i)
!        enddo
!*/
      endif
      end

      subroutine find_closest(xs,ys, icls, nrim, xrim,yrim)
      implicit none
      integer xs,ys, icls, nrim, xrim(0:nrim), yrim(0:nrim), i
      real(kind=8) min_dist,dist
      icls=1
      min_dist=(dble(xs)-dble(xrim(1)))**2
     &        +(dble(ys)-dble(yrim(1)))**2
      do i=1,nrim
        dist=(dble(xs)-dble(xrim(i)))**2
     &      +(dble(ys)-dble(yrim(i)))**2
        if (dist < min_dist) then
          min_dist=dist ; icls=i
        endif
      enddo
      end

      subroutine find_intsec(x1,y1, x2,y2, x3,y3, x4,y4, x,y)

! Finds coordinates of intersection of diagonals of a quadrilateral
! polygon (x1,y1), (x2,y2), (x3,y3), (x4,y4) labelled consecutively.

      implicit none
      integer(kind=4) x1,y1,x2,y2,x3,y3,x4,y4           !  4 <-- 3
      real(kind=8) x,y, a11,a12,a21,a22, d1,d2, det     !  !  X  !
                                                        !  1 --> 2
      a11=dble(y1-y3) ; a12=dble(x3-x1)
      a21=dble(y2-y4) ; a22=dble(x4-x2)

      d1=dble(x3)*dble(y1)-dble(x1)*dble(y3)
      d2=dble(x4)*dble(y2)-dble(x2)*dble(y4)

      det=1./(a11*a22 -a12*a21) ! Gabriel Cramer
      x=det *(d1 *a22 -a12*d2 ) ! rule of 1750
      y=det *(a11* d2 -d1 *a21)
      end
#endif

! Notes to core algorithm of "gshhs_to_mask_thread" below:

! (1) Longitudinal bounds of vertices GSHHS polygons are defined on
! polygon-to-polygon basis, rather than, say 0...360 or -180...+180 for
! all polygons. However, the definition is always the same (or intended
! to be the same, however see "check_sanity" above which occasionally
! finds some exclusions from this rule) for all vertices of the same
! polygon.  Thus, if the contour needs to cross Greenwich meridian or
! dateline (or both), longitude of the vertices of polygon always
! changes continuously through such crossings (no +360 or -360 degree
! jumps from one vertex to the next are allowed).  This also means that
! if both meridians are crossed, as in the case of Eurasian continent,
! some points of its coastline have longitude defined as -180...+180
! and some as 0...360.  It is for the reason of per-polygon definition
! point "xs=x1,x2,x3,x4" must be adjusted by +360 or -360 degrees
! individually for each polygon according to its extreme bounds, rather
! than do it once for all. Once the adjustment is made, no further
! logic is needed to take care about 360-degree periodicity.


! (2) The code below assumes that GSHHS polygons redundantly repeat
! their first and last points (hence the ranges for loops over vertices
! are istr+1:iend with i and i-1 used inside and no special care for
! loop ends).  Antarctica coastline is closed through periodicity: the
! contour starts at +180 degrees longitude, proceeds toward -180, so
! the first and the last point have the same latitude, but longitudes
! are +180 and -180 respectively, hence the contour is continuous from
! its first to its last point.

! (3) The algorithm is to count the number of crossings of the ray
! originating from the point (xs,ys) and going to the north, CPP-switch
! MERIDIAN is activated, or to the east, otherwise.  Each coastline
! segment is expected to be a closed (first and last points are the
! same) oriented contour enclosing land in counter-clock-wise direction
! (when going along the contour in the direction of index increase,
! land is expected to be on the left, water on the right).
! The crossings counted positive (+2 or +1) if, as seen from the point
! (xs,ys), the direction of increasing index is to the left (that is
! to the west as seen from the south), or negative (-2 or -1) is to
! the right (east). Increment 2 is for normal crossing when crossing
! strictly from the east side to strictly west (or vice versa);
! increment 1 is for half crossing, that is on of the edges happens
! to be exactly at the ray (in this case two or more consecutive
! segments are involved, and in the end is will add up to +2, -2,
! or 0). If everything goes normally as planned, all mask(ic,jc) will
! end up having values of either 0 or +2 (but not +1 or -1 or +4 or
! anything else).  However, arriving at values 0 or +2 may involve
! more than one count, e.g., a point belonging to a lake inside land
! is enclosed by two contours, one of which, continental coastline,
! be counted as +2 for enclosing the point in counter-clock-wise
! direction, and the other, lake coast, as -2 for clock-wise, so
! the point ends up +2-2=0, water.  If a point fall exactly on the
! contour, it is considered to be outside of that contour.
! Choosing to point the ray to the north (not, say, to the east,
! whatever) is motivated by allowing to handle Antarctica coastline
! without any special treatment whatsoever.

! (4) WARNING: The code below operates with large integers so only
! comparisons operations between x-ses and y-s are allowed.  Any
! multiplication of them will most likely cause overflow (e.g., using
! (xm-xs)*(xs-xp)<=0 as logical condition that "xs" falls somewhere
! between "xm" and "xp", while "xm" can be more, less, or equal to
! "xp" is very likely result in logical error because of integer
! overflow.

! (5) Full resolution GSHHS dataset contains approximately 190,000
! coastline polygons totalling to 11,000,000 vertices put together,
! so building mask for an (nx,ny)-dimensioned grid leads to overall
! 11,000,000*nx*ny operations -- basically each point on ROMS grid
! must be checked against all the polygons, whether it is inside or
! outside -- therefore for the reason of performance
!
!   (5a) reject polygons whose north-south and east-west bounds
!        (as deduced from GSHHS segment header) do not cover the
!        coordinates of point being evaluated;
!
!   (5b) the algorithm is split into two stages: at first, search
!        and just record indices of polygon segments which cross the
!        meridian, while deferring any further analysis to keep the
!        number of operations inside loops to absolute minimum;
!        secondly, process the recorded segments to count the number
!        of crossings (subroutine "count_crossings") to determine
!        whether the point is inside or outside of a polygon;
!
!   (5c) Manually unroll the search loops in (5b) to process 9 points
!        of ROMS grid arranged 3x3 square (CPP-switch FAST_MODE) at the
!        same time, essentially changing
!
!           do jc=1,ny     !<-- two outer loops are
!            do ic=1,ny    !<-- over ROMS grid points
!             do m=1,nsegm    !<-- loop over polygons in GSHHS dataset
!              do i=istr,iend         !<-- search loop over vertices
!                 check whether north-bound             of polygon "m"
!                 ray starting from point(ic,jc)
!                 crosses one or more edges of
!                 polygon "m" and record indices
!                 of the crossings;
!              enddo
!             enddo
!            enddo
!           enddo
!
!        into
!
!           do jc=1,ny,+3
!            do ic=1,ny,+3
!             do m=1,nsegm
!              do i=istr,iend
!                 * * *     simultaneously process nine
!                 * * *       points [ ic, ic+1, ic+2 ]
!                 * * *            x [ jc, jc+1, jc+2 ]
!              enddo
!             enddo
!            enddo
!           enddo
!
!        Mathematically speaking, to check whether an individual grid
!        point is enclosed by a polygon, the entire polygon (in the
!        worst case consisting on 1,000,000 vertices) must be walked
!        through, resulting in out-of-cache stalls, because there are
!        simply not enough operations in the search loop to offset the
!        cost of main memory reads associated with polygon vertices.
!        In practical terms, unrolling the search loop 4 times (in 2x2
!        arrangement) leads to a factor of 3 gain in speed of this
!        algorithm. Unrolling 9 times (in 3x3) leads to 4+ gain.
!
! (6) Multi-threading approach here involves splitting jc-loop of
! the grid into 3-point-wide strips which are processed by different
! threads using dynamic scheduling: each thread takes out the next
! available three rows of points, jc, jc+1, and jc+2, to work on and
! increments counter "comm_j_indx" by 3 (a shared variable being
! accessed and modified only inside a critical region) to signal the
! others that these rows is taken.  Once the thread completes its rows
! it goes back to the beginning of the do-while loop and checks out
! whatever there are next rows which are not taken by others, and so
! on until nothing left to process.  Timing variable "net_clk_bak" is
! needed to restrict printing out wall-clock time elapsed from the
! beginning of this subroutine, but no more often than once in
! approximately 5 seconds.

!      module comm_vars_gshhs_mask             ! shared variable
!       integer, save :: comm_j_indx=0         ! for synchronization
!      end module comm_vars_gshhs_mask         ! of Open MP threads

      subroutine gshhs_to_mask_thread(nsegm,npts_total, isegm,segmhdr,
     &                                     segmpt, nx,ny, xc, yc,mask)
!      use comm_vars_gshhs_mask
      implicit none
      integer :: comm_j_indx
      integer nsegm, npts_total, isegm(nsegm), segmhdr(11,nsegm),
     &        segmpt(2,npts_total), nx,ny,  xc(nx,ny), yc(nx,ny)
      integer(kind=2) mask(nx,ny)
!-->
      integer west,east,north,south, level, version, greenwich,
     &        source, river, area_scale, numthreads, trd, ic,jc,
     &                             xm,xp, m,i, istr,iend, count
      integer, parameter :: max_ncrs=256
      integer x1,y1, ncr1, icr1(max_ncrs)
#ifdef FAST_MODE
      integer x2,y2, ncr2, icr2(max_ncrs),
     &        x3,y3, ncr3, icr3(max_ncrs), ic1,
     &        x4,y4, ncr4, icr4(max_ncrs), ic2,
     &        x5,y5, ncr5, icr5(max_ncrs), jc1,
     &        x6,y6, ncr6, icr6(max_ncrs), jc2,
     &        x7,y7, ncr7, icr7(max_ncrs),
     &        x8,y8, ncr8, icr8(max_ncrs),
     &        x9,y9, ncr9, icr9(max_ncrs)
#endif
#ifndef MERIDIAN
      integer ym,yp
#endif
      real(kind=8) area
#ifdef USE_C_READER
      integer ierr, decode_flag !<-- function
#endif
C$    integer omp_get_num_threads, omp_get_thread_num
#ifdef TIMING
      integer(kind=4) iclk(2), nclk, clk_rate, clk_max
      integer(kind=8) :: inc_clk, net_clk, net_clk_bak
      real(kind=4) tstart, run_time
      call cpu_time(tstart)
      nclk=1 ; net_clk=0 ;  net_clk_bak=0
      call system_clock(iclk(nclk), clk_rate, clk_max)
#endif
      comm_j_indx=0
      trd=0 ; numthreads=1
C$    trd=omp_get_thread_num() ; numthreads=omp_get_num_threads()

      do while(.true.)
#ifdef TIMING
        nclk=3-nclk
        call system_clock(iclk(nclk), clk_rate,clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        if (inc_clk < 0) inc_clk=inc_clk+clk_max
        net_clk=net_clk+inc_clk
        if (net_clk-net_clk_bak > 5*clk_rate) net_clk_bak=net_clk
#endif
C$OMP CRITICAL(progress_sync)
        jc=comm_j_indx+1
#ifdef FAST_MODE
        comm_j_indx=jc+2 !<-- three rows at at time
#else
        comm_j_indx=jc
#endif
#ifdef TIMING
        if (trd == 0 .and. net_clk == net_clk_bak) then
         write(*,'(F14.2,2x,I6,2x,A,I6,4x,A,I3)')
     &         dble(net_clk)/dble(clk_rate), jc,'out of',ny,'trd =',trd
        else
#endif
         write(*,'(16x,I6,2x,A,I6,4x,A,I3)') jc,'out of',ny,'trd =',trd
#ifdef TIMING
        endif
#endif
C$OMP END CRITICAL(progress_sync)
        if (jc > ny) exit !goto 98         !--> DONE

#ifdef FAST_MODE
        jc1=min(jc+1,ny) ; jc2=min(jc+2,ny)   ! less than three rows
        do ic=1,nx,3                          ! may be left at the end
          ic1=min(ic+1,nx) ; ic2=min(ic+2,nx)
          y1=yc(ic,jc ) ; y2=yc(ic1,jc ) ; y3=yc(ic2,jc )
          y4=yc(ic,jc1) ; y5=yc(ic1,jc1) ; y6=yc(ic2,jc1)
          y7=yc(ic,jc2) ; y8=yc(ic1,jc2) ; y9=yc(ic2,jc2)
#else
        do ic=1,nx
          y1=yc(ic,jc)
#endif
          do m=1,nsegm
#ifdef USE_C_READER
            ierr=decode_flag(segmhdr(3,m),  level, version, greenwich,
#else
            call dec_segm_flag(segmhdr(3,m), level, version, greenwich,
#endif
     &                                    source,  river,  area_scale)
            area=dble(segmhdr(8,m))/dble(10**area_scale)
c**         if (level == 1 .and. area > 0.D0) then
              south=segmhdr(6,m) ; north=segmhdr(7,m)
#ifdef FAST_MODE
              if (south <= max(y1,y2,y3,y4, y6,y7,y8,y9)  .and.
     &                     min(y1,y2,y3,y4, y6,y7,y8,y9) <= north) then
#else
              if (south <= y1 .and. y1 <= north) then
#endif
                west=segmhdr(4,m) ; east=segmhdr(5,m)  ! Try to bring
                x1=xc(ic,jc)                           ! longitude of
                xm=x1-360000000 ; xp=x1+360000000      ! grid point
                if (west <= xm .and. xm <= east) x1=xm ! into polygon
                if (west <= xp .and. xp <= east) x1=xp ! bounds.
#ifdef FAST_MODE
                x2=xc(ic1,jc)
                xm=x2-360000000 ; xp=x2+360000000
                if (west <= xm .and. xm <= east) x2=xm
                if (west <= xp .and. xp <= east) x2=xp

                x3=xc(ic2,jc)
                xm=x3-360000000 ; xp=x3+360000000
                if (west <= xm .and. xm <= east) x3=xm
                if (west <= xp .and. xp <= east) x3=xp

                x4=xc(ic,jc1)
                xm=x4-360000000 ; xp=x4+360000000
                if (west <= xm .and. xm <= east) x4=xm
                if (west <= xp .and. xp <= east) x4=xp

                x5=xc(ic1,jc1)
                xm=x5-360000000 ; xp=x5+360000000
                if (west <= xm .and. xm <= east) x5=xm
                if (west <= xp .and. xp <= east) x5=xp

                x6=xc(ic2,jc1)
                xm=x6-360000000 ; xp=x6+360000000
                if (west <= xm .and. xm <= east) x6=xm
                if (west <= xp .and. xp <= east) x6=xp

                x7=xc(ic,jc2)
                xm=x7-360000000 ; xp=x7+360000000
                if (west <= xm .and. xm <= east) x7=xm
                if (west <= xp .and. xp <= east) x7=xp

                x8=xc(ic1,jc2)
                xm=x8-360000000 ; xp=x8+360000000
                if (west <= xm .and. xm <= east) x8=xm
                if (west <= xp .and. xp <= east) x8=xp

                x9=xc(ic2,jc2)
                xm=x9-360000000 ; xp=x9+360000000
                if (west <= xm .and. xm <= east) x9=xm
                if (west <= xp .and. xp <= east) x9=xp

                if (west <= max(x1,x2,x3,x4, x6,x7,x8,x9)  .and.
     &                      min(x1,x2,x3,x4, x6,x7,x8,x9) <= east) then
#else
                if (west <= x1 .and. x1 <= east) then
#endif
                  istr=isegm(m) ; iend=istr+segmhdr(2,m)-1
                  ncr1=0
#ifdef FAST_MODE
                  ncr2=0 ; ncr3=0 ; ncr4=0 ; ncr5=0
                  ncr6=0 ; ncr7=0 ; ncr8=0 ; ncr9=0
#endif
#ifdef MERIDIAN
                  xm=segmpt(1,istr)
                  do i=istr+1,iend,+1   !--> recursive
                    xp=segmpt(1,i)
                    if (xm <= x1 .and. x1 <= xp .or.
     &                  xp <= x1 .and. x1 <= xm) then
                       ncr1=ncr1+1 ; icr1(ncr1)=i
                    endif
# ifdef FAST_MODE
                    if (xm <= x2 .and. x2 <= xp .or.
     &                  xp <= x2 .and. x2 <= xm) then
                       ncr2=ncr2+1 ; icr2(ncr2)=i
                    endif
                    if (xm <= x3 .and. x3 <= xp .or.
     &                  xp <= x3 .and. x3 <= xm) then
                       ncr3=ncr3+1 ; icr3(ncr3)=i
                    endif
                    if (xm <= x4 .and. x4 <= xp .or.
     &                  xp <= x4 .and. x4 <= xm) then
                       ncr4=ncr4+1 ; icr4(ncr4)=i
                    endif
                    if (xm <= x5 .and. x5 <= xp .or.
     &                  xp <= x5 .and. x5 <= xm) then
                       ncr5=ncr5+1 ; icr5(ncr5)=i
                    endif
                    if (xm <= x6 .and. x6 <= xp .or.
     &                  xp <= x6 .and. x6 <= xm) then
                       ncr6=ncr6+1 ; icr6(ncr6)=i
                    endif
                    if (xm <= x7 .and. x7 <= xp .or.
     &                  xp <= x7 .and. x7 <= xm) then
                       ncr7=ncr7+1 ; icr7(ncr7)=i
                    endif
                    if (xm <= x8 .and. x8 <= xp .or.
     &                   xp <= x8 .and. x8 <= xm) then
                       ncr8=ncr8+1 ; icr8(ncr8)=i
                    endif
                    if (xm <= x9 .and. x9 <= xp .or.
     &                  xp <= x9 .and. x9 <= xm) then
                       ncr9=ncr9+1 ; icr9(ncr9)=i
                    endif
# endif
                    xm=xp               !<-- recursion
                  enddo
#else
                  ym=segmpt(2,istr)
                  do i=istr+1,iend,+1
                    yp=segmpt(2,i)
                    if (ym <= y1 .and. y1 <= yp .or.
     &                  yp <= y1 .and. y1 <= ym) then
                       ncr1=ncr1+1 ; icr1(ncr1)=i
                    endif
# ifdef FAST_MODE
                    if (ym <= y2 .and. y2 <= yp .or.
     &                  yp <= y2 .and. y2 <= ym) then
                       ncr2=ncr2+1 ; icr2(ncr2)=i
                    endif
                    if (ym <= y3 .and. y3 <= yp .or.
     &                  yp <= y3 .and. y3 <= ym) then
                       ncr3=ncr3+1 ; icr3(ncr3)=i
                    endif
                    if (ym <= y4 .and. y4 <= yp .or.
     &                  yp <= y4 .and. y4 <= ym) then
                       ncr4=ncr4+1 ; icr4(ncr4)=i
                    endif
                    if (ym <= y5 .and. y5 <= yp .or.
     &                  yp <= y5 .and. y5 <= ym) then
                       ncr5=ncr5+1 ; icr5(ncr5)=i
                    endif
                    if (ym <= y6 .and. y6 <= yp .or.
     &                  yp <= y6 .and. y6 <= ym) then
                       ncr6=ncr6+1 ; icr6(ncr6)=i
                    endif
                    if (ym <= y7 .and. y7 <= yp .or.
     &                  yp <= y7 .and. y7 <= ym) then
                       ncr7=ncr7+1 ; icr7(ncr7)=i
                    endif
                    if (ym <= y8 .and. y8 <= yp .or.
     &                  yp <= y8 .and. y8 <= ym) then
                       ncr8=ncr8+1 ; icr8(ncr8)=i
                    endif
                    if (ym <= y9 .and. y9 <= yp .or.
     &                  yp <= y9 .and. y9 <= ym) then
                       ncr9=ncr9+1 ; icr9(ncr9)=i
                    endif
# endif
                    ym=yp
                  enddo
#endif
                  if ( ncr1 > max_ncrs
#ifdef FAST_MODE
     &                  .or. ncr2 > max_ncrs .or. ncr3 > max_ncrs
     &                  .or. ncr4 > max_ncrs .or. ncr5 > max_ncrs
     &                  .or. ncr6 > max_ncrs .or. ncr7 > max_ncrs
     &                  .or. ncr8 > max_ncrs .or. ncr9 > max_ncrs
#endif
     &               ) then
                    write(*,*) '### ERROR: gshhs_to_mask_thread :: ',
     &                  'insufficient size of parameter max_ncrs =',
     &                                                     max_ncrs
                    stop
                  endif
                  if (ncr1 > 0) then
                    call count_crossings(npts_total, segmpt,
     &                                     x1,y1,ncr1,icr1, count)
                    if (count /= 0) mask(ic,jc)=mask(ic,jc)+count
                  endif
#ifdef FAST_MODE
                  if (ncr2 > 0 .and. ic < nx) then
                    call count_crossings(npts_total, segmpt,
     &                                       x2,y2,ncr2,icr2, count)
                    if (count /= 0) mask(ic1,jc)=mask(ic1,jc)+count
                  endif
                  if (ncr3 > 0 .and. ic1 < nx) then
                    call count_crossings(npts_total, segmpt,
     &                                       x3,y3,ncr3,icr3, count)
                    if (count /= 0) mask(ic2,jc)=mask(ic2,jc)+count
                  endif
                  if (ncr4 > 0               .and. jc < ny) then
                    call count_crossings(npts_total, segmpt,
     &                                       x4,y4,ncr4,icr4, count)
                    if (count /= 0) mask(ic,jc1)=mask(ic,jc1)+count
                  endif
                  if (ncr5 > 0 .and. ic < nx .and. jc < ny) then
                    call count_crossings(npts_total, segmpt,
     &                                        x5,y5, ncr5,icr5, count)
                    if (count /= 0) mask(ic1,jc1)=mask(ic1,jc1)+count
                  endif
                  if (ncr6 > 0 .and. ic1 < nx .and. jc < ny) then
                    call count_crossings(npts_total, segmpt,
     &                                         x6,y6,ncr6,icr6, count)
                    if (count /= 0) mask(ic2,jc1)=mask(ic2,jc1)+count
                  endif
                  if (ncr7 > 0                .and. jc1 < ny) then
                    call count_crossings(npts_total, segmpt,
     &                                       x7,y7,ncr7,icr7, count)
                    if (count /= 0) mask(ic,jc2)=mask(ic,jc2)+count
                  endif
                  if (ncr8 > 0 .and. ic < nx .and. jc1 < ny) then
                    call count_crossings(npts_total, segmpt,
     &                                      x8,y8,ncr8,icr8, count)
                    if (count /= 0) mask(ic1,jc2)=mask(ic1,jc2)+count
                  endif
                  if (ncr9 > 0 .and. ic1 < nx .and. jc1 < ny) then
                    call count_crossings(npts_total, segmpt,
     &                                         x9,y9,ncr9,icr9, count)
                    if (count /= 0) mask(ic2,jc2)=mask(ic2,jc2)+count
                  endif
#endif
                endif  !<-- east-west bounds
              endif  !<-- north-south bounds
c      endif
c**         endif  !<-- level, area
          enddo  !<-- m=1,nsegm
        enddo  !<-- ic
      enddo  !<-- while(.true.), jc-loop

!  98  continue
C$OMP BARRIER
#ifdef TIMING
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      call cpu_time(run_time) ; run_time=run_time-tstart
       if (trd == 0) then
        write(*,'(/2x,A,F11.2,1x,A/2x,A,I4,4x,A,F11.2,1x,A/)')
     & 'cpu time =',run_time, 'sec', 'number of threads =', numthreads,
     & 'elapsed wall-clock time =', dble(net_clk)/dble(clk_rate), 'sec'
      endif
#endif
      end


      subroutine count_crossings(npts_total,segmpt, xs,ys, ncrs,icrs,
     &                                                        count)
      implicit none
      integer npts_total, xs,ys, ncrs, count
      integer segmpt(2,npts_total), icrs(ncrs)
      integer i,k, xm,ym, xp,yp
      real(kind=8) crss

      count=0
      do k=1,ncrs
        i=icrs(k)
        xp=segmpt(1,i) ; xm=segmpt(1,i-1)
        yp=segmpt(2,i) ; ym=segmpt(2,i-1)

#ifdef MERIDIAN
        if (xp /= xm) then
          crss=( dble(ym)*(dble(xp)-dble(xs))
     &          +dble(yp)*(dble(xs)-dble(xm))
     &                  )/(dble(xp)-dble(xm))
# ifdef VERBOSE
          if ((dble(yp)-crss)*(crss-dble(ym)) < 0.D0) write(*,*)
     &      '### ERROR: Search algorithm failure at x,y=', xs,ys
# endif
          if (crss > dble(ys)) then
            if (xs == xm .or. xs == xp) then
              i=1
            else
              i=2
            endif
            if (xp < xm) then
              count=count+i
            elseif (xp > xm) then
              count=count-i
            endif
          endif
        endif
#else
        if (yp /= ym) then
          crss=( dble(xm)*(dble(yp)-dble(ys))
     &          +dble(xp)*(dble(ys)-dble(ym))
     &                  )/(dble(yp)-dble(ym))
# ifdef VERBOSE
          if ((dble(xp)-crss)*(crss-dble(xm)) < 0.D0) write(*,*)
     &      '### ERROR: Search algorithm failure at x,y=', xs,ys
# endif
          if (crss > dble(xs)) then
            if (ys == ym .or. ys == yp) then
              i=1
            else
              i=2
            endif
            if (yp > ym) then
              count=count+i
            elseif (yp < ym) then
              count=count-i
            endif
          endif
        endif
#endif
      enddo
      end
