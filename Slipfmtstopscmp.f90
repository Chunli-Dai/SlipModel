        Program Read_caltech
!       Modified from 
!~/Osuwork/GRACE/Seismic2gravityByLei/org_dai/psgrn+pscmp-2008a/Read_caltech.f90
!       Purpose: read slip format  :
!                transform the format to be the input for pscmp08a.exe
!       Modification: from version bp1, dcl, June 2014
!       Purpose: to be consistent with the pscmp08a.exe, output the
!                on the left top point instead of the center of each
!                patch.
!       See also: GRVD2FAULTSApll/workMw/Fault2M0Loc.f
!!!!       
!       Format WangRJ: read GPS_SlipModelA.dat sent by Wang RJ 
!                      Example:
!   lat_deg     lon_deg    depth_km  x_local_km  y_local_km   length_km
!   width_km  slp_strk_m  slp_ddip_m    slp_am_m  strike_deg     dip_deg
!   rake_deg sig_stk_MPa sig_ddi_MPa sig_nrm_MPa
!   ...
!!!!       
!       Format 1     : Caltech or USGS format, eg. Sumatra12/SM12Mw86CJi.txt
!                      Example:
!#Total number of fault_segments=     1
!#Fault_segment =   1 nx(Along-strike)=  25 Dx= 15.00km ny(downdip)=   8 Dy=  5.00km
!#Boundary of Fault_segment     1. EQ in cell 13,4. Lon: 93.0725   Lat: 2.3484
!#Lon.  Lat.  Depth
!       93.64040        3.92750        6.58910
!       92.55420        0.75240        6.58910
!       92.49720        0.77190       44.33420
!       93.58330        3.94700       44.33420
!       93.64040        3.92750        6.58910
!#Lat. Lon. depth slip rake strike dip t_rup t_rise mo
!...
!!!!       
!       Format 2     : format as in caltech/caltech.out 
!                      Example:
!% lat1     lon1      wid(km) length(km)  dep(km)  dip(degree) strike(degree) rake(degree) slip(cm)
!...
!!!!
!       Format 3     : format as in ../Sumatra/SM04Caltech.txt 
!                      Example:
!       Input slip model format:
!slip model of the Sumatra-Andaman Mw9.1 event from Chlieh et al. BSSA 2007
!
!Long.      Lat.    slip(cm)  rake     Azimuth    Dip  Depth(km)
!...
!!!!
!       Format 4     : Input slip model format: Ozawa_coseismic_fault_para.dat
!                      Example: 
!lat.          longi.       dislo(slip in meter)   slip(rake)   depth azim(strike) dip   Lenth  Width
!  42.00130  138.99920        0.00144      158.83607      205.59300   182.216    34.002   32.0   32.0
!  ...
!99999
!!!!       
!       Format GCMTndk:
!                      Example:/0/home/dai.56/share/Osuworks/UltraSH/Slepian/IFILES/CMT/jan76_feb10.ndk 
!!!!
!       Format SIV  :Finite source model from SIV project, http://equake-rc.info/SRCMOD/
!                       Example: Sumatra12/s2012SUMATR01YUEx.fsp
!%
!-------------------------------------------------------------------------------------------------- 
!%    LAT     LON     X==EW     Y==NS     Z     SLIP     RAKE    
!%
!-------------------------------------------------------------------------------------------------- 
!    2.0070   93.9122   94.6167  -33.6747    5.8510    0.0000  170.0000  
!    2.0566   93.7390   75.3915  -28.1619    5.8510    2.9440  170.0000  
!!!!
!       Format Vigny : Finite source model from Vigny,
!                      http://geologie.ens.fr/~vigny/MAULE/vigny_etal_2011_model.dat
!                     Example: Chile/vigny_etal_2011_model.dat
!# file format :
!# lon lat depth length-along-strike width-along-dip dip strike sslip(cm) dslip opening
!# where lon & lat   are at the lower left corner of each patch

        implicit none
        integer i,n,err
        double precision lat,lon,wid,length,dep,dip,strike,rake,slip
        double precision slip_strike,slip_downdip
        character string*200
        character ifile*80,ofile*80
        double precision pi
        double precision x,y,slip_am,tmp,slip_op
        double precision oceanthi,depmo
        integer NP
        integer ns,nx,ny,nt
        double precision dx,dy
        character fmttype*40
        integer j,k,faultflag,ie,cnt,ifs,nfs,ntype
        integer ilat,ilon,idep,islip,irake,istrike,idip,iwid,ilen
        character strp*20,strtmp*20
        character*20,allocatable:: strtype(:)
        double precision,allocatable :: dat(:)
        logical center
        double precision d2r
        character cmtcode*80, fname*160
        double precision mu,M0
        integer nsbfs,nh

        center=.false. ! true: output the center of each patch;
                      ! false: output the left top corner of each patch
        
        pi=datan(1d0)*4d0
        d2r=PI/180D0
        write(*,*)'pi=',pi
                                                                                           
        write(*,*)'Input ifile ofile name, format type, ocean thickness &
     & Ozawa_coseismic_fault_pscmp.dat 4 3.9433'
        read(*,*)ifile,ofile,fmttype,oceanthi
        write(*,*)trim(ifile),' ',trim(ofile),' ',trim(fmttype),' ', oceanthi

        if(trim(fmttype).eq.'WangRJ')then
!       write(*,*)'Ex: WangRJGPS/GPS_SlipModelA.dat GPS_SlipModelA_pscmp.out 3.9433'
          open(1,file=ifile)
          read(1,*)string
          i=0
          do while(.true.)
            read(1,fmt='(a200)', iostat=err)string
            if(err.lt.0)exit
            i=i+1
          end do
          nt=i
          write(*,*)'Number of lines:',nt
          rewind(1)
          do i=1,1
            read(1,fmt='(a200)')string
          enddo

          i=0
          open(2,file=ofile)
          write(2,fmt='(I8)')nt
          do while(.true.)
            read(1,fmt='(a200)', iostat=err)string
            if(err.lt.0)exit
            i=i+1
            read(string,fmt=*)lat,lon,dep,x,y,length,wid,slip_strike, &
     & slip_downdip,slip_am,strike,dip,rake,tmp,tmp,tmp
            slip_op=0d0

            if (.not.center)then !If output the top left of each patch
            ! compute the top left coordinates from center of each patch
              call center2topleft(length,wid,strike,dip,lat,lon,dep)
            endif

            depmo=dep-oceanthi
            if(depmo.le.0d0)then
              write(*,*)'Warning: Depth above ocean bottom, &
     &  force it to be below!'
              depmo=0.001d0
            endif

            write(2,fmt='(I8,7(1x,E23.15),2(1x,I3),1x,f3.1)')  &
                i,lat,lon,depmo,length,wid,strike,dip,1,1,0d0
            write(2,fmt='(8x,2(1x,E23.15),3(1x,E23.15))') &
                length/2d0,wid/2d0,slip_strike,slip_downdip,slip_op
          enddo
          close(1)
          close(2)

        elseif(trim(fmttype).eq.'1')then
!       Works for 1s1 1s2 1
! 1s1   write(*,*)'Ex: ../Sumatra12/SM12Mw86p1USGS.txt ../Sumatra12/SM12Mw86p1USGS_pscmp.out 4.0'
! 1s2   write(*,*)'Ex: ../Sumatra12/SM12Mw86CJi.txt ../Sumatra12/SM12Mw86CJi_pscmp.out 4.0'
! 1     write(*,*)'Ex: caltech/static_caltech.txt caltech/static_caltech_pscmp.out 3.9433'
          open(1,file=ifile)
          faultflag=0
          do while(.true.)
            read(1,fmt='(a200)', iostat=err)string
            if(err.lt.0)exit
            do i=1,80
              ie=i+30-1
              if (string(i:ie).eq.'Total number of fault_segments')then
                strp=string(ie+2:ie+2+6) !=     1
                read(strp,*)nfs
                faultflag=1
                exit
              endif
            enddo
            if (faultflag.eq.1) exit
          end do
          if(faultflag.eq.0)then
            write(*,*)'Error: found no fault_segments!'
          elseif (faultflag.eq.1)then
            !Read the fault_segment
            cnt=0
            i=0
            open(2,file=ofile)
            write(2,fmt='(I8)')0
            do ifs=1,nfs
              read(1,fmt='(a200)', iostat=err)string
              if(string(73:74).ne.'Dy') &
     & write(*,*)'Format might be read wrong'
              read(string,fmt='(16x,I4,18x,I4,4x,f6.2,15x,i4,4x,f6.2)')ns,nx,dx,ny,dy
              if(ns.ne.ifs)write(*,*)'Format might be read wrong'
              length=dx;wid=dy
              cnt=cnt+nx*ny
              do k=1,7
                read(1,fmt='(a200)')string
              enddo
              read(1,fmt='(a200)')string
              ntype=7
              allocate(strtype(ntype),dat(ntype))
              read(string(2:200),fmt=*)(strtype(k),k=1,ntype)
!c----- set order of the record types -----
!#Lat. Lon. depth slip rake strike dip t_rup t_rise mo
!#Lat.       Lon.     depth     slip(cm)     rake     strike   dip  
!#Lon. Lat. depth slip rake strike dip
              ilat=99
              ilon=99
              idep=99
              islip=99
              irake=99
              istrike=99
              idip=99
              do k=1,ntype
                strtmp=strtype(k)
                if(strtmp(1:3).eq."Lat") ilat = k
                if(strtmp(1:3).eq."Lon") ilon = k
                if(strtmp(1:5).eq."depth") idep = k
                if(strtmp(1:4).eq."slip") islip = k
                if(strtmp(1:4).eq."rake") irake = k
                if(strtmp(1:6).eq."strike") istrike = k
                if(strtmp(1:3).eq."dip") idip = k
              enddo
              
!             i=0
              do j=1,nx*ny
                read(1,fmt='(a200)', iostat=err)string
                if(err.lt.0)exit
                i=i+1
!               write(*,*)string
!               read(string,fmt=*)lat,lon,wid,length,dep,dip,strike,rake,slip
!               read(string,fmt=*)lon,lat,dep,slip,rake,strike,dip
!               read(string,fmt=*)lat,lon,dep,slip,rake,strike,dip
                read(string,fmt=*)(dat(k),k=1,ntype)
                lat=dat(ilat);lon=dat(ilon);dep=dat(idep);slip=dat(islip)
                rake=dat(irake);strike=dat(istrike);dip=dat(idip)
                slip=slip*1d-2        !cm to meter
                slip_strike=slip*dcos(rake*pi/180d0)
                slip_downdip=-slip*dsin(rake*pi/180d0)
                slip_op=0d0
      
            if (.not.center)then !If output the top left of each patch
            ! compute the top left coordinates from center of each patch
              call center2topleft(length,wid,strike,dip,lat,lon,dep)
            endif

                depmo=dep-oceanthi
                if(depmo.le.0d0)then
                  write(*,*)'Warning: Depth above ocean bottom, &
     & force it to be below!'
                  depmo=0.001d0
                endif
      
                write(2,fmt='(I8,7(1x,E23.15),2(1x,I3),1x,f3.1)')  &
                    i,lat,lon,depmo,length,wid,strike,dip,1,1,0d0
                write(2,fmt='(8x,2(1x,E23.15),3(1x,E23.15))') &
                   length/2d0,wid/2d0,slip_strike,slip_downdip,slip_op
              enddo  !do j=1,nx*ny
              deallocate(strtype,dat)
            enddo !do ifs=1,nfs
            close(2)
            read(1,fmt='(a200)', iostat=err)string
            if(.not.err.lt.0)write(*,*)'Supposed to be the end of file, &
     & but not!'
       !    recl is expressed in bytes (characters) for formatted form
       !         is expressed in 4-bytes unit for unformatted
            open(2,file=ofile,form='formatted',access='direct',recl=8)
!           do j=1,8
!             read(2,fmt='(a8)',rec=j,iostat=err)strtmp
!             write(*,*)err,strtmp
!           enddo
            write(2,fmt='(I8)',rec=1)cnt
            close(2)
          endif !elseif (faultflag.eq.1)
          close(1)

        elseif(trim(fmttype).eq.'2')then
          !write(*,*)'Ex: caltech/caltech.out caltech/caltech_pscmp.out 1.0'
          open(1,file=ifile)
          read(1,fmt='(a200)')string
          ntype=9
          allocate(strtype(ntype),dat(ntype))
          read(string(2:200),fmt=*)(strtype(k),k=1,ntype)
!c----- set order of the record types -----
!% lat1     lon1      wid(km) length(km)  dep(km)  dip(degree) strike(degree) rake(degree) slip(cm)
          ilat=99
          ilon=99
          iwid=99;ilen=99;
          idep=99
          islip=99
          irake=99
          istrike=99
          idip=99
          do k=1,ntype
            strtmp=strtype(k)
            if(strtmp(1:3).eq."lat") ilat = k
            if(strtmp(1:3).eq."lon") ilon = k
            if(strtmp(1:3).eq."wid") iwid = k
            if(strtmp(1:3).eq."len") ilen = k
            if(strtmp(1:3).eq."dep") idep = k
            if(strtmp(1:4).eq."slip") islip = k
            if(strtmp(1:4).eq."rake") irake = k
            if(strtmp(1:6).eq."strike") istrike = k
            if(strtmp(1:3).eq."dip") idip = k
          enddo
          open(2,file=ofile)
          write(2,fmt='(I8)')0
          i=0
          do while(.true.)
            read(1,fmt='(a200)', iostat=err)string
            if(err.lt.0)exit
            i=i+1
!           write(*,*)string
!           read(string,fmt=*)lat,lon,wid,length,dep,dip,strike,rake,slip
            read(string,fmt=*)(dat(k),k=1,ntype)
            lat=dat(ilat);lon=dat(ilon);dep=dat(idep);slip=dat(islip)
            rake=dat(irake);strike=dat(istrike);dip=dat(idip)
            length=dat(ilen);wid=dat(iwid)
            slip=slip*1d-2        !cm to meter
            slip_strike=slip*dcos(rake*pi/180d0)
            slip_downdip=-slip*dsin(rake*pi/180d0)
            slip_op=0d0
  
            if (.not.center)then !If output the top left of each patch
            ! compute the top left coordinates from center of each patch
              call center2topleft(length,wid,strike,dip,lat,lon,dep)
            endif

            depmo=dep-oceanthi
            if(depmo.le.0d0)then
              write(*,*)'Warning: Depth above ocean bottom, &
     & force it to be below!'
              depmo=0.001d0
            endif
            write(2,fmt='(I8,7(1x,E23.15),2(1x,I3),1x,f3.1)')  &
              i,lat,lon,depmo,length,wid,strike,dip,1,1,0d0
            write(2,fmt='(8x,2(1x,E23.15),3(1x,E23.15))') &
              length/2d0,wid/2d0,slip_strike,slip_downdip,slip_op
          enddo
          close(1)
          close(2)
          nt=i
          open(2,file=ofile,form='formatted',access='direct',recl=8)
          write(2,fmt='(I8)',rec=1)nt
          close(2)

        elseif(trim(fmttype).eq.'3')then
!       write(*,*)'Ex: ../Sumatra/SM04Caltech.txt ../Sumatra/SM04Caltech_pscmp.out 1.5'

        length=20d0;wid=16d0 !Length width see readdcl
        write(*,*)'Plese check length width=',length,wid

        open(1,file=ifile)
        do i=1,3
          read(1,fmt='(a200)')string
        enddo
        i=0
        do while(.true.)
          read(1,fmt='(a200)', iostat=err)string
          if(err.lt.0)exit
          i=i+1
        end do
        nt=i
        rewind(1)
        do i=1,3
          read(1,fmt='(a200)')string
        enddo

        open(2,file=ofile)
        write(2,fmt='(I8)')nt
        i=0
        do while(.true.)
          read(1,fmt='(a200)', iostat=err)string
          if(err.lt.0)exit
          i=i+1
          read(string,fmt=*)lon,lat,slip,rake,strike,dip,dep
          slip=slip*1d-2        !cm to meter
          slip_strike=slip*dcos(rake*pi/180d0)
          slip_downdip=-slip*dsin(rake*pi/180d0)
          slip_op=0d0

          if (.not.center)then !If output the top left of each patch
          ! compute the top left coordinates from center of each patch
            call center2topleft(length,wid,strike,dip,lat,lon,dep)
          endif

          depmo=dep-oceanthi
          if(depmo.le.0d0)then
            write(*,*)'Warning: Depth above ocean bottom,&
     & force it to be below!'
            depmo=0.001d0
          endif
          write(2,fmt='(I8,7(1x,E23.15),2(1x,I3),1x,f3.1)')  &
              i,lat,lon,depmo,length,wid,strike,dip,1,1,0d0
          write(2,fmt='(8x,2(1x,E23.15),3(1x,E23.15))') &
                length/2d0,wid/2d0,slip_strike,slip_downdip,slip_op
        enddo
        close(1)
        close(2)

        elseif(trim(fmttype).eq.'4')then
!         write(*,*)'Ozawa_coseismic_fault_para.dat Ozawa_coseismic_fault_pscmp.dat 3.9433'
          open(1,file=ifile)
          do i=1,1
            read(1,fmt='(a200)')string
          enddo
          i=0
          do while(.true.)
            read(1,fmt='(a200)', iostat=err)string
            if(string(2:6).eq.'99999')exit
            i=i+1
          end do
          nt=i
          write(*,*)'Number of lines:',nt
          rewind(1)
          do i=1,1
            read(1,fmt='(a200)')string
          enddo
  
          open(2,file=ofile)
          write(2,fmt='(I8)')nt
          i=0
          do while(.true.)
            read(1,fmt='(a200)', iostat=err)string
            if(string(2:6).eq.'99999')exit
        !   if(err.lt.0)exit
            i=i+1
            read(string,fmt=*)lat,lon,slip,rake,dep,strike,dip,length,wid
            slip_strike=slip*dcos(rake*pi/180d0)
            slip_downdip=-slip*dsin(rake*pi/180d0)
            slip_op=0d0
  
            if (.not.center)then !If output the top left of each patch
            ! compute the top left coordinates from center of each patch
              call center2topleft(length,wid,strike,dip,lat,lon,dep)
            endif

            depmo=dep-oceanthi
            if(depmo.le.0d0)then
              write(*,*)'Warning: Depth above ocean bottom, &
     & force it to be below!'
              depmo=0.001d0
            endif
            write(2,fmt='(I8,7(1x,E23.15),2(1x,I3),1x,f3.1)')  &
                i,lat,lon,depmo,length,wid,strike,dip,1,1,0d0
            write(2,fmt='(8x,2(1x,E23.15),3(1x,E23.15))') &
                length/2d0,wid/2d0,slip_strike,slip_downdip,slip_op
          enddo
          close(1)
          close(2)

        elseif(trim(fmttype).eq.'GCMTndk')then
          !refers to SimonsSoftware/cmtsol.m 
          length=1d0; wid=1d0 
!         length=1d-3; wid=1d-3 
!         cmtcode='M201103110546A';
          cmtcode=trim(ifile)
          fname='/0/home/dai.56/share/Osuworks/UltraSH/Slepian/IFILES/CM&
     &T/jan76_feb10.ndk' 
          ! file name of the .ndk data
          call cmtsol(cmtcode,fname,lat,lon,dep,M0,strike,dip,rake); 

          mu=30e9; ! in Pa; 1Pa = 1 N/m^2; Value in pscmp08a_code_in/pscdisc.f
          slip=M0/(length*1e3*wid*1e3*mu)  !in meter

          slip_strike=slip*dcos(rake*pi/180d0)
          slip_downdip=-slip*dsin(rake*pi/180d0)
          slip_op=0d0
      
          if (.not.center)then !If output the top left of each patch
          ! compute the top left coordinates from center of each patch
            call center2topleft(length,wid,strike,dip,lat,lon,dep)
          endif

          depmo=dep-oceanthi
          if(depmo.le.0d0)then
            write(*,*)'Warning: Depth above ocean bottom, &
     & force it to be below!'
                depmo=0.001d0
          endif
  
          cnt=1;i=1
          open(2,file=ofile)
          write(2,fmt='(I8)')cnt
          write(2,fmt='(I8,7(1x,E23.15),2(1x,I3),1x,f3.1)')  &
            i,lat,lon,depmo,length,wid,strike,dip,1,1,0d0
          write(2,fmt='(8x,2(1x,E23.15),3(1x,E23.15))') &
             length/2d0,wid/2d0,slip_strike,slip_downdip,slip_op
          close(2)

        elseif(trim(fmttype).eq.'SIV')then !Jan 23, 2015
!       write(*,*)'Sumatra12/s2012SUMATR01YUEx.fsp fault_pscmp.out 3.9433'
          open(1,file=ifile)
          faultflag=0
          do while(.true.)
            read(1,fmt='(a200)', iostat=err)string
            if(err.lt.0)exit   !exit do while
            do i=1,80
              ie=i+17-1
              if (string(i:i+4-1).eq.'Mech')then
                do j=1,80
                  if (string(j:j+4-1).eq.'STRK')then
                    strp=string(j+4+2:j+4+2+16) 
                    read(strp,*)strike
                  elseif(string(j:j+3-1).eq.'DIP')then
                    strp=string(j+3+2:j+3+2+16)
                    read(strp,*)dip
                  endif
                enddo !do j
                exit  !do i
              elseif (string(i:i+10-1).eq.'Invs :  Dx')then
                do j=1,80
                  if (string(j:j+2-1).eq.'Dx')then
                    strp=string(j+2+3:j+2+3+16) 
                    read(strp,*)length
                  elseif(string(j:j+2-1).eq.'Dz')then
                    strp=string(j+2+3:j+2+3+16)
                    read(strp,*)wid
                  endif
                enddo !do j
                exit  !do i
              elseif (string(i:ie).eq.'of fault segments')then
                strp=string(27:30) !=     1
                read(strp,*)nfs
                faultflag=1
                exit  !do i
              endif
            enddo  !do i
!           if (faultflag.eq.1) exit    !exit do while
            if(string(3:25).eq.'SOURCE MODEL PARAMETERS')then
              exit  !exit do while
            endif
          end do  !do while

          if(faultflag.eq.0)then
            write(*,*)'Error: found no fault_segments!'
          elseif (faultflag.eq.1)then
            !Read the fault_segment
            cnt=0
            i=0
            open(2,file=ofile)
            write(2,fmt='(I8)')0
            do ifs=1,nfs
              !Read header of each segment
              do while(.true.)
                read(1,fmt='(a200)', iostat=err)string
                if(err.lt.0)exit   !exit do while
                if (string(1:10).eq.'% SEGMENT ')then
                  do j=1,80
                    if (string(j:j+6-1).eq.'STRIKE')then
                      strp=string(j+6+2:j+6+2+16) 
                      read(strp,*)strike
                    elseif(string(j:j+3-1).eq.'DIP')then
                      strp=string(j+3+2:j+3+2+16)
                      read(strp,*)dip
                    endif
                  enddo !do j
                  cycle ! skip the rest,go to do while
                endif
                !scrutinize the string
                do k=1,80
                  if (string(k:k+4-1).eq.'Dx =')then
                    do j=1,80
                      if (string(j:j+2-1).eq.'Dx')then
                        strp=string(j+2+2:j+2+2+16) 
                        read(strp,*)length
                      elseif(string(j:j+2-1).eq.'Dz')then
                        strp=string(j+2+2:j+2+2+16)
                        read(strp,*)wid
                      endif
                    enddo !do j
                    exit  !do k
                  elseif (string(k:k+5-1).eq.'Nsbfs')then
                    strp=string(k+5+2:k+5+2+4) !=     1
                    read(strp,*)nsbfs
                    exit  !do k
                  endif
                enddo  !do k
                if(string(1:8).eq.'%    LAT')then
                  exit  !exit do while
                endif
              enddo  !do while
              read(1,fmt=*)
              !End of header of each segment

              cnt=cnt+nsbfs
              ntype=7
              allocate(strtype(ntype),dat(ntype))
              read(string(2:200),fmt=*)(strtype(k),k=1,ntype)
!c----- set order of the record types -----
!%    LAT     LON     X==EW     Y==NS     Z     SLIP     RAKE
              ilat=99
              ilon=99
              idep=99
              islip=99
              irake=99
              istrike=99
              idip=99
              do k=1,ntype
                strtmp=strtype(k)
                if(strtmp(1:3).eq."LAT") ilat = k
                if(strtmp(1:3).eq."LON") ilon = k
                if(strtmp(1:1).eq."Z") idep = k
                if(strtmp(1:4).eq."SLIP") islip = k
                if(strtmp(1:4).eq."RAKE") irake = k
                if(strtmp(1:6).eq."strike") istrike = k
                if(strtmp(1:3).eq."dip") idip = k
              enddo
              
              do j=1,nsbfs
                read(1,fmt='(a200)', iostat=err)string
                if(err.lt.0)exit
                i=i+1
!               write(*,*)string
!               read(string,fmt=*)lat,lon,X,Y,dep,slip,rake
                read(string,fmt=*)(dat(k),k=1,ntype)
                lat=dat(ilat);lon=dat(ilon);dep=dat(idep);slip=dat(islip)
                rake=dat(irake); !strike=dat(istrike);dip=dat(idip)
!               slip=slip*1d-2        !cm to meter
                slip_strike=slip*dcos(rake*pi/180d0)
                slip_downdip=-slip*dsin(rake*pi/180d0)
                slip_op=0d0
      
                if (.not.center)then !If output the top left of each patch
            ! compute the top left coordinates from top-center of each patch
!             call center2topleft(length,wid,strike,dip,lat,lon,dep)
                  x=-length/2d0; y=0d0;
                  call rf2any(x,y,strike,dip,lat,lon,dep,1)
                endif

                depmo=dep-oceanthi
                if(depmo.le.0d0)then
                  write(*,*)'Warning: Depth above ocean bottom, &
     & force it to be below!'
                  depmo=0.001d0
                endif
      
                write(2,fmt='(I8,7(1x,E23.15),2(1x,I3),1x,f3.1)')  &
                    i,lat,lon,depmo,length,wid,strike,dip,1,1,0d0
                write(2,fmt='(8x,2(1x,E23.15),3(1x,E23.15))') &
                   length/2d0,wid/2d0,slip_strike,slip_downdip,slip_op
              enddo  !do j=1,nsbfs
              deallocate(strtype,dat)
            enddo !do ifs=1,nfs
            close(2)
       !    recl is expressed in bytes (characters) for formatted form
       !         is expressed in 4-bytes unit for unformatted
            open(2,file=ofile,form='formatted',access='direct',recl=8)
            write(2,fmt='(I8)',rec=1)cnt
            close(2)
          endif !elseif (faultflag.eq.1)
          close(1)

        elseif(trim(fmttype).eq.'Vigny')then
!       write(*,*)'Ex: Chile/vigny_etal_2011_model.dat fault_pscmp.out 1.3267 '
          open(1,file=ifile)
          nh=0
          do while(.true.)
            read(1,fmt='(a200)', iostat=err)string
            if(err.lt.0)exit
            if (string(1:1).eq.'#')then
                nh=nh+1
            else
                exit
            endif
          end do
          write(*,*)'Lines of header:',nh
          rewind(1)
          do i=1,nh
            read(1,fmt='(a200)')string
          enddo
          i=0
          do while(.true.)
            read(1,fmt='(a200)', iostat=err)string
            if(err.lt.0)exit
            i=i+1
          end do
          nt=i
          write(*,*)'Number of lines:',nt
          rewind(1)
          do i=1,nh
            read(1,fmt='(a200)')string
          enddo

          i=0
          open(2,file=ofile)
          write(2,fmt='(I8)')nt
          do while(.true.)
            read(1,fmt='(a200)', iostat=err)string
            if(err.lt.0)exit
            i=i+1
            read(string,fmt=*)lon,lat,dep,length,wid,dip,strike, &
     & slip_strike,slip_downdip,slip_op
            !slip in cm
            slip_strike=slip_strike*1d-2 !cm to meter
            slip_downdip=-slip_downdip*1d-2 !dip direction given updip
            slip_op=slip_op*1d-2

            if (.not.center)then !If output the top left of each patch
            ! compute the top left coordinates from center of each patch
!             call center2topleft(length,wid,strike,dip,lat,lon,dep)
              x=0d0; y=-wid;
              call rf2any(x,y,strike,dip,lat,lon,dep,1)
            endif

            depmo=dep-oceanthi
            if(depmo.le.0d0)then
              write(*,*)'Warning: Depth above ocean bottom, &
     &  force it to be below!'
              depmo=0.001d0
            endif

            write(2,fmt='(I8,7(1x,E23.15),2(1x,I3),1x,f3.1)')  &
                i,lat,lon,depmo,length,wid,strike,dip,1,1,0d0
            write(2,fmt='(8x,2(1x,E23.15),3(1x,E23.15))') &
                length/2d0,wid/2d0,slip_strike,slip_downdip,slip_op
          enddo
          close(1)
          close(2)
        endif

        stop
        end
        include 'rf2any.f90'

        subroutine center2topleft(length,wid,strike,dip,lat,lon,dep)
        !output the top left of each patch
        !compute the top left coordinates from center of each patch
        !Input: length, width, of this patch, in km
        !       dip angle of this patch, in degree
        !       latitude, longitude of the center of this patch,in degrees
        !       depth of the center, in degree
        !Output: latitude, longitude, depth of the top left point, in
        !        degree and km
        !       
        double precision length,wid,strike,dip,lat,lon,dep
        double precision Cl,Cw,dipangle,c,alfa,colat,z1
        double precision pi,d2r,aa
        aa=6371.d0
        pi=datan(1d0)*4d0
        d2r=PI/180D0
        
        Cl=length
        Cw=wid
        dipangle=dip*d2r; z1=strike*d2r
        colat=pi/2d0-lat*d2r
        c=dsqrt(((0.5d0)*Cw*dcos(dipangle))**2+((0.5d0)*Cl)**2)  !km
        alfa=datan((0.5d0)*Cw*dcos(dipangle)/((0.5d0)*Cl))      !radian
        colat=colat+c*dcos(z1+alfa)/aa      !radian
        lon=lon-c*dsin(z1+alfa)/(aa*dsin(colat))/d2r  !degree
        dep=dep-(0.5d0)*Cw*dsin(dipangle)       !km
        lat=90d0-colat/d2r  !degree

        return
        end


        subroutine cmtsol(cmtcode,fname,lat,lon,dep,M0,strike,dip,rake)
        implicit none
        character cmtcode*80, fname*160,string*200,MomentType*5
        integer lencmt,faultflag,err,iexp
        double precision HalfDuration,time,stdt
        double precision lat,stdlat,lon,stdlon,dep,stddep
        double precision Mt(6),std(6),M0,strike,dip,rake
        integer iStrike, idip, irake,i
        logical alive
        
        inquire(file=fname,exist=alive)
        if(.not.alive)then
          write(*,*)fname," doesn't exist."
          stop
        endif

        lencmt=len(trim(cmtcode))

        faultflag=0
        open(1,file=fname)
        do while(.true.)
          read(1,fmt='(a200)', iostat=err)string
          if(err.lt.0)exit
          if (string(1:lencmt).eq.trim(cmtcode))then
            faultflag=1
            exit
          endif
        end do

        if(faultflag.eq.0)then
          write(*,*)'Error: found no matched CMT data!'
        elseif (faultflag.eq.1)then
!Second line
!M201103110546A   B:  0    0   0 S:  0    0   0 M:100  271 300 CMT: 1 TRIHD: 70.0
          read(string,fmt='(69x,a5,1x,f5.1)')MomentType,HalfDuration

!Third line
!CENTROID:     69.8 0.2  37.52 0.01  143.05 0.02  20.0  0.0 FIX  Q-20110311234057
          read(1,fmt='(a200)', iostat=err)string
          if(string(1:8).ne.'CENTROID')write(*,*)'Format read wrong'
          read(string,fmt='(9x,f9.1,f4.1,f7.2,f5.2,f8.2,f5.2,f6.1,f5.1)')time,stdt,lat,stdlat,lon,stdlon,dep,stddep

!Fourth line
!29  1.730 0.006 -0.281 0.005 -1.450 0.005  2.120 0.068  4.550 0.065 -0.657 0.004
          read(1,fmt='(a200)', iostat=err)string
!         read(string,fmt='(I2,6(f7.3,f6.3))')iexp,Mrr,std(1),Mtt,std(2),Mpp,std(3),Mrt,std(4),Mrp,std(5),Mtp,std(6)
          read(string,fmt='(I2,6(f7.3,f6.3))')iexp,(Mt(i),std(i),i=1,6)
          !Mrr, Mtt, Mpp, Mrt, Mrp, Mtp
          !where r is up, t is south, and p is east. See Aki and
          !Richards for conversions to other coordinate systems.
          !Mrr*10**iexp [dyne-cm]; 1 dyne*cm =1e-7 N*m

!Fifth line
!V10   5.305 55 295   0.014  0 205  -5.319 35 115   5.312 203 10   88  25 80   90
          read(1,fmt='(a200)', iostat=err)string
          read(string,fmt='(48x,f8.3,I4,I3,I5)')M0,iStrike, idip, irake 
          Mt=Mt*10d0**dble(iexp-7d0); M0=M0*10d0**dble(iexp-7d0)  !in N*m
          strike=istrike; dip=idip; rake=irake

        endif !elseif (faultflag.eq.1)
        close(1)

        return
        end
