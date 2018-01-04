        subroutine rf2any(x,y,strike,dip,lat,lon,dep,method)
        !Compute the long lat of ANY point, given the lat lon of
        !ReFerence point, distance and azimuth
        !Input: 
        !     x, y: the coordinate the unknown point refers to the
        !           referece point, in km. x is the strike direction,
        !           and y is the dip direction, origin is at the
        !           reference point
        !     strike: strike angle of this fault plane, in degree
        !     dip:  dip angle of the fault plane, in degree
        !     lat,lon: latitude, longitude of the reference point, in degrees
        !     dep:  the depth of the reference point, in km 
        !     method: 1, the algorithm in Sun Wenke's code GravityInt.f; 
        !             2, the algorithm in Wang Rongjiang's code PSCMP
        !Output: latitude, longitude, depth of the unknown point, in
        !        degree and km
        !       
        implicit none
        integer method
        double precision x,y,strike,dip,lat,lon,dep
        double precision Cl,Cw,dipangle,c,alfa,colat,z1
        double precision pi,d2r,aa
        double precision ys,azim,sma,smb,smc,bga,bgc
        aa=6371.d0
        pi=datan(1d0)*4d0
        d2r=PI/180D0
        
        dipangle=dip*d2r; z1=strike*d2r
        colat=pi/2d0-lat*d2r
        ys=y*dcos(dipangle)  !y projection to Earth's surface
        c=dsqrt(ys*ys+x*x)  !surface distance, km
        alfa=datan2(ys,x)      !radian, [-pi, pi]
        if(alfa.lt.0)alfa=alfa+pi*2d0  ![0,2pi]
        azim=alfa+z1  ! azimuth, radian
        !Compute the long lat of ANY point, given the lat lon of
        !ReFerence point, distance and azimuth
        if (method.eq.1)then
          colat=colat-c*dcos(azim)/aa      !radian
          lon=lon+c*dsin(azim)/(aa*dsin(colat))/d2r  !degree
        elseif(method.eq.2)then
!c             spherical triangle:
!c             A = pole, B = source position, C = reference position
!c
              sma=c/aa
              smb=colat
              bgc=azim
!             smc=dacos(dcos(sma)*dcos(smb)                             &   
!    &           +dsin(sma)*dsin(smb)*dcos(bgc))
          smc=dacos(dcos(sma)*dcos(smb)+dsin(sma)*dsin(smb)*dcos(bgc))
              bga=dasin(dsin(sma)*dsin(bgc)/dsin(smc))
!c
!c             geographic coordinate of the point source
!c
              colat=smc
              lon=lon+bga/d2r  !degree
        endif
        lat=90d0-colat/d2r  !degree
        dep=dep+y*dsin(dipangle)       !km
  
        return
        end
