      implicit none
      integer ns,nt,nd
      double precision, allocatable :: source(:,:),targ(:,:)
      double precision, allocatable :: charge(:,:),dipvec(:,:,:)
      double precision, allocatable :: pot(:,:)
      double precision, allocatable :: potex(:,:)
      double precision, allocatable :: pottarg(:,:)

      double precision eps
      integer i,j,k,idim
      double precision hkrand,thresh,erra,ra
      


c
cc      initialize printing routine
c
      call prini(6,13)
      write(*,*)
      write(*,*)
      write(*,*) "================================="
      call prin2("This code is an example fortran driver*",i,0)
      call prin2("On output, the code prints sample pot,pottarg*",i,0)
      write(*,*)
      write(*,*)

      ns = 13
      nt = 19
      
      nd = 1

      allocate(source(3,ns))
      allocate(targ(3,nt))
      allocate(charge(nd,ns))
      allocate(dipvec(nd,3,ns))
      allocate(pot(nd,ns))
      allocate(potex(nd,ns))
      allocate(pottarg(nd,nt))

      thresh =1.0d-16

      eps = 0.5d-9

      write(*,*) "=========================================="

c
c   
c       example demonstrating use of 
c        source to source+targ, charges+dipoles, pot
c



c
cc      generate sources uniformly in the unit cube 
c
c
      do i=1,ns
        source(1,i) = hkrand(0)**2
        source(2,i) = hkrand(0)**2*0
        source(3,i) = hkrand(0)**2

        do idim=1,nd
          charge(idim,i) = hkrand(0) 
          dipvec(idim,1,i) = hkrand(0) 
          dipvec(idim,2,i) = hkrand(0) 
          dipvec(idim,3,i) = hkrand(0) 
          pot(idim,i) = 0
          potex(idim,i) = 0
        enddo
      enddo

      do i=1,nt
        targ(1,i) = hkrand(0)**3
        targ(2,i) = hkrand(0)**3
        targ(3,i) = hkrand(0)**3

        do idim=1,nd
          pottarg(idim,i) = 0
        enddo
      enddo


C       call lfmm3d_st_cd_p_vec(nd,eps,ns,source,charge,
C     1      dipvec,pot,nt,targ,pottarg)

       call lfmm3d_s_cd_p_vec(nd,eps,ns,source,charge,
     1      dipvec,pot)

       call l3ddirectcdp(nd,source,charge,
     1      dipvec,ns,source,ns,potex,thresh)

cc       call prin2_long("potential at sources=*",pot,20)
cc       call prin2_long("potential at sources=*",potex,20)
C       print *, pot
C       print *, potex


       erra = 0
       ra = 0
       do i=1,ns
         ra = ra + potex(1,i)**2
         erra = erra + (potex(1,i)-pot(1,i))**2

       enddo

       erra = sqrt(erra/ra)
       call prin2("error=*",erra,1)


      stop
      end
c----------------------------------------------------------
c
cc
c
c
