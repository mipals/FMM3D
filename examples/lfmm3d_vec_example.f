      implicit none
      integer ns,nt,nd
      double precision, allocatable :: source(:,:),targ(:,:)
      double precision, allocatable :: charge(:,:),dipvec(:,:,:)
      double precision, allocatable :: pot(:,:),grad(:,:,:)
      double precision, allocatable :: potex(:,:),gradex(:,:,:)
      double precision, allocatable :: pottarg(:,:),gradtarg(:,:,:)
      double precision, allocatable :: pottargex(:,:),gradtargex(:,:,:)

      double precision eps
      integer i,j,k,idim,ntest
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

      ns = 100000
      nt = 19

      ntest = 10
      
      nd = 1

      allocate(source(3,ns))
      allocate(targ(3,nt))
      allocate(charge(nd,ns))
      allocate(dipvec(nd,3,ns))

      allocate(pot(nd,ns))
      allocate(potex(nd,ntest))

      allocate(grad(nd,3,ns))
      allocate(gradex(nd,3,ntest))
      
      
      allocate(pottarg(nd,nt))
      allocate(pottargex(nd,ntest))
      
      allocate(gradtarg(nd,3,nt))
      allocate(gradtargex(nd,3,ntest))

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
        source(2,i) = hkrand(0)**2
        source(3,i) = hkrand(0)**2

        do idim=1,nd
          charge(idim,i) = hkrand(0) 
          dipvec(idim,1,i) = hkrand(0) 
          dipvec(idim,2,i) = hkrand(0) 
          dipvec(idim,3,i) = hkrand(0) 
          pot(idim,i) = 0
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



       call lfmm3d_s_cd_g_vec(nd,eps,ns,source,charge,
     1      dipvec,pot,grad)

       do i=1,ntest
         do idim=1,nd
           potex(1,i) = 0
           gradex(idim,1,i) = 0
           gradex(idim,2,i) = 0
           gradex(idim,3,i) = 0

           pottargex(idim,i) = 0
           gradtargex(idim,1,i) = 0
           gradtargex(idim,2,i) = 0
           gradtargex(idim,3,i) = 0
         enddo

       enddo

       call l3ddirectcdg(nd,source,charge,
     1      dipvec,ns,source,ntest,potex,gradex,thresh)

       call prin2("potential at sources=*",pot,6)
       call prin2("potential at sources=*",potex,6)


       erra = 0
       ra = 0
       do i=1,ntest
         ra = ra + potex(1,i)**2
         erra = erra + (potex(1,i)-pot(1,i))**2
       enddo

       erra = sqrt(erra/ra)
       call prin2("error pot=*",erra,1)

      erra = 0
      ra = 0
      do i=1,ntest
        do j=1,3
          erra = erra + (gradex(1,j,i)-grad(1,j,i))**2
          ra = ra + (gradex(1,j,i))**2
        enddo
      enddo

      erra = sqrt(erra/ra)

      call prin2('error grad=*',erra,1)

      call prin2_long('grad=*',grad,3*ntest)
      call prin2_long('gradex=*',gradex,3*ntest)




      stop
      end
c----------------------------------------------------------
c
cc
c
c
