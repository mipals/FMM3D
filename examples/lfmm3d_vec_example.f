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

      ns = 200000
      nt = 0

      ntest = 1000
      
      nd = 8

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
          grad(idim,1,i) = 0
          grad(idim,2,i) = 0
          grad(idim,3,i) = 0
        enddo
      enddo

      if(.false.) then
         source(1,1) = 0.00d0
         source(2,1) = 0.00d0
         source(3,1) = 0.00d0

         source(1,2) = 1.00d0
         source(2,2) = 0.00d0
         source(3,2) = 0.00d0

         source(1,3) = 0.00d0
         source(2,3) = 1.00d0
         source(3,3) = 0.00d0

         source(1,4) = 1.00d0
         source(2,4) = 1.00d0
         source(3,4) = 0.00d0

         source(1,5) = 0.00d0
         source(2,5) = 0.00d0
         source(3,5) = 1.00d0

         source(1,6) = 1.00d0
         source(2,6) = 0.00d0
         source(3,6) = 1.00d0

         source(1,7) = 0.00d0
         source(2,7) = 1.00d0
         source(3,7) = 1.00d0

         source(1,8) = 1.00d0
         source(2,8) = 1.00d0
         source(3,8) = 1.00d0

cccccccccccccccccccccc
         source(1,9) = 0.51d0
         source(2,9) = 0.01d0
         source(3,9) = 0.01d0

         source(1,10) = 0.51d0
         source(2,10) = 0.49d0
         source(3,10) = 0.01d0

         source(1,11) = 0.99d0
         source(2,11) = 0.49d0
         source(3,11) = 0.01d0

         source(1,12) = 0.51d0
         source(2,12) = 0.01d0
         source(3,12) = 0.49d0

         source(1,13) = 0.99d0
         source(2,13) = 0.01d0
         source(3,13) = 0.49d0

         source(1,14) = 0.51d0
         source(2,14) = 0.49d0
         source(3,14) = 0.49d0

         source(1,15) = 0.99d0
         source(2,15) = 0.49d0
         source(3,15) = 0.49d0

C         source(1,16) = 0.74d0
C         source(2,16) = 0.01d0
C         source(3,16) = 0.25d0
      endif




      do i=1,nt
        targ(1,i) = hkrand(0)**3
        targ(2,i) = hkrand(0)**3
        targ(3,i) = hkrand(0)**3

        do idim=1,nd
          pottarg(idim,i) = 0
          gradtarg(idim,1,i) = 0
          gradtarg(idim,2,i) = 0
          gradtarg(idim,3,i) = 0
        enddo
      enddo



       call lfmm3d_s_c_p_vec(nd,eps,ns,source,charge,
     1      pot,nt)

       do i=1,ntest
         do idim=1,nd
           potex(idim,i) = 0
           gradex(idim,1,i) = 0
           gradex(idim,2,i) = 0
           gradex(idim,3,i) = 0

           pottargex(idim,i) = 0
           gradtargex(idim,1,i) = 0
           gradtargex(idim,2,i) = 0
           gradtargex(idim,3,i) = 0
         enddo

       enddo

       call l3ddirectcp(nd,source,charge,
     1      ns,source,ntest,potex,thresh)


C       call l3ddirectcdg(nd,source,charge,
C     1      dipvec,ns,targ,ntest,pottargex,gradtargex,thresh)

       call prin2("potential at sources=*",pot,6)
       call prin2("potential at sources=*",potex,6)


       erra = 0
       ra = 0
       do i=1,ntest
         do idim=1,nd
           ra = ra + potex(idim,i)**2
           erra = erra + (potex(idim,i)-pot(idim,i))**2
         enddo
       enddo

       erra = sqrt(erra/ra)
       call prin2("error pot src=*",erra,1)

      erra = 0
      ra = 0
      do i=1,ntest
        do j=1,3
          do idim=1,nd
            erra = erra + (gradex(idim,j,i)-grad(idim,j,i))**2
            ra = ra + (gradex(idim,j,i))**2
          enddo
        enddo
      enddo

      erra = sqrt(erra/ra)

      call prin2('error grad src=*',erra,1)



c       erra = 0
c       ra = 0
c       do i=1,ntest
c         do idim=1,nd
c           ra = ra + pottargex(idim,i)**2
c           erra = erra + (pottargex(idim,i)-pottarg(idim,i))**2
c         enddo
c       enddo
c
c       erra = sqrt(erra/ra)
c       call prin2("error pot targ=*",erra,1)
c
c      erra = 0
c      ra = 0
c      do i=1,ntest
c        do j=1,3
c          do idim=1,nd
c            erra = erra + (gradtargex(idim,j,i)-gradtarg(idim,j,i))**2
c            ra = ra + (gradtargex(idim,j,i))**2
c          enddo
c        enddo
c      enddo
c
c      erra = sqrt(erra/ra)
c
c      call prin2('error grad targ=*',erra,1)




      stop
      end
c----------------------------------------------------------
c
cc
c
c
