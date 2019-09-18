c
c       This file tests various options for
c       computing M2M translations
c
c        Currently, three options are being tested
c
c         1)  Form the matrix, and use blas
c         2)  Use vectorized mpmp translation operator
c             to compute several translations at the same time
c         3) Serially compute all translations
c
      implicit real *8 (a-h,o-z)
      integer nterms,nsize
      complex *16, allocatable :: mpmpmat(:,:)
      complex *16, allocatable :: mpolein(:,:,:)
      complex *16, allocatable :: mpolein0(:,:,:)
      complex *16, allocatable :: mpoleout0(:,:,:)
      complex *16, allocatable :: mpoleout(:,:,:)
      complex *16, allocatable :: mpolein2(:,:)
      complex *16, allocatable :: mpoleout2(:,:)
      complex *16 zk,ima,cra1,cra2,cra3,alpha,beta
      real *8 rscale0, rscale1, center0(3),center1(3)
      real *8 xnodes(10000),wts(10000)
      integer, parameter :: nd = 100000
      character *1 transa,transb 
      data ima/(0.0d0,1.0d0)/

      call prini(6,13)

c
c
c       generate data for testing
c

      nterms = 39
      nsize = (nterms+1)**2

      zk = 1.1d0 + ima*0.1d0
      radius = 0.9d0

      rscale0 = hkrand(0)
      rscale1 = hkrand(0)
      

      center0(1) = 0
      center0(2) = 0
      center0(3) = 0

      center1(1) = hkrand(0)
      center1(2) = hkrand(0)
      center1(3) = hkrand(0)


      allocate(mpolein(nd,0:nterms,-nterms:nterms))
      allocate(mpolein0(0:nterms,-nterms:nterms,nd))
      allocate(mpoleout(nd,0:nterms,-nterms:nterms))
      allocate(mpoleout0(0:nterms,-nterms:nterms,nd))
      allocate(mpolein2(nsize,nd))
      allocate(mpoleout2(nsize,nd))

      do idim=1,nd
        ii = 1 
        do i=0,nterms
          do j=-i,i
            mpolein(idim,i,j) = hkrand(0) + ima*hkrand(0)
            mpolein2(ii,idim) = mpolein(idim,i,j)
            mpolein0(i,j,idim) = mpolein(idim,i,j)

            mpoleout2(ii,idim) = 0
            mpoleout(idim,i,j) = 0
            mpoleout0(i,j,idim) = 0
            ii = ii +1
          enddo
        enddo
      enddo

c
c
c       Strategy 1:
c

      allocate(mpmpmat(nsize,nsize))
c
c        get translation operator matrix
c
      call getmpmpmat(nterms,nsize,zk,rscale0,rscale1,center0,center1,
     1       radius,mpmpmat)

c
c        test translation operator matrix on a random rhs
c
      call testmpmpmat(nterms,nsize,zk,rscale0,rscale1,center0,center1,
     1       radius,mpmpmat,err)
      print *, "error in mpmp mat=", err


      transa = 'n'
      transb = 'n'
      alpha = 1
      beta = 0
      print *, "Starting strategy 1"
     
      
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call zgemm(transa,transb,nsize,nd,nsize,alpha,mpmpmat,nsize,
     1       mpolein2,nsize,beta,mpoleout2,nsize)
      call cpu_time(t2) 
C$      t2 = omp_get_wtime()     
      
      print *, "Time for blas=",(t2-t1)/nd

c
c
c       strategy 2: vectorized mpmp operators
c
      ifinit2 = 1
      nquad2 = nterms*2.5
      call legewhts(nquad2,xnodes,wts,ifinit2)

      t1 = 0
      t2 = 0

      print *, "Starting strategy 2"
      
      call cpu_time(t1) 
C$      t1 = omp_get_wtime()      
      call h3dmpmp(nd,zk,rscale0,center0,mpolein,nterms,
     1      rscale1,center1,mpoleout,nterms,radius,xnodes,
     2      wts,nquad2)
      call cpu_time(t2) 
C$      t2 = omp_get_wtime()     
 
      print *, "Time mpmp vec=",(t2-t1)/nd

c
c
c      Strategy 3: apply mpmp operator 1 expansion at a time
c
      t1 = 0
      t2 = 0
      print *, "Starting strategy 3"
      call cpu_time(t1)
C$       t1 = omp_get_wtime()      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim)     
      do idim=1,nd
      call h3dmpmp(1,zk,rscale0,center0,mpolein0(0,-nterms,idim),
     1   nterms,rscale1,center1,mpoleout0(0,-nterms,idim),nterms,
     2   radius,xnodes,wts,nquad2)
      enddo
C$OMP END PARALLEL DO
      call cpu_time(t2)
C$       t2 = omp_get_wtime()      
      print *, "time mpmp=",(t2-t1)/nd


c
c       Compute sum of quantities and write first three output
c       vectors to file
c
c       This is done to ensure that all quantities are
c       computed post compiler optimizations
c
      cra1 = 0
      cra2 = 0
      cra3 = 0

      err0 = 0
      err1 = 0

 1000 format(3(2x,e11.5))
      do idim=1,nd
        ii = 1 
        do i=0,nterms
          do j=-i,i
            cra1 = cra1+ mpoleout2(ii,idim) 
            cra2 = cra2 + mpoleout(idim,i,j) 
            cra3 = cra3 + mpoleout0(i,j,idim)
            if(idim.le.3)
     1       write(33,1000) real(mpoleout2(ii,idim)),
     1         real(mpoleout(idim,i,j)),real(mpoleout0(i,j,idim))
            ii = ii +1
          enddo
        enddo
      enddo
      print *, cra1
      print *, cra2
      print *, cra3
      




      stop
      end
c------------------------------------------------------------------

      subroutine getmpmpmat(nterms,nsize,zk,rscale0,rscale1,center0,
     1             center1,radius,xmat)
      implicit real *8 (a-h,o-z)
      complex *16 xmat(nsize,nsize),zk
      real *8 center0(3),center1(3),rscale0,rscale1
      real *8 xnodes(10000),wts(10000)
      complex *16, allocatable :: mpolein(:,:)
      complex *16, allocatable :: mpoleout(:,:)

      allocate(mpolein(0:nterms,-nterms:nterms))
      allocate(mpoleout(0:nterms,-nterms:nterms))

      ifinit2 = 1
      nquad2 = nterms*2.5
      call legewhts(nquad2,xnodes,wts,ifinit2)

      icol = 1
      do ii=0,nterms
        do jj=-ii,ii

           do i=0,nterms
             do j=-i,i
               mpolein(i,j) = 0
               mpoleout(i,j) = 0
             enddo
           enddo
           mpolein(ii,jj) = 1
           ier = 0
           call h3dmpmp(1,zk,rscale0,center0,mpolein,nterms,
     1           rscale1,center1,mpoleout,nterms,radius,xnodes,
     2           wts,nquad2)

           irow = 1
           do i=0,nterms
             do j=-i,i
               xmat(irow,icol) = mpoleout(i,j)
               irow = irow+1
             enddo
           enddo
           icol = icol + 1
        enddo
      enddo

      return
      end

c-------------------------------------------------------------
      subroutine testmpmpmat(nterms,nsize,zk,rscale0,rscale1,
     1     center0,center1,radius,xmat,err)
      implicit real *8 (a-h,o-z)
      complex *16 xmat(nsize,nsize),zk
      real *8 center0(3),center1(3),rscale0,rscale1
      real *8 xnodes(10000),wts(10000)
      complex *16, allocatable :: mpolein(:,:),mpolein2(:)
      complex *16, allocatable :: mpoleout(:,:),mpoleout2(:)
      complex *16 ima,alpha,beta
      data ima/(0.0d0,1.0d0)/
      character *1 transa

      allocate(mpolein(0:nterms,-nterms:nterms))
      allocate(mpoleout(0:nterms,-nterms:nterms))
      allocate(mpolein2(nsize),mpoleout2(nsize))

      ifinit2 = 1
      nquad2 = nterms*2.5
      call legewhts(nquad2,xnodes,wts,ifinit2)

      ii = 1
      do i=0,nterms
        do j=-i,i
          mpolein(i,j) = hkrand(0) + ima*hkrand(0)
          mpolein2(ii) = mpolein(i,j)
          ii = ii + 1
        enddo
      enddo

      transa = 'n'
      incx = 1
      incy = 1

      alpha = 1
      beta = 0

      call zgemv(transa,nsize,nsize,alpha,xmat,nsize,mpolein2,incx,
     1   beta,mpoleout2,incy)
      

      call h3dmpmp(1,zk,rscale0,center0,mpolein,nterms,
     1       rscale1,center1,mpoleout,nterms,radius,xnodes,
     2       wts,nquad2,ier)

      ii = 1
      erra = 0
      do i=0,nterms
        do j=-i,i
          erra = abs(mpoleout(i,j)-mpoleout2(ii))**2

cc          call prin2('mpoleout=*',mpoleout(i,j),2)
cc          call prin2('mpoleout2=*',mpoleout2(ii),2)
          ii = ii +1
        enddo
      enddo

      erra = sqrt(erra)

      return
      end
c---------------------------------------------------
      
