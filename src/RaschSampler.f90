subroutine sampler(n,k,inputmat,tfixed,burn_in,n_eff,step,seed,outputvec,ier)

      ! sample binary matrices with given marginals.
      ! input: n = number of rows (integer*4)
      !        k = number of columns (integer*4)
      !        inputmat: input binary matrix (n*k) (integer*4)
      !        tfixed: main diagonal is fixed if true; has only effect if n.eq.k
      !        step: if the matrix #i is an effective matrix, then the next effective matrix is #(i+step)
      !        burn_in: number of burn_in matrices in units of step
      !        n_eff: number of effective matrices, to be written in the output vector
      ! I/O    seed (integer*4). Seed of the random generator
      !          if seed <> 0: seed is unaltered and |seed| is used as the seed of the random generator.
      !          if seed.eq.0: a seed is generated from the system clock, and its value is returned.
      !          ATTENTION: currently seed.eq.0 is deactivated and only seed<>0 is allowed
      !                    (see lines 103-113) (rh 2006-10-25)
      ! output:outputvec (integer*4 vector): n_eff binary matrices, stored bitwise in the following way:
      !           if(k<=32) one row of a matrix is stored in one position (number) of four bytes, the first element
      !                     is bit 0, the second is bit 1, etc.
      !           if(32<k<=64) two positions are used per row: the frst 32 in position 1, the remainder in position 2,
      !                     again starting with element 33 in bit 0, element 34 in bit 1, etc.
      !           not used bits are set to zero.
      !        ier: output error code
      !             ier = 0: O.K.
      !                   1: n > nmax = 1024 = 2**10   !!!!!! changed to 2**12 in 0.8-3
      !                   2: k > kmax = 64 = 2**6      !!!!!! changed to 2**7  in 0.8-3
      !                   4: n_eff > n_effmax = 8191 = 2**13 - 1
      !                   8: burn_in < 0
      !                  16: step <= 0
      !                1-31: sums of the foregoing codes
      !                  32: input matrix contains values other than one or zero
      !                  64: the input matrix has a Guttman form
      ! if tfixed and n.eq.k, the main diagonal is condidered as fixed.
      ! if tfixed and n!=k, the case is treated as a rectangular matrix without constraints
      ! the Markov chain transition matrix used is Q**step

      integer(kind=4), dimension(n*(k+31)/32*(n_eff+1)),intent(out)::outputvec
      integer(kind=4), intent(out)     :: ier
      integer(kind=4), intent(in)      :: n,k,burn_in,n_eff,step
      integer(kind=4), dimension(n,k),intent(in) :: inputmat
      integer(kind=4)                 :: seed
      logical(kind=4), intent(in)     :: tfixed

      character(len=10)               :: timevec

      !integer(kind=4),parameter       :: nmax=1024,kmax=64,n_effmax=8191  !!!!!! kmax changed to 2**7 nmax changed to 2**12
      integer(kind=4),parameter       :: nmax=4096,kmax=256,n_effmax=8191
      integer(kind=4), allocatable    :: a(:),b(:),aold(:),bold(:),a_kol(:),b_kol(:),iwork(:)
      integer(kind=4)                 :: i,j,m,kk2,kk3,it,krand,k2,k3,nhex,k2old,k3old,nhexold
      integer(kind=4)                 :: x1,x2  ! x1 and x2 are reserved for the random generator
      integer(kind=4)                 :: words_per_row, offset
      integer(kind=4),allocatable     :: hexmat(:,:),hexmatold(:,:)

      real(kind=4)                    :: tijd

      logical(kind=1),parameter       :: t=.true.,f=.false.
      logical(kind=1),allocatable     :: t_in(:,:)
      logical(kind=1),dimension(3,3)  :: hexa,hexb
      logical(kind=1),allocatable     :: twa(:),twb(:),tw(:),tng(:),tngold(:),col1(:),col2(:) !tng = non-guttman pair
      logical(kind=1)                 :: t_eff,tfixnow

      data hexa/f,f,t,t,f,f,f,t,f/,hexb/f,t,f,f,f,t,t,f,f/

      !check error codes 1, 2 4, 8 and 16
      ier=0
      if(n.le.0 .or. n.gt.nmax)ier=ier+1
      if(k.le.0 .or. k.gt.kmax)ier=ier+2
      if(n_eff.le.0 .or. n_eff.gt.n_effmax)ier=ier+4
      if(burn_in.lt.0)ier=ier+8
      if(step.le.0)ier=ier+16
      if(ier.ne.0)return

      ! allocate the necessary arrays

      kk2=k*(k-1)/2
      allocate(t_in(n,k))
      allocate(a(kk2),b(kk2),aold(kk2),bold(kk2),a_kol(kk2),b_kol(kk2))
      allocate(twa(n),twb(n),tw(n),tng(kk2),tngold(kk2),col1(n),col2(n))
      allocate (iwork(n))

      tfixnow=tfixed
      if(n.ne.k)tfixnow=.false.
      if(tfixnow)then
        kk3=kk2*(k-2)/3
        allocate(hexmat(3,kk3),hexmatold(3,kk3))
      endif

      ! check error code 32
      ier=count(inputmat.gt.1)+count(inputmat.lt.0)
      if(ier.ne.0)then
        ier=32
        return
      endif

      ! copy input matrix to t_in
      !!!t_in=inputmat
      !!!replaced by
      t_in=btest(inputmat,0)
      if(tfixnow) then
        forall (i=1:n) t_in(i,i)=f
      endif
!      !select seed for random number generation
!      if(seed.eq.0)then
!        call date_and_time(TIME=timevec)
!        read(timevec,'(f10.3)')tijd
!        x1=tijd*1000.
!        x1=0.
!        x1=x1+536870911 ! =x1 + 2**29 - 1
!        seed=x1
!      else
!        x1=abs(seed)
!      endif
      x1=abs(seed) ! added from upper else clause (to be removed if random seed enabled)
      x2=x1
      call rand1(x1)
      krand=1+mod(k,25)  ! KRAND selects a random generator (used in RAND_INTEGER42)

      !fill the arrays a_kol, b_kol, a, b, and b_pairs for the input matrices t_in
      !determine the weight k2 (= #neighbour column pairs of the input matrix)
      it=0
      do i=2,k
        do j=1,i-1
          it=it+1
          a_kol(it)=i
          b_kol(it)=j
          call findab(t_in(1:n,i),t_in(1:n,j),i,j,a(it),b(it))
          tng(it)=a(it)*b(it).gt.0
        end do
      end do
      k2=count(tng(1:kk2))
      k3=k2
      if(tfixnow)then
        call hexagon
        if(nhex.gt.0)k3=k2+1
      endif
      ! check on the Guttman condition (error code 64)
      if(k3.eq.0)then
        ier=64
        return
      endif

      ! pack the input matrix and put it in the outputvector
      offset=0 ! offset is the number of words defined until now
      words_per_row=(k+31)/32
      call pack_matrix(outputvec) ! offset is updated within the subroutine

      ! compute the total number of matrices to be computed
      n_tot=burn_in+n_eff
      ! do the sampling
      do i=1,n_tot
        t_eff=i.gt.burn_in
        do j=1,step  ! generate step matrices before computing the statistic
          call rand_integer42(it,k3,x1,x2,krand)
          !{*** from here to ***} only applies to alternating hexagons
          if(k3.gt.k2.and.it.eq.k3)then  ! there are restrictions on the main diagonal
                                         ! and there exists at least one alternating hexagon (k3>k2)
                                         ! an alternating hexagon has to be changed (it=k3)
            ! first save necessary elements for Metropolis-Hastings
            k2old=k2
            k3old=k3
            aold=a
            bold=b
            nhexold=nhex
            hexmatold(:,1:nhex)=hexmat(:,1:nhex)
            tngold=tng
            ! make a new matrix by switching one alternating hexagon
            call rand_integer42(it,nhex,x1,x2,krand)
            call make_matrix3(it)
            call update_pairs3(hexmat(1,it),hexmat(2,it),hexmat(3,it))
            call update_pairs3(hexmat(2,it),hexmat(1,it),hexmat(3,it))
            call update_pairs3(hexmat(3,it),hexmat(1,it),hexmat(2,it))
            call hexagon
            k2=count(tng(1:kk2))
            k3=k2+1 ! there is at least one alternating hexagon, vz, the one that just changed into its complement
            if(k3old.lt.k3) then ! check if process possibly remains in the same state (Metropolis-Hastings)
              call rand_integer42(m,k3,x1,x2,krand)
              if(m.gt.k3old)then ! process remains in the same state
                call make_matrix3(it) ! this restores the matrix
                k3=k3old
                k2=k2old
                a=aold
                b=bold
                nhex=nhexold
                hexmat(:,1:nhex)=hexmatold(:,1:nhex)
                tng=tngold
              endif
            endif
            cycle
          endif
          !***}
          call pick_a_pair(it)
          col1=t_in(1:n,a_kol(it))
          col2=t_in(1:n,b_kol(it))
          k2old=k2
          k3old=k3
          aold=a
          bold=b
          tngold=tng
          if(tfixnow)then
            nhexold=nhex
            hexmatold(:,1:nhex)=hexmat(:,1:nhex)
          endif
          call make_matrix(it)
          call update_pairs(a_kol(it),b_kol(it))
          call update_pairs(b_kol(it),a_kol(it))
          k2=count(tng(1:kk2))
          k3=k2
          if(tfixnow)then
            call hexagon
            if(nhex.gt.0)k3=k3+1
          endif
          if(k3old.lt.k3)then  ! apply Metropolis-Hastings
            call rand_integer42(m,k3,x1,x2,krand)
            if(m.gt.k3old)then ! process remains in the current state
              t_in(1:n,a_kol(it))=col1(1:n)
              t_in(1:n,b_kol(it))=col2(1:n)
              k2=k2old
              k3=k3old
              a=aold
              b=bold
              tng=tngold
              if(tfixnow)then
                nhex=nhexold
                hexmat(:,1:nhex)=hexmatold(:,1:nhex)
              endif
            endif
          endif
        end do
        ! if this is an 'effective' sample, pack it and store it in outputvec
        if(t_eff)call pack_matrix(outputvec)
      end do ! here ends the sampling procedure

      deallocate(a,b,aold,bold,a_kol,b_kol,iwork,twa,twb,tw,tng,tngold,col1,col2,t_in)
      if(tfixnow) deallocate(hexmat,hexmatold)


      contains

      subroutine pack_matrix(vec)
        integer(kind=4) vec(*)
        integer(kind=4) :: i,j,ib,ie,it,iw
        do i=1,n
          ib=1
          do iw=1,words_per_row
            offset=offset+1
            vec(offset)=0
            ie=min(ib+31,k)
            it=-1
            do j=ib,ie
              it=it+1
              if(.not.t_in(i,j))cycle
              vec(offset)=ibset(vec(offset),it)
            end do
            ib=ie+1
          end do
        end do
      end subroutine pack_matrix

      subroutine findab(ta,tb,i,j,a,b)
 !!!!!  logical(kind=1):: ta(n),tb(n)
        logical(kind=1):: ta(n),tb(n),test(n)
        integer(kind=4)::a,b,i,j
        tw=(ta.neqv.tb)
        if(tfixnow)then
          tw(i)=.false.
          tw(j)=.false.
        endif
             test=ta.and.tw
             a=count(test)
             test=tb.and.tw
             b=count(test)
 !!!!!       a=count(ta.and.tw)
 !!!!!       b=count(tb.and.tw)
      end subroutine findab


      subroutine pick_a_pair(it)
        integer(kind=4),intent(out)::it
        integer(kind=4) ::i,m
        call rand_integer42(it,k2,x1,x2,krand)
        m=count(tng(1:it))
        if(m.eq.it)return
        do i=it+1,kk2
          if(.not.tng(i))cycle
          m=m+1
          if(m.eq.it)then
            it=i
            return
          endif
        end do
      end subroutine pick_a_pair

      subroutine make_matrix(it)
        integer(kind=4),intent(in)::it
        integer(kind=4)           ::m,i,j,ii,jj
!!!!!       logical(kind=1),allocatable     :: test(:)
        logical(kind=1) :: test(n)

        ii=a_kol(it)
        jj=b_kol(it)
        if(a(it)*b(it).eq.1)then ! columns ii and jj contain a single tetrad.
                                 ! no sampling is necessary: the tetrad is complemented
          j=0
          do i=1,n
            if(tfixnow)then
              if(i.eq.ii)cycle
              if(i.eq.jj)cycle
            endif
            if(t_in(i,ii).eqv.t_in(i,jj))cycle
            j=j+1
            t_in(i,ii)=t_in(i,jj)
            t_in(i,jj)=.not.t_in(i,ii)
            if(j.eq.2)return
          end do
        endif
        do                     ! a random binomial operation is applied
          ! copy the two selected colums in the logical vectors twa and twb
          twa(1:n)=t_in(1:n,ii)
          twb(1:n)=t_in(1:n,jj)
          m=a(it)+b(it)
          call combine(m,a(it),iwork,tw) ! generate the random combination of a(it) objets out of m
          ! insert the combination into the vectors twa and twb
          j=0
          do i=1,n
            if(tfixnow)then
              if(i.eq.ii)cycle
              if(i.eq.jj)cycle
            endif
            if(twa(i).eqv.twb(i))cycle
            j=j+1
            if(tw(j))then
              twa(i)=.true.
              twb(i)=.false.
            else
              twa(i)=.false.
              twb(i)=.true.
            endif
            if(j.eq.m)exit
          end do
          ! check whether matrix has changed
 !!!!!
 !!!!!    test=twa
          test=twa(1:n).eqv.t_in(1:n,ii)
          m=count(test)
 !!!!!         m=count(twa(1:n).eqv.t_in(1:n,ii))
          if(m.ne.n)exit ! the matrix has changed
        end do  ! the matrix has not changed; a new combination is tried.
        ! the changes are inserted in the matrix t_in
        t_in(1:n,ii)=twa(1:n)
        t_in(1:n,jj)=twb(1:n)

      end subroutine make_matrix

      subroutine make_matrix3(it)
        integer(kind=4), intent(in)::it
        integer(kind=4)            ::i,j,ii,jj
        do i=1,2
          ii=hexmat(i,it)
          do j=i+1,3
            jj=hexmat(j,it)
            t_in(ii,jj)=.not.t_in(ii,jj)
            t_in(jj,ii)=.not.t_in(jj,ii)
          end do
        end do
      end subroutine make_matrix3

      subroutine combine(n,k,ix,tx)
        ! generate a random combination of k objects out of n
        ! the result is stored in the logical n-vector tx
        ! ix is a working array
        integer(kind=4) ::n,k,kk,ii,iu,nnu
        integer(kind=4) ::ix(n)
        logical(kind=1) ::tx(n)

        ix(1:n)=(/(ii,ii=1,n)/)
        tx(1:n)=.false.
        nnu=n
        kk=min(k,n-k)
        do ii=1,kk
          call rand_integer42(iu,nnu,x1,x2,krand)
          tx(ix(iu))=.true.
          if(iu.lt.nnu)ix(iu:nnu-1)=ix(iu+1:nnu)
          nnu=nnu-1
        end do
        if(kk.lt.k)tx=.not.tx
      end subroutine combine

      subroutine update_pairs(i,j)
        integer(kind=4),intent(in) :: i,j
        integer(kind=4):: jt,m
        do m=1,i-1
          if(m.eq.j)cycle
          jt=(i-1)*(i-2)/2+m
          call findab(t_in(1:n,i),t_in(1:n,m),i,m,a(jt),b(jt))
          tng(jt)=a(jt)*b(jt).gt.0
        end do
        do m=i+1,k
          if(m.eq.j)cycle
          jt=(m-1)*(m-2)/2+i
          call findab(t_in(1:n,m),t_in(1:n,i),m,i,a(jt),b(jt))
          tng(jt)=a(jt)*b(jt).gt.0
        end do
      end subroutine update_pairs

      subroutine update_pairs3(i,j,l)
        integer(kind=4),intent(in) :: i,j,l
        integer(kind=4):: jt,m
        do m=1,i-1
          if(m.eq.j)cycle
          if(m.eq.l)cycle
          jt=(i-1)*(i-2)/2+m
          call findab(t_in(1:n,i),t_in(1:n,m),i,m,a(jt),b(jt))
          tng(jt)=a(jt)*b(jt).gt.0
        end do
        do m=i+1,k
          if(m.eq.j)cycle
          if(m.eq.l)cycle
          jt=(m-1)*(m-2)/2+i
          call findab(t_in(1:n,m),t_in(1:n,i),m,i,a(jt),b(jt))
          tng(jt)=a(jt)*b(jt).gt.0
        end do
      end subroutine update_pairs3

      subroutine hexagon
        logical(kind=1),dimension(3,3):: c
        integer(kind=4),dimension(3)  :: v
        integer(kind=4)               :: i,j,m
        nhex=0
        do i=1,n-2
          v(1)=i
          do j=i+1,n-1
            if(t_in(i,j).eqv.t_in(j,i))cycle
            v(2)=j
            do m=j+1,n
              v(3)=m
              c=t_in(v,v)
              if(all(c.eqv.hexa).or.all(c.eqv.hexb))then
                nhex=nhex+1
                hexmat(1:3,nhex)=v
              endif
            end do
          end do
        end do
      end subroutine hexagon

      subroutine rand1(x)
        ! see Ripley, B.D., Stochastic Simulation. New-York:Wiley, 1987, pp. 37-39. (generator 4)
        integer(kind=4) x,p,q,r,a,b,c

        data a,b,c,p/127773,16807,2836,2147483647/     ! generator 4
        q=x/a
        r=mod(x,a)
        x=b*r-c*q
        if(x.lt.0)x=x+p
      end subroutine rand1

      subroutine rand_integer42(iu,a,x1,x2,k)
        ! draw a uniformly distributed integer from {1, 2, ..., A}
        ! the rejection rate is (P - BOUND)/P
        ! where P equals 2**31-2 (= the number of different values RAND1 can take)
        ! and BOUND equals A*(P/A)
        ! example: for A = 1000, the rejection rate is 3.008E-7
        ! Notice that zero as result of the draw leads to rejection (see !***)
        ! The routine calls RAND2 with the K-th set of coefficients

        integer(kind=4) iu,a,x1,x2,bound,k
        integer(kind=4),parameter :: p = 2147483646
                            ! P is the number of values X can take: 2**31 - 2,
                            ! because 0 is excluded
        bound=(p/a)*a
        do
          call rand2(x1,x2,k)
          if(x2.eq.0)cycle
          if(x2.le.bound)exit        !*** LE not LT
        end do
        iu=1+mod(x2,a)

      end subroutine rand_integer42

      subroutine rand2(x1,x2,k)

        ! Random number generators MRG (multiple recursive generator)
        ! Lih-Yuan Deng and Dennis K.J. Lin (2000).Random Number Generation for the New Century,
        ! The American Statistician, vol 54, no. 2, pp. 145-150
        ! To compute the formulae, the method in the referenced article is used. B(K) is the table as published
        ! A(K)=P/B(K) and C(K)=P-A(K)*B(K), P = 2**31 -1

        integer (kind=4)::k
        integer (kind=4)::x1,x2,p,q,r,y
        integer (kind=4), dimension(25)::b,c,a

        data b/26403,33236,36673,40851,43693,27149,33986,36848,40961, &
               44314,29812,34601,37097,42174,44530,30229,36098,37877, &
               42457,45670,31332,36181,39613,43199,46338/
        data a/81334,64613,58557,52568,49149,79099,63187,58279,52427, &
               48460,72034,62064,57888,50919,48225,71040,59490,56696, &
               50580,47021,68539,59353,54211,49711,46343/
        data c/22045, 5979,22786,28279,16390,24896,10265,19055,21300, &
               27207, 6039, 7183,12511,25741,24397,15487,13627, 9255, &
                8587,34577,19699,32754,23304,18158,41713/
        data p/2147483647/

        q=x1/a(k)
        r=mod(x1,a(k))
        y=b(k)*r-c(k)*q
        if(y.ge.x2-p)then
          y=y-x2
        else
          y=y+(p-x2)
        endif
        if(y.lt.0)y=y+p
        x1=x2
        x2=y

      end subroutine rand2

      end subroutine sampler



      subroutine unpack(vec,words_per_row,t_out,n,k)


      integer(kind=4) offset,words_per_row,n,k
      integer(kind=4),dimension(n*words_per_row)::vec
      integer(kind=4) i,j,it,ib,ie,ioff
      ! matrix t_in is not needed
      !logical(kind=1),dimension(n,k) :: t_in
      integer(kind=4),dimension(n,k) :: t_out

      t_out=0 ! intialize t_out
      ioff=0
      do i=1,n
        ib=1

        do iw=1,words_per_row
          ioff=ioff+1
          ie=min(ib+31,k)
          it=-1
          do j=ib,ie
            it=it+1
            !t_in(i,j)=btest(vec(ioff),it)
            !replace the preceding statement by
            if(btest(vec(ioff),it)) t_out(i,j)=1
          end do
          ib=ie+1
        end do

      end do

      end subroutine unpack
