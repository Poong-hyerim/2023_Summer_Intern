program FDM_2DF

        implicit none

        integer :: nx, nz, nt, nf, ix, iz, it, ifreq, ixz, G, k, node, lindex, lvalue, nnt
        real :: dx, dz, dt, df, fmax, tmax, c, pi, alpha, ndt !c=min.velocity
        complex :: ww

        real, allocatable, dimension(:) :: source, wk, nwk
        integer, allocatable, dimension(:) :: iwk, niwk
        complex, allocatable, dimension(:) :: csource, u, cf
        complex, allocatable, dimension(:,:) :: green

        integer, allocatable, dimension(:) :: indexx, col
        complex, allocatable, dimension(:) :: w, x, value
        integer keep(20), info(40), icntl(20)
        real(8) cntl(10), rinfo(20)          !real*8 or real(8)  


        pi=acos(-1.)
        nx=200
        nz=100
        dx=0.005
        dz=0.005
        tmax=1.
        G=16
        c=1.5
        fmax=c/G/dx
        alpha=alog(100.)/tmax

        dt=0.1/fmax
        nt=nint(tmax/dt)+1
        ndt=0.01/fmax
        nnt=nint(tmax/ndt)+1

        df=1./tmax
        nf=nint(fmax/df)+1

        allocate(source(nt), csource(nt), u(nnt))
        allocate(iwk(6*nt+150), wk(6*nt+150))
        allocate(niwk(6*nnt+150), nwk(6*nnt+150))
        allocate(green(nx*nz,nf))

        lindex=200*nx*nz
        lvalue=200*nx*nz
        allocate(indexx(lindex), value(lvalue), cf(nx*nz), w(4*nx*nz), x(nx*nz),col(lindex))



        open(10,file='2D_f.dat',access='direct',form='unformatted',recl=4*nnt)

        call fdgaus(source,fmax,dt,nt)
        csource(:)=cmplx(source(:),0.)

        call fftcc(csource,nt,iwk,wk)
        green(:,:)=0.

        do ifreq=2,nf
                ww=cmplx(2.*pi*(ifreq-1)*df,alpha)
                node=0
                indexx(:)=0
                col(:)=0
                do ix=1, nx
                        do iz=1,nz
                                if(ix.gt.1) then
                                        node=node+1
                                        indexx(node)=(ix-1)*nz+iz
                                        col(node)=(ix-2)*nz+iz
                                        value(node)=cmplx(-1/dx**2,0.)
                                endif
                                if(iz.gt.1) then
                                        node=node+1
                                        indexx(node)=(ix-1)*nz+iz
                                        col(node)=(ix-1)*nz+iz-1
                                        value(node)=cmplx(-1/dz**2,0.)
                                endif
                                node=node+1
                                indexx(node)=(ix-1)*nz+iz
                                col(node)=(ix-1)*nz+iz
                                value(node)=-ww**2/c**2+cmplx(2/dx**2+2/dz**2,0.)
                                if(iz.lt.nz) then
                                        node=node+1
                                        indexx(node)=(ix-1)*nz+iz
                                        col(node)=(ix-1)*nz+iz+1
                                        value(node)=cmplx(-1/dz**2,0.)
                                endif
                                if(ix.lt.nx) then
                                        node=node+1
                                        indexx(node)=(ix-1)*nz+iz
                                        col(node)=ix*nz+iz
                                        value(node)=cmplx(-1/dx**2,0.)
                                endif
                        enddo
                enddo
                indexx(node+1:2*node)=col(1:node)
                cf(:)=cmplx(0.,0.)
                cf((nx/2-1)*nz+2)=cmplx(1.,0.)    !x=nx/2, z=2에서 소스
                call umc21i(keep,cntl,icntl)
                        icntl(3)=2
                        icntl(6)=1
                call umc2fa(nx*nz,node,0,.false.,lvalue,lindex,value,indexx,keep,cntl,icntl,info,rinfo)
                
                if(lindex.lt.info(19).or.lvalue.lt.info(21)) then
                        write(*,*) 'dimension error'
                        write(*,*) 'lindex=',lindex,',lvalue=',lvalue
                        write(*,*) 'info(19)=',info(19),',info(21)=',info(21)
                        stop
                endif
                call umc2so(nx*nz,0,.false.,lvalue,lindex,value,indexx,keep,cf,x,w,cntl,icntl,info,rinfo)

                do ixz=1,nx*nz
                        green(ixz,ifreq)=x(ixz)
                enddo
        enddo


        !write
        do ix=1,nx
                u(:)=0
                do ifreq=2,nf
                        u(ifreq)=green((ix-1)*nz+2,ifreq)*csource(ifreq)
                enddo

                do ifreq=2,nf
                        u(nnt-ifreq+2)=conjg(u(ifreq))
                enddo

                do it=1,nnt
                        u(it)=conjg(u(it))
                enddo

                call fftcc(u,nnt,niwk,nwk)

                write(10,rec=ix) (real(u(it))*exp(alpha*it*ndt)/nnt,it=1,nnt)
        enddo

        print*, nnt
        close(10)


        deallocate(source,csource,u,iwk,wk,niwk,nwk,green,indexx,value,cf,w,x,col)

stop
end program



!
!
    subroutine fdgaus(w,cutoff,dt,nt)
    dimension w(1)
    phi=4*atan(1.)
!****I changed a, a=phi*cutoff**2
    a=phi*(5.*cutoff/8.)**2
    amp=sqrt(a/phi)
    do 10 i=1,nt
        t=(i-1)*dt
        arg=-a*t**2
        if(arg.lt.-32.) arg=-32.
10  w(i)=amp*exp(arg)
    do 20 i=1,nt
        if(w(i).lt.0.001*w(1)) then
        icut=i
        t0=(icut-1)*dt
        go to 30
        endif
20    continue
30    continue
    do 40 i=1,nt
        t=(i-1)*dt
        t=t-t0
        arg=-a*t**2
        if(arg.lt.-32.) arg=-32.
40    w(i)=-2.*sqrt(a)*a*t*exp(arg)/sqrt(phi)
        smax=0.
    do i=1,nt
        if(abs(w(i)).gt.smax) smax=abs(w(i)) 
    enddo
    do i=1,nt
            w(i)=w(i)/smax
    enddo
    return
    end

!
!*********************************************************************
!   imsl routine name   - fftcc
!*********************************************************************
!-----------------------------------------------------------------------
!
!   latest revision     - january 1, 1978
!
!   purpose             - compute the fast fourier transform of a
!                           complex valued sequence
!
!   usage               - call fftcc (a,n,iwk,wk)
!
!   arguments    a      - complex vector of length n. on input a
!                           contains the complex valued sequence to be
!                           transformed. on output a is replaced by the
!                           fourier transform.
!                n      - input number of data points to be
!                           transformed. n may be any positive
!                           integer.
!                iwk    - integer work vector of length 6*n+150.
!                           (see programming notes for further details)
!                wk     - real work vector of length 6*n+150.
!                           (see programming notes for further details)
!
!   reqd. imsl routines - none required
!
!   notation            - information on special notation and
!                           conventions is available in the manual
!                           introduction or through imsl routine uhelp
!
!   remarks  1.  fftcc computes the fourier transform, x, according
!                to the following formula;
!
!                  x(k+1) = sum from j = 0 to n-1 of
!                           a(j+1)*cexp((0.0,(2.0*pi*j*k)/n))
!                  for k=0,1,...,n-1 and pi=3.1415...
!
!                note that x overwrites a on output.
!            2.  fftcc can be used to compute
!
!                  x(k+1) = (1/n)*sum from j = 0 to n-1 of
!                           a(j+1)*cexp((0.0,(-2.0*pi*j*k)/n))
!                  for k=0,1,...,n-1 and pi=3.1415...
!
!                by performing the following steps;
!
!                     do 10 i=1,n
!                        a(i) = conjg(a(i))
!                  10 continue
!                     call fftcc (a,n,iwk,wk)
!                     do 20 i=1,n
!                        a(i) = conjg(a(i))/n
!                  20 continue
!
!-----------------------------------------------------------------------
!
    subroutine fftcc (a,n,iwk,wk)
!                                  specifications for arguments
    integer            n,iwk(1)
    real               wk(1)
    complex            a(n)
!                                  specifications for local variables
    integer            i,iam,iap,ibm,ibp,ic,icc,icf,ick,id,idm1,ii,  &
                    ija,ikb,ikt,ill,im,ird,isf,isk,isp,iss,ita,itb,&
                    j,ja,jf,jj,jk,k,k0,k1,k2,k3,ka,kb,kd2,kf,kh,kn,&
                    kt,ktp,l,l1,m,mm,mm1,mp
    real               cm,sm,c1,c2,c3,s1,s2,s3,c30,rad,a0,a1,a4,b4,  &
                    a2,a3,b0,b1,b2,b3,zero,half,one,two,z0(2),     &
                    z1(2),z2(2),z3(2),z4(2)
    complex            za0,za1,za2,za3,za4,ak2
    equivalence        (za0,z0(1)),(za1,z1(1)),(za2,z2(1)),          &
                    (za3,z3(1)),(a0,z0(1)),(b0,z0(2)),(a1,z1(1)),  &
                    (b1,z1(2)),(a2,z2(1)),(b2,z2(2)),(a3,z3(1)),   &
                    (b3,z3(2)),(za4,z4(1)),(z4(1),a4),(z4(2),b4)
    data               rad/6.283185/,                                &
                    c30/.8660254/
    data               zero,half,one,two/0.0,0.5,1.0,2.0/
!                                  first executable statement
    if (n .eq. 1) go to 9005
    k = n
    m = 0
    j = 2
    jj = 4
    jf = 0
!                                  determine the square factors of n
    iwk(1) = 1
5 i = k/jj
    if (i*jj .ne. k) go to 10
    m = m+1
    iwk(m+1) = j
    k = i
    go to 5
10 j = j + 2
    if (j .eq. 4) j = 3
    jj = j * j
    if (jj .le. k) go to 5
    kt = m
!                                  determine the remaining factors of n
    j = 2
15 i = k / j
    if (i*j .ne. k) go to 20
    m = m + 1
    iwk(m+1) = j
    k = i
    go to 15
20 j = j + 1
    if (j .eq. 3) go to 15
    j = j + 1
    if(j.le.k) go to 15
    k = iwk(m+1)
    if (iwk(kt+1) .gt. iwk(m+1)) k = iwk(kt+1)
    if(kt.le.0) go to 30
    ktp = kt + 2
    do 25  i = 1,kt
        j = ktp - i
        m = m+1
        iwk(m+1) = iwk(j)
25 continue
30 mp = m+1
    ic = mp+1
    id = ic+mp
    ill = id+mp
    ird = ill+mp+1
    icc = ird+mp
    iss = icc+mp
    ick = iss+mp
    isk = ick+k
    icf = isk+k
    isf = icf+k
    iap = isf+k
    kd2 = (k-1) / 2 + 1
    ibp = iap + kd2
    iam = ibp + kd2
    ibm = iam + kd2
    mm1 = m-1
    i = 1
35 l = mp - i
    j = ic - i
    iwk(ill+l) = 0
    if ((iwk(j-1) + iwk(j)) .eq. 4) iwk(ill+l) = 1
    if (iwk(ill+l) .eq. 0) go to 40
    i = i + 1
    l = l - 1
    iwk(ill+l) = 0
40 i = i + 1
    if(i.le.mm1) go to 35
    iwk(ill+1) = 0
    iwk(ill+mp) = 0
    iwk(ic) = 1
    iwk(id) = n
    do 45  j = 1,m
        k = iwk(j+1)
        iwk(ic+j) = iwk(ic+j-1) * k
        iwk(id+j) = iwk(id+j-1) / k
        wk(ird+j) = rad/iwk(ic+j)
        c1 = rad/k
        if (k .le. 2) go to 45
        wk(icc+j) = cos(c1)
        wk(iss+j) = sin(c1)
45 continue
    mm = m
    if (iwk(ill+m) .eq. 1) mm = m - 1
    if (mm .le. 1) go to 50
    sm = iwk(ic+mm-2) * wk(ird+m)
    cm = cos(sm)
    sm = sin(sm)
50 kb = 0
    kn = n
    jj = 0
    i = 1
    c1 = one
    s1 = zero
    l1 = 1
55 if (iwk(ill+i+1) .eq. 1) go to 60
    kf = iwk(i+1)
    go to 65
60 kf = 4
    i = i+1
65 isp = iwk(id+i)
    if (l1 .eq. 1) go to 70
    s1 = jj * wk(ird+i)
    c1 = cos(s1)
    s1 = sin(s1)
!                                  factors of 2, 3, and 4 are
!                                  handled separately.
70 if (kf .gt. 4) go to 140
    go to (75,75,90,115), kf
75 k0 = kb + isp
    k2 = k0 + isp
    if (l1 .eq. 1) go to 85
80 k0 = k0 - 1
    if (k0 .lt. kb) go to 190
    k2 = k2 - 1
    za4 = a(k2+1)
    a0 = a4*c1-b4*s1
    b0 = a4*s1+b4*c1
    a(k2+1) = a(k0+1)-za0
    a(k0+1) = a(k0+1)+za0
    go to 80
85 k0 = k0 - 1
    if (k0 .lt. kb) go to 190
    k2 = k2 - 1
    ak2 = a(k2+1)
    a(k2+1) = a(k0+1)-ak2
    a(k0+1) = a(k0+1)+ak2
    go to 85
90 if (l1 .eq. 1) go to 95
    c2 = c1 * c1 - s1 * s1
    s2 = two * c1 * s1
95 ja = kb + isp - 1
    ka = ja + kb
    ikb = kb+1
    ija = ja+1
    do 110 ii = ikb,ija
        k0 = ka - ii + 1
        k1 = k0 + isp
        k2 = k1 + isp
        za0 = a(k0+1)
        if (l1 .eq. 1) go to 100
        za4 = a(k1+1)
        a1 = a4*c1-b4*s1
        b1 = a4*s1+b4*c1
        za4 = a(k2+1)
        a2 = a4*c2-b4*s2
        b2 = a4*s2+b4*c2
        go to 105
100    za1 = a(k1+1)
        za2 = a(k2+1)
105    a(k0+1) = cmplx(a0+a1+a2,b0+b1+b2)
        a0 = -half * (a1+a2) + a0
        a1 = (a1-a2) * c30
        b0 = -half * (b1+b2) + b0
        b1 = (b1-b2) * c30
        a(k1+1) = cmplx(a0-b1,b0+a1)
        a(k2+1) = cmplx(a0+b1,b0-a1)
110 continue
    go to 190
115 if (l1 .eq. 1) go to 120
    c2 = c1 * c1 - s1 * s1
    s2 = two * c1 * s1
    c3 = c1 * c2 - s1 * s2
    s3 = s1 * c2 + c1 * s2
120 ja = kb + isp - 1
    ka = ja + kb
    ikb = kb+1
    ija = ja+1
    do 135 ii = ikb,ija
        k0 = ka - ii + 1
        k1 = k0 + isp
        k2 = k1 + isp
        k3 = k2 + isp
        za0 = a(k0+1)
        if (l1 .eq. 1) go to 125
        za4 = a(k1+1)
        a1 = a4*c1-b4*s1
        b1 = a4*s1+b4*c1
        za4 = a(k2+1)
        a2 = a4*c2-b4*s2
        b2 = a4*s2+b4*c2
        za4 = a(k3+1)
        a3 = a4*c3-b4*s3
        b3 = a4*s3+b4*c3
        go to 130
125    za1 = a(k1+1)
        za2 = a(k2+1)
        za3 = a(k3+1)
130    a(k0+1) = cmplx(a0+a2+a1+a3,b0+b2+b1+b3)
        a(k1+1) = cmplx(a0+a2-a1-a3,b0+b2-b1-b3)
        a(k2+1) = cmplx(a0-a2-b1+b3,b0-b2+a1-a3)
        a(k3+1) = cmplx(a0-a2+b1-b3,b0-b2-a1+a3)
135 continue
    go to 190
140 jk = kf - 1
    kh = jk/2
    k3 = iwk(id+i-1)
    k0 = kb + isp
    if (l1 .eq. 1) go to 150
    k = jk - 1
    wk(icf+1) = c1
    wk(isf+1) = s1
    do 145 j = 1,k
        wk(icf+j+1) = wk(icf+j) * c1 - wk(isf+j) * s1
        wk(isf+j+1) = wk(icf+j) * s1 + wk(isf+j) * c1
145 continue
150 if (kf .eq. jf) go to 160
    c2 = wk(icc+i)
    wk(ick+1) = c2
    wk(ick+jk) = c2
    s2 = wk(iss+i)
    wk(isk+1) = s2
    wk(isk+jk) = -s2
    do 155 j = 1,kh
        k = jk - j
        wk(ick+k) = wk(ick+j) * c2 - wk(isk+j) * s2
        wk(ick+j+1) = wk(ick+k)
        wk(isk+j+1) = wk(ick+j) * s2 + wk(isk+j) * c2
        wk(isk+k) = -wk(isk+j+1)
155 continue
160 k0 = k0 - 1
    k1 = k0
    k2 = k0 + k3
    za0 = a(k0+1)
    a3 = a0
    b3 = b0
    do 175 j = 1,kh
        k1 = k1 + isp
        k2 = k2 - isp
        if (l1 .eq. 1) go to 165
        k = kf - j
        za4 = a(k1+1)
        a1 = a4*wk(icf+j)-b4*wk(isf+j)
        b1 = a4*wk(isf+j)+b4*wk(icf+j)
        za4 = a(k2+1)
        a2 = a4*wk(icf+k)-b4*wk(isf+k)
        b2 = a4*wk(isf+k)+b4*wk(icf+k)
        go to 170
165    za1 = a(k1+1)
        za2 = a(k2+1)
170    wk(iap+j) = a1 + a2
        wk(iam+j) = a1 - a2
        wk(ibp+j) = b1 + b2
        wk(ibm+j) = b1 - b2
        a3 = a1 + a2 + a3
        b3 = b1 + b2 + b3
175 continue
    a(k0+1) = cmplx(a3,b3)
    k1 = k0
    k2 = k0 + k3
    do 185 j = 1,kh
        k1 = k1 + isp
        k2 = k2 - isp
        jk = j
        a1 = a0
        b1 = b0
        a2 = zero
        b2 = zero
        do 180  k = 1,kh
        a1 = a1 + wk(iap+k) * wk(ick+jk)
        a2 = a2 + wk(iam+k) * wk(isk+jk)
        b1 = b1 + wk(ibp+k) * wk(ick+jk)
        b2 = b2 + wk(ibm+k) * wk(isk+jk)
        jk = jk + j
        if (jk .ge. kf) jk = jk - kf
180    continue
        a(k1+1) = cmplx(a1-b2,b1+a2)
        a(k2+1) = cmplx(a1+b2,b1-a2)
185 continue
    if (k0 .gt. kb) go to 160
    jf = kf
190 if ( i .ge. mm ) go to 195
    i = i + 1
    go to 55
195 i = mm
    l1 = 0
    kb = iwk(id+i-1) + kb
    if (kb .ge. kn) go to 215
200 jj = iwk(ic+i-2) + jj
    if (jj .lt. iwk(ic+i-1)) go to 205
    i = i - 1
    jj = jj - iwk(ic+i)
    go to 200
205 if (i .ne. mm) go to 210
    c2 = c1
    c1 = cm * c1 - sm * s1
    s1 = sm * c2 + cm * s1
    go to 70
210 if (iwk(ill+i) .eq. 1) i = i + 1
    go to 55
215 i = 1
    ja = kt - 1
    ka = ja + 1
    if(ja.lt.1) go to 225
    do 220  ii = 1,ja
        j = ka - ii
        iwk(j+1) = iwk(j+1) - 1
        i = iwk(j+1) + i
220 continue
!                                  the result is now permuted to
!                                  normal order.
225 if (kt .le. 0) go to 270
    j = 1
    i = 0
    kb = 0
230 k2 = iwk(id+j) + kb
    k3 = k2
    jj = iwk(ic+j-1)
    jk = jj
    k0 = kb + jj
    isp = iwk(ic+j) - jj
235 k = k0 + jj
240 za4 = a(k0+1)
    a(k0+1) = a(k2+1)
    a(k2+1) = za4
    k0 = k0 + 1
    k2 = k2 + 1
    if (k0 .lt. k) go to 240
    k0 = k0 + isp
    k2 = k2 + isp
    if (k0 .lt. k3) go to 235
    if (k0 .ge. k3 + isp) go to 245
    k0 = k0 - iwk(id+j) + jj
    go to 235
245 k3 = iwk(id+j) + k3
    if (k3 - kb .ge. iwk(id+j-1)) go to 250
    k2 = k3 + jk
    jk = jk + jj
    k0 = k3 - iwk(id+j) + jk
    go to 235
250 if (j .ge. kt) go to 260
    k = iwk(j+1) + i
    j = j + 1
255 i = i + 1
    iwk(ill+i) = j
    if (i .lt. k) go to 255
    go to 230
260 kb = k3
    if (i .le. 0) go to 265
    j = iwk(ill+i)
    i = i - 1
    go to 230
265 if (kb .ge. n) go to 270
    j = 1
    go to 230
270 jk = iwk(ic+kt)
    isp = iwk(id+kt)
    m = m - kt
    kb = isp/jk-2
    if (kt .ge. m-1 ) go to 9005
    ita = ill+kb+1
    itb = ita+jk
    idm1 = id-1
    ikt = kt+1
    im = m+1
    do 275 j = ikt,im
        iwk(idm1+j) = iwk(idm1+j)/jk
275 continue
    jj = 0
    do 290 j = 1,kb
        k = kt
280    jj = iwk(id+k+1) + jj
        if (jj .lt. iwk(id+k)) go to 285
        jj = jj - iwk(id+k)
        k = k + 1
        go to 280
285    iwk(ill+j) = jj
        if (jj .eq. j) iwk(ill+j) = -j
290 continue
!                                  determine the permutation cycles
!                                  of length greater than or equal
!                                  to two.
    do 300  j = 1,kb
        if (iwk(ill+j) .le. 0) go to 300
        k2 = j
295    k2 = iabs(iwk(ill+k2))
        if (k2 .eq. j) go to 300
        iwk(ill+k2) = -iwk(ill+k2)
        go to 295
300 continue
!                                  reorder a following the
!                                  permutation cycles
    i = 0
    j = 0
    kb = 0
    kn = n
305 j = j + 1
    if (iwk(ill+j) .lt. 0) go to 305
    k = iwk(ill+j)
    k0 = jk * k + kb
310 za4 = a(k0+i+1)
    wk(ita+i) = a4
    wk(itb+i) = b4
    i = i + 1
    if (i .lt. jk) go to 310
    i = 0
315 k = -iwk(ill+k)
    jj = k0
    k0 = jk * k + kb
320 a(jj+i+1) = a(k0+i+1)
    i = i + 1
    if (i .lt. jk) go to 320
    i = 0
    if (k .ne. j) go to 315
325 a(k0+i+1) = cmplx(wk(ita+i),wk(itb+i))
    i = i + 1
    if (i .lt. jk) go to 325
    i = 0
    if (j .lt. k2) go to 305
    j = 0
    kb = kb + isp
    if (kb .lt. kn) go to 305
9005 return
    end
