program dispersion2D
        implicit none
        
        integer n, it
        real v, pi, dx, fmax, k, ramda, vph_v, G_r, p, vgr_v
        !real, allocatable, dimension(:) ::

        n=300
        v=2.0
        fmax=1.5
        pi=acos(-1.)
        ramda=v/fmax
        k=2.*pi/ramda

        p=0.9
       !---------V_ph---------------------! 
        open(10, file='disper1d09.dat')
        do it=1,n
                G_r=it/1000.
                dx=G_r*ramda
                vph_v=1./G_r/pi/p*asin(p*sin(k*dx/2.))

                write(10,*) G_r, vph_v
        enddo

        close(10)
        


        !--------0.99 data------------------!
        !open(11, file='99per.dat')
        !do it=1,n
        !        G_r=it/1000.
        !        write(11,*) G_r, 0.99
        !enddo
        !close(11)


        !--------V_gr-----------------------!
!        open(12, file='gr_disper1d09.dat')
!        do it=1,n
!                G_r=it/1000.
!                vgr_v=cos(pi*G_r)/cos(asin(p*sin(pi*G_r)))
!
!                write(12,*) G_r, vgr_v
!        enddo
!
!        close(12)



end program
