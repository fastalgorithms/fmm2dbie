      subroutine legepolz(z,n,zpols)
c
c     this routine is written to produce the same
c     output as legepols if z is real-valued.
c      
c     n is the highest degree term to compute,
c     so zpols should be length n+1
c      
c     stability is not guaranteed for arbitrary
c     complex-valued z, especially for large n
c      
      implicit real *8 (a-h,o-z)
      complex *16 z, zpols(0:*)
      complex *16 zim1,zim2, zi

      zpols(0) = 1
      if (n .eq. 0) return
      zpols(1) = z
      if (n .eq. 1) return
      
      zim1 = 1
      zi = z

      do i = 2,n
         zim2 = zim1
         zim1 = zi
         zi = ((2*i-1)*z*zim1-(i-1)*zim2)/i
         zpols(i) = zi
      enddo

      return
      end
c 
