      subroutine scfftb (n,c,wsave)
      real c(1), wsave(1)
c
      if (n .eq. 1) return
c
      iw1 = n+n+1
      iw2 = iw1+n+n
      call scftb1 (n,c,wsave,wsave(iw1),wsave(iw2))
c
      return
      end
