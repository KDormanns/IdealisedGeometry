      subroutine deffti (n,wsave)
      double precision wsave(1)
c
      if (n .eq. 1) return
c
      call defft1 (n,wsave(2*n+1),wsave(3*n+1))
c
      return
      end
