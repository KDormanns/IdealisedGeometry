      subroutine scffti (n,wsave)
      real wsave(1)
c
      if (n .eq. 1) return
c
      iw1 = n+n+1
      iw2 = iw1+n+n
      call scfti1 (n,wsave(iw1),wsave(iw2))
c
      return
      end
