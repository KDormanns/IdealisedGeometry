      subroutine srffti (n,wsave)
      real wsave(1)
c
      if (n .eq. 1) return
c
      call srfti1 (n,wsave(n+1),wsave(2*n+1))
c
      return
      end
