      subroutine dsinqb (n,x,wsave)
      double precision x(1), wsave(1), xhold
c
      if (n .gt. 1) go to 101
      x(1) = 4.d0*x(1)
      return
c
  101 ns2 = n/2
      do 102 k=2,n,2
         x(k) = -x(k)
  102 continue
c
      call dcosqb (n,x,wsave)
c
      do 103 k=1,ns2
         kc = n-k
         xhold = x(k)
         x(k) = x(kc+1)
         x(kc+1) = xhold
  103 continue
c
      return
      end
