program tester
      implicit none
double precision :: q(2), o(2), f
      q(1)=1d0
      q(2)=2d0
      o(1)=3d0
      o(2)=4d0

      f = dot_product(q(:),o(:)*q(:))
write(*,*) f
end program tester
