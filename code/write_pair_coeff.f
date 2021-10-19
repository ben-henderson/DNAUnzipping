      program writepaircoeff
      implicit none
      integer*4 i,j,ntype,ntypemax
      integer*4 ntypeprot,nchromosomes,ntypeDNAprotein
      parameter(ntypemax=10000)
      integer*4 cellcycle,soft
      integer*4 type1,type2
      double precision sigma,epsilon,cutoff,newepsilon
      double precision size(ntypemax)


!      write(6,*) 'n types = '
!      read(5,*) ntypeprot
!     cellcycle=0 interphase; cellcycle=1 mitosis; 
!     cellcycle=2 mitosis with topo II;
!     cellcycle=-1 only excluded volume
!     ntypeprot=number of protein types
!     nchromosomes=number of chains/polymers/chromosomes
!     ntypeDNA=number of types of DNA beads which stick to protein
      ntype=200!ntypeprot
      do i=1,ntype
         size(i)=1.d0
      enddo

      do i=1,ntype
         do j=i,ntype
            sigma=(size(i)+size(j))/2.d0
         enddo
      enddo

      do i=1,ntype
            epsilon=1.4d0
            cutoff=2.5d0!(2.d0)**(1.d0/6.d0)*sigma
            write(6,*) 'pair_coeff ',i,i,epsilon,sigma,cutoff
      enddo

      stop
      end
