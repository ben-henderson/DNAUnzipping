      program initial_condition
      implicit none
      integer*4 i,j,k,nbead,nbeadmax,ntype
      integer*4 nbondmax,nanglemax,nbond,nangle,ntypemax,nmax
      parameter(nbeadmax=100000)      
      parameter(nbondmax=200000,nanglemax=200000)
      integer*4 seed
      double precision lx,ly,lz
      double precision theta,a1,ran2,pi
      double precision polymer(nbeadmax,3)
      integer*4 bond(nbondmax,2),angle(nbondmax,3)
      integer*4 atom_type(2*nbeadmax)
      external ran2

      pi=4.d0*datan(1.d0)
      lx=100.d0
      ly=100.d0
      lz=100.d0

      write(6,*) 'lx,ly,lz = ' 
      read(5,*) lx,ly,lz
      write(6,*) 'n beads in each strand (max = ',nbeadmax,' )'
      read(5,*) nbead

      seed=-1

      do k=1,3
         polymer(1,k)=0.d0
      enddo
      
      atom_type(1)=1

!     initialise double stranded DNA as two 2D random walks,
!     one on top of each other along z

      do i=2,nbead
         atom_type(i)=i
 5       a1=ran2(seed)
         theta=a1*2.d0*pi
         polymer(i,1)=polymer(i-1,1)+dcos(theta)
         polymer(i,2)=polymer(i-1,2)+dsin(theta)
         polymer(i,3)=polymer(i-1,3)
         if(dabs(polymer(i,1)).ge.lx/2.d0) goto 5
         if(dabs(polymer(i,2)).ge.ly/2.d0) goto 5
      enddo

      do i=nbead+1,2*nbead
         atom_type(i)=i-nbead
         do k=1,2
            polymer(i,k)=polymer(i-nbead,k)
         enddo
         polymer(i,k)=polymer(i-nbead,k)+1.d0
      enddo

      nbond=0
      nangle=0

      do i=1,nbead-1
         nbond=nbond+1
         bond(nbond,1)=i
         bond(nbond,2)=i+1
      enddo

      do i=nbead+1,2*nbead-1
         nbond=nbond+1
         bond(nbond,1)=i
         bond(nbond,2)=i+1
      enddo

      do i=1,nbead-2
         nangle=nangle+1
         angle(nangle,1)=i
         angle(nangle,2)=i+1
         angle(nangle,3)=i+2
      enddo

      do i=nbead+1,2*nbead-2
         nangle=nangle+1
         angle(nangle,1)=i
         angle(nangle,2)=i+1
         angle(nangle,3)=i+2
      enddo

      open(unit=7,file='lammps_input',status='unknown')

      write(7,*) 'LAMMPS data file from restart file'

      write(7,*) 2*nbead, '  atoms'
      write(7,*) nbond, '  bonds'
      write(7,*) nangle, '  angles'

      write(7,*) ''
      write(7,*) nbead,' atom types'
      write(7,*) '1 bond types'
      write(7,*) '1 angle types'
      write(7,*) ''

      write(7,*) -lx/2.d0,lx/2.d0,' xlo xhi'
      write(7,*) -ly/2.d0,ly/2.d0,' ylo yhi'
      write(7,*) -lz/2.d0,lz/2.d0,' zlo zhi'

      write(7,*) ''
      write(7,*) 'Masses'
      write(7,*) ''
      
      do j=1,nbead
         write(7,*) j,1
      enddo

      write(7,*) ''
      write(7,*) 'Atoms'
      write(7,*) ''

      do i=1,2*nbead
         write(7,*) i,1,atom_type(i),
     1           polymer(i,1),polymer(i,2),
     1           polymer(i,3),0,0,0
      enddo
      
      write(7,*) ''
      write(7,*) 'Velocities'
      write(7,*) ''

      do i=1,2*nbead
         write(7,*) i,0,0,0
      enddo

      write(7,*) ''
      write(7,*) 'Bonds'
      write(7,*) ''

      do i=1,nbond
         write(7,*) i,1,bond(i,1),bond(i,2)
      enddo

      write(7,*) ''
      write(7,*) 'Angles'
      write(7,*) ''

      do i=1,nangle
         write(7,*) i,1,angle(i,1),angle(i,2),angle(i,3)
      enddo

      stop
      end


      FUNCTION ran2(idum)
      implicit double precision (a-h,o-z)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *     IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *     NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1
                   k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
 11      continue                              
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
       ran2=min(AM*iy,RNMX)
      return
      END
C     (C) Copr. 1986-92 Numerical Recipes Software ')0.      




