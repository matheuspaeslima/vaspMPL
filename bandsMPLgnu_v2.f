      program bandsMPLgnu
      implicit none
      integer             :: nbands, nkp, nspin
      real*8, allocatable :: ener(:,:,:), kpoints(:,:), dkpoints(:,:)
      real*8              :: a0, cell(3,3), kcell(3,3), volume, pi2
      real*8              :: dk, Efermi
      real*8              :: dummy1, dummy2, dummy3, dummy4, dummy5
      character*128       :: cdummy
      character*128       :: filename, fname_out

c  .: working variables
      integer             :: ik, ibands, ispin

      print*, "provide the name of the output file"
      read(*,*) fname_out
      print*, "provide the Fermi energy:"
      read(*,*) Efermi

c  .: calculating reciprocal lattice vectors from POSCAR
      print*, "=>Calculating reciprocal lattice vectors from POSCAR"
      open(28,file='POSCAR',status="unknown")
      read(28,*) cdummy 
      read(28,*) a0
      read(28,*) cell(1,:)
      read(28,*) cell(2,:)
      read(28,*) cell(3,:)
      close(28)
      volume=cell(1,1)*(cell(2,2)*cell(3,3)-cell(2,3)*cell(3,2))
     .      +cell(1,2)*(cell(2,3)*cell(3,1)-cell(2,1)*cell(3,3))
     .      +cell(1,3)*(cell(2,1)*cell(3,2)-cell(2,2)*cell(3,1))
      pi2=2*dacos(-1.d0)
      kcell(1,1)=(cell(2,2)*cell(3,3)-cell(2,3)*cell(3,2))
      kcell(1,2)=(cell(2,3)*cell(3,1)-cell(2,1)*cell(3,3))
      kcell(1,3)=(cell(2,1)*cell(3,2)-cell(2,2)*cell(3,1))

      kcell(2,1)=(cell(3,2)*cell(1,3)-cell(3,3)*cell(1,2))
      kcell(2,2)=(cell(3,3)*cell(1,1)-cell(3,1)*cell(1,3))
      kcell(2,3)=(cell(3,1)*cell(1,2)-cell(3,2)*cell(1,1))

      kcell(3,1)=(cell(1,2)*cell(2,3)-cell(1,3)*cell(2,2))
      kcell(3,2)=(cell(1,3)*cell(2,1)-cell(1,1)*cell(2,3))
      kcell(3,3)=(cell(1,1)*cell(2,2)-cell(1,2)*cell(2,1))

      kcell(:,:)=(pi2/volume)*kcell(:,:)

      print'(3f14.8)',kcell(1,:)
      print'(3f14.8)',kcell(2,:)
      print'(3f14.8)',kcell(3,:)

c  .: reading EIGENVAL
      print*, "=>reading EIGENVAL"
      filename="EIGENVAL"

      open(29,file=filename,status="unknown")
     
      read(29,*) dummy1, dummy2, dummy3, nspin
      read(29,*) dummy1
      read(29,*) dummy1
      read(29,*) cdummy
      read(29,*) cdummy
      read(29,*) dummy1, nkp, nbands

c  .: allocating variables
      allocate(ener(nkp,nbands,nspin), kpoints(3,nkp),dkpoints(3,nkp))

      do ik=1, nkp
      read(29,*) 
      read(29,*) kpoints(:,ik), dummy1
c     write(*,*) kpoints(1,ik), kpoints(2,ik), kpoints(3,ik)
          do ibands=1,nbands
              read(29,*) dummy1, ener(ik,ibands,:) 
          enddo
          dkpoints(1,ik)=kpoints(1,ik)*kcell(1,1)
     .                  +kpoints(2,ik)*kcell(2,1)
     .                  +kpoints(3,ik)*kcell(3,1)
          dkpoints(2,ik)=kpoints(1,ik)*kcell(1,2)
     .                  +kpoints(2,ik)*kcell(2,2)
     .                  +kpoints(3,ik)*kcell(3,2)
          dkpoints(3,ik)=kpoints(1,ik)*kcell(1,3)
     .                  +kpoints(2,ik)*kcell(2,3)
     .                  +kpoints(3,ik)*kcell(3,3)
      enddo

c  .: writing EIGENVAL
      print*,"=>writting out the band structure file"
      open(30,file=fname_out,status='unknown')
      write(30,'((a),3I9)') "# nbands, nkp, nspin", nbands, nkp, nspin

      do ispin=1,nspin
      do ibands=1,nbands
          ik=1
c  .: bands with lattice vectors 
c  .:     dk is the starting point for the kpoint axis
          dk=0.0        
          write(30,'(2f18.8)') dk, ener(ik,ibands,ispin)-Efermi
          do ik=2,nkp
              dk=dk+sqrt((dkpoints(1,ik-1)-dkpoints(1,ik))**2 
     .                  +(dkpoints(2,ik-1)-dkpoints(2,ik))**2 
     .                  +(dkpoints(3,ik-1)-dkpoints(3,ik))**2) 
              write(30,'(2f18.8)') dk, ener(ik,ibands,ispin)-Efermi
          enddo
          write(30,*) "            "
      enddo
      enddo           
      close(30)
      print*, '=> exit successfully'

      end program bandsMPLgnu



