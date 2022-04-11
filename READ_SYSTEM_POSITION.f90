!    Copyright (C) 2012 The SPEED FOUNDATION
!    Author: Ilario Mazzieri
!
!    This file is part of SPEED.
!
!    SPEED is free software; you can redistribute it and/or modify it
!    under the terms of the GNU Affero General Public License as
!    published by the Free Software Foundation, either version 3 of the
!    License, or (at your option) any later version.
!
!    SPEED is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with SPEED.  If not, see <http://www.gnu.org/licenses/>.

!> @brief Reads oscillator position and writes SYSLST.input file
!! @author Aline Herlin, Srihari Sangaraju
!> @date February, 2021
!> @version 1.0

subroutine READ_SYSTEM_POSITION()

  use speed_par
  use SDOF_SYSTEM, only : SDOFnum

  implicit none

  include 'SPEED.MPI'

  if (sys_lst.eq.1) then
     if (nelem_loc .gt. 0) then
        allocate (highest_sys_lst_loc(nelem_loc))
    endif
    call GET_HIGHEST_NODE(nnod_loc, nelem_loc, zz_spx_loc,&
                         local_node_num, con_nnz_loc, con_spx_loc, &
                         nmat, tag_mat, sdeg_mat, highest_sys_lst_loc)
  endif


  if (mpi_id.eq.0) then
     write(*,'(A)')
     write(*,'(A)') '--------------------SDOF System location--------------------'
  endif


!*****************************************************************************************
!                      OSCILLATOR POSITION
!*****************************************************************************************

  SDOFnum = 0

  if (sys_lst .eq. 1) then

    file_SYS = 'SYS.input'
  	call READ_DIME_FILEPG(file_SYS,SDOFnum)

    if (SDOFnum.gt.0) then

      allocate(n_system_lst(SDOFnum),el_system_lst(SDOFnum),dist_system_lst(SDOFnum))
  	  allocate(system_label(SDOFnum),x_system_lst(SDOFnum),y_system_lst(SDOFnum),z_system_lst(SDOFnum))
  	  allocate(xr_system_lst(SDOFnum),yr_system_lst(SDOFnum),zr_system_lst(SDOFnum))

      call READ_FILESYS(file_SYS,SDOFnum,system_label,x_system_lst,y_system_lst,z_system_lst)

      !Modified, so that In a simulation both types of Point sources (usual point sources and  SDOF Oscillators) can be used (-SS)
      ! do i=1,SDOFnum
      !   do j=1,nload_poiX_el
      !     if(system_label(i).eq.building_id_x(j)) then
      !       val_poiX_el(j,1) = x_system_lst(i); val_poiX_el(j,2) = y_system_lst(i); val_poiX_el(j,3) = z_system_lst(i);
      !     endif
      !   enddo
      !   do j=1,nload_poiY_el
      !     if(system_label(i).eq.building_id_y(j)) then
      !       val_poiY_el(j,1) = x_system_lst(i); val_poiY_el(j,2) = y_system_lst(i); val_poiY_el(j,3) = z_system_lst(i);
      !     endif
      !   enddo
      !   do j=1,nload_poiZ_el
      !     if(system_label(i).eq.building_id_z(j)) then
      !       val_poiZ_el(j,1) = x_system_lst(i); val_poiZ_el(j,2) = y_system_lst(i); val_poiZ_el(j,3) = z_system_lst(i);
      !     endif
      !   enddo
      ! enddo

      if (file_sys_lst.eq.0) then ! NO input file with the position of oscillators
        allocate(x_system_real(SDOFnum), y_system_real(SDOFnum), z_system_real(SDOFnum))

  		  do i = 1,SDOFnum

          call GET_NEAREST_NODE_PGM(nnod_loc, xx_spx_loc, yy_spx_loc, zz_spx_loc,&
  					x_system_lst(i), y_system_lst(i), z_system_lst(i),&
  					n_system_lst(i), dist_system_lst(i), depth_search_sys_lst)

  			  call GET_PNT_POS_PGM(nelem_loc,&
  				  alfa11,alfa12,alfa13,alfa21,alfa22,alfa23,&
  					alfa31,alfa32,alfa33,beta11,beta12,beta13,&
  	 				beta21,beta22,beta23,beta31,beta32,beta33,&
  					gamma1,gamma2,gamma3,delta1,delta2,delta3,&
  					x_system_lst(i),y_system_lst(i),z_system_lst(i),&
  					el_system_lst(i),xr_system_lst(i),yr_system_lst(i),zr_system_lst(i),&
  					highest_sys_lst_loc, depth_search_sys_lst)


          x_system_real(i) = xx_spx_loc(n_system_lst(i))
          y_system_real(i) = yy_spx_loc(n_system_lst(i))
          z_system_real(i) = zz_spx_loc(n_system_lst(i))
          n_system_lst(i) = local_node_num(n_system_lst(i))
          el_system_lst(i) = local_el_num(el_system_lst(i))

  		  enddo

        j = 1

        do while(j .le. SDOFnum)
      	  allocate(dist_system_glo(mpi_np),n_system_glo(mpi_np), el_system_glo(mpi_np))
          allocate(xr_system_glo(mpi_np), yr_system_glo(mpi_np), zr_system_glo(mpi_np))
          allocate(x_system_glo_real(mpi_np), y_system_glo_real(mpi_np),z_system_glo_real(mpi_np))

          call MPI_BARRIER(mpi_comm, mpi_ierr)
          call MPI_ALLGATHER(dist_system_lst(j), 1, SPEED_DOUBLE, dist_system_glo, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
          call MPI_ALLGATHER(xr_system_lst(j), 1, SPEED_DOUBLE, xr_system_glo, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
          call MPI_ALLGATHER(yr_system_lst(j), 1, SPEED_DOUBLE, yr_system_glo, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
          call MPI_ALLGATHER(zr_system_lst(j), 1, SPEED_DOUBLE, zr_system_glo, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
          call MPI_ALLGATHER(n_system_lst(j), 1, SPEED_INTEGER, n_system_glo, 1, SPEED_INTEGER, mpi_comm, mpi_ierr)
          call MPI_ALLGATHER(el_system_lst(j), 1, SPEED_INTEGER, el_system_glo, 1, SPEED_INTEGER, mpi_comm, mpi_ierr)
          call MPI_ALLGATHER(x_system_real(j), 1, SPEED_DOUBLE, x_system_glo_real, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
          call MPI_ALLGATHER(y_system_real(j), 1, SPEED_DOUBLE, y_system_glo_real, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)
          call MPI_ALLGATHER(z_system_real(j), 1, SPEED_DOUBLE, z_system_glo_real, 1, SPEED_DOUBLE, mpi_comm, mpi_ierr)

          call GET_MINVALUES(n_system_glo, dist_system_glo, mpi_np, n_system_lst(j), 1, mpi_np)

          ! Do we really need this step? the partition containing nearest node is already saved in n_system_lst
          call GET_INDLOC_FROM_INDGLO(n_system_glo, mpi_np, n_system_glo(n_system_lst(j)), ic)

          x_system_real(j) = x_system_glo_real(ic)
          y_system_real(j) = y_system_glo_real(ic)
          z_system_real(j) = z_system_glo_real(ic)
          xr_system_lst(j) = xr_system_glo(ic)
          yr_system_lst(j) = yr_system_glo(ic)
          zr_system_lst(j) = zr_system_glo(ic)
          el_system_lst(j) = el_system_glo(ic)

          j=j+1

          deallocate(dist_system_glo, n_system_glo, el_system_glo, xr_system_glo, yr_system_glo, zr_system_glo, &
                     x_system_glo_real, y_system_glo_real, z_system_glo_real)
        enddo

        deallocate(dist_system_lst)

        if(mpi_id.eq. 0) then

          file_SYSLST = 'SLST.input'
  			  call WRITE_FILE_MPGM(file_SYSLST, SDOFnum, n_system_lst, el_system_lst, &
                               xr_system_lst, yr_system_lst, zr_system_lst)

          sys_filename  = 'SLST.position'
  			  call WRITE_FILE_MPGM(sys_filename, SDOFnum, n_system_lst, el_system_lst, &
                                                x_system_real, y_system_real, z_system_real)

        endif

        deallocate(x_system_real, y_system_real, z_system_real)

  	  else ! YES, it exists an input file with the position of LST monitors

  		  file_SYSLST = 'SLST.input'
  		  call READ_FILE_MPGM(file_SYSLST, SDOFnum, n_system_lst, el_system_lst, &
  				                  xr_system_lst, yr_system_lst, zr_system_lst)

  	  endif

	    deallocate(highest_sys_lst_loc)

    endif
  else
    write(*,'(A)') 'SYSLST key not found!'

  endif

  if (mpi_id.eq.0) then
    write(*,'(A)')
    write(*,'(A,I10)') 'SYStem positionLST: ',SDOFnum
    write(*,'(A,I2)') 'File SYSLST : ',file_sys_lst
  endif

  call MPI_BARRIER(mpi_comm, mpi_ierr)

  !-----------------------------------------------------------------------------------
  ! Mapping SDOF oscillators to nearest LGL nodes
  !------------------------------------------------------------------------------------
  call MAKE_SYSTEM_POSITION_LGLNODES(nnod_loc , local_node_num, con_nnz_loc, con_spx_loc, &
                                            xx_spx_loc, yy_spx_loc, zz_spx_loc, SDOFnum, system_label, &
                                            x_system_lst, y_system_lst, z_system_lst, &
                                            mpi_id, mpi_comm, mpi_np, node_counter_sdof)

  !----------------------------------------------------------------------------------
  !	WRITING SDOF_SYSTEM.INFO
  !----------------------------------------------------------------------------------

  if(mpi_id .eq. 0) then

    allocate(system_files(mpi_np))
    system_files = 0

    do i = 1, SDOFnum
       !write(*,*) 'el_system_lst(i) ', el_system_lst(i), 'elem_domain ', elem_domain(el_system_lst(i)), 'system files ', system_files(elem_domain(el_system_lst(i))+1)
       system_files(elem_domain(el_system_lst(i))+1) = system_files(elem_domain(el_system_lst(i))+1) + 1
    enddo

    if(len_trim(monitor_file) .ne. 70) then
       monitor_file_new = monitor_file(1:len_trim(monitor_file)) // '/SDOF_SYSTEM.INFO'
    else
       monitor_file_new = 'SDOF_SYSTEM.INFO'
    endif


    open(unit=50,file=monitor_file_new)
    write(50,*) mpi_np
    write(50,*) SDOFnum
    do i = 1, mpi_np
       write(50,*) system_files(i)
    enddo
    close(50)

    deallocate(system_files)
  endif


end subroutine READ_SYSTEM_POSITION
