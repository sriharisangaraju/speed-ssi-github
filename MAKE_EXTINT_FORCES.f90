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

!> @brief Computes external loads.
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0
!> @param[in] nnod_loc local node number
!> @param[in] xs_loc x-coordinate of GLL nodes
!> @param[in] ys_loc y-coordinate of GLL nodes
!> @param[in] zs_loc z-coordinate of GLL nodes   
!> @param[in] local_n_num local node numeration
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc local spectral connectivity vector
!> @param[in] nm number of materials
!> @param[in] tag_mat material labels
!> @param[in] type_mat material type (dummy)
!> @param[in] sdeg_mat polynomial degree vector
!> @param[in] tref_mat dummy
!> @param[in] prop_mat  material properties (rho, lambda, mu, gamma)
!> @param[in] ne_loc number of local elements
!> @param[in] local_el_num local element numeration  
!> @param[in] alfa11 costant values for the bilinear map
!> @param[in] alfa12 costant values for the bilinear map
!> @param[in] alfa13 costant values for the bilinear map
!> @param[in] alfa21 costant values for the bilinear map
!> @param[in] alfa22 costant values for the bilinear map
!> @param[in] alfa23 costant values for the bilinear map
!> @param[in] alfa31 costant values for the bilinear map
!> @param[in] alfa32 costant values for the bilinear map
!> @param[in] alfa33 costant values for the bilinear map
!> @param[in] beta11 costant values for the bilinear map
!> @param[in] beta12 costant values for the bilinear map
!> @param[in] beta13 costant values for the bilinear map
!> @param[in] beta21 costant values for the bilinear map
!> @param[in] beta22 costant values for the bilinear map
!> @param[in] beta23 costant values for the bilinear map
!> @param[in] beta31 costant values for the bilinear map
!> @param[in] beta32 costant values for the bilinear map
!> @param[in] beta33 costant values for the bilinear map
!> @param[in] gamma1 costant values for the bilinear map
!> @param[in] gamma2 costant values for the bilinear map
!> @param[in] gamma3 costant values for the bilinear map
!> @param[in] delta1 costant values for the bilinear map
!> @param[in] delta2 costant values for the bilinear map
!> @param[in] delta3 costant values for the bilinear map
!> @param[in] cs_nnz_bc_loc length of cs_bc_loc
!> @param[in] cs_bc_loc local spectral boundary connectivity vector
!> @param[in] nl_***  number of load of type ***
!> @param[in] val_***  value for load type ***
!> @param[in] fun_*** func value for load type ***
!> @param[in] tag_***  tag for load type ***
!> @param[in] n_test  number of functions in test mode
!> @param[in] fun_test function for test mode
!> @param[in] nfunc  number of functions
!> @param[in] tag_func label for functions
!> @param[in] func_type function type
!> @param[in] func_indx index for functions 
!> @param[in] func_data function data
!> @param[in] nfunc_data number of function data for each function
!> @param[in] nhexa number of hexahedras
!> @param[in] con_hexa connectivity matrix
!> @param[in] length_cns  lenght check seismic nodes
!> @param[in] max_num max number of seismic nodes
!> @param[in] nl_sism number of seismic loads
!> @param[in] sour_ns infos about  seismic nodes  
!> @param[in] facsmom  seismic moment factor
!> @param[in] tausmom rise time for seismic moment
!> @param[in] length_cne  lenght check explosive nodes
!> @param[in] max_num_ne max number of explosive nodes
!> @param[in] nl_expl number of explosive loads
!> @param[in] sour_ne infos about explosive nodes  
!> @param[in] facsexpl explosive moment factor
!> @param[in] mpi_comm  mpi common world
!> @param[in] mpi_np   number of mpi process
!> @param[in] testmode 1 if test mode is active, 0 otherwise
!> @param[out] fmat vector of external applied loads


       subroutine MAKE_EXTINT_FORCES(nnod_loc,xs_loc,ys_loc,zs_loc,local_n_num,cs_nnz_loc,cs_loc,&
                          nm,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
                          ne_loc,local_el_num, &
                          alfa11,alfa12,alfa13,alfa21,alfa22,alfa23,&
                          alfa31,alfa32,alfa33,beta11,beta12,beta13,&
                          beta21,beta22,beta23,beta31,beta32,beta33,&
                          gamma1,gamma2,gamma3,delta1,delta2,delta3,&
                          cs_nnz_bc_loc,cs_bc_loc,&
                          nl_dirX,val_dirX,fun_dirX,tag_dirX,&
                          nl_dirY,val_dirY,fun_dirY,tag_dirY,&
                          nl_dirZ,val_dirZ,fun_dirZ,tag_dirZ,&
                          nl_neuX,val_neuX,fun_neuX,tag_neuX,&
                          nl_neuY,val_neuY,fun_neuY,tag_neuY,&
                          nl_neuZ,val_neuZ,fun_neuZ,tag_neuZ,&
                          nl_neuN,val_neuN,fun_neuN,tag_neuN,&                 
                          nl_poiX,val_poiX,fun_poiX,&
                          nl_poiY,val_poiY,fun_poiY,&
                          nl_poiZ,val_poiZ,fun_poiZ,&
                          nl_traX,val_traX,fun_traX,&
                          nl_traY,val_traY,fun_traY,&
                          nl_traZ,val_traZ,fun_traZ,&
                          nl_plaX,val_plaX,fun_plaX,tag_plaX,&                 
                          nl_plaY,val_plaY,fun_plaY,tag_plaY,&                 
                          nl_plaZ,val_plaZ,fun_plaZ,tag_plaZ,&                 
                          srcmodflag,szsism, &
                          nl_sism,val_sism,fun_sism,tag_sism,&                 
                          nl_expl,val_expl,fun_expl,tag_expl,&                 
                          nl_forX,val_forX,fun_forX,&
                          nl_forY,val_forY,fun_forY,&
                          nl_forZ,val_forZ,fun_forZ,&
                          nl_forc,val_forc,fun_forc,&
                          nl_pres,val_pres,fun_pres,&
                          nl_shea,val_shea,fun_shea,&
                          n_test,fun_test, & !val_fun_test,  &
                          nfunc,tag_func,func_type,func_indx,func_data,nfunc_data, &
                          fmat,&
                          con_hexa, nhexa,&                                 
                          length_cns,&                                        
                          sour_ns,max_num_ns,num_ns,&                        
                          facsmom,node_index_seq,&                        
                          tausmom,&                                        
                          length_cne,&                                        
                          sour_ne,max_num_ne,num_ne,&                        
                          facsexpl,&                                        
                          mpi_comm, mpi_np, mpi_id,testmode, &
                          nmat_nhe, rho_nhe, lambda_nhe, mu_nhe)
            
      use speed_exit_codes

      implicit none

      include 'SPEED.MPI'
     
      character*12 :: name_prop
      character*70 :: filename

      integer*4 :: nnod_loc,cs_nnz_loc,nm,ne_loc,cs_nnz_bc_loc
      integer*4 :: nl_dirX,nl_dirY,nl_dirZ,nl_neuX,nl_neuY,nl_neuZ
      integer*4 :: nl_neuN                                                  
      integer*4 :: nl_poiX,nl_poiY,nl_poiZ,nl_forX,nl_forY,nl_forZ
      integer*4 :: nl_plaX,nl_plaY,nl_plaZ,nl_traX,nl_traY,nl_traZ                                 
      integer*4 :: nl_sism, srcmodflag, szsism
      integer*4 :: nl_expl                                                 
      integer*4 :: nl_forc,nl_pres,nl_shea, n_test
      integer*4 :: nfunc, nfunc_data
      integer*4 :: im,ifun,ie,ip,il,nface_loc,nn,fn
      integer*4 :: is,in,id
      integer*4 :: i,j,k, val_st
      integer*4 :: ipl, nb_tra_load                                                                
      integer*4 :: nhexa                                                        
      integer*4 :: sum_node_bottom,sum_node_bottom1,sum_node_bottom2                
      integer*4 :: sum_node_bottom_first3                                        
      integer*4 :: last_node_bottom,last_node_bottom1,last_node_bottom2         
      integer*4 :: C,sit                                                        
      integer*4 :: isism,ipsism                                                 
      integer*4 :: length_cns, sum_length_cns                                        
      integer*4 :: max_num_ns                                                         
      integer*4 :: tt, ie_glob, ie_surf   
      integer*4 :: ic, ic1, ic2, ic3, ic4
      integer*4 :: mpi_comm, mpi_np, mpi_ierr, mpi_id
      integer*4 :: iexpl,ipexpl                                         
      integer*4 :: length_cne, sum_length_cne                           
      integer*4 :: max_num_ne                                                 
      integer*4 :: nnode_neuN,nelem_neuN                                
      integer*4 :: ne1,ne2,ne3,ne4                                        
      integer*4 :: index_vector,index, testmode, prova, nmat_nhe

      integer*4, dimension(:), allocatable :: i4count        
      integer*4, dimension(:), allocatable :: i4normal
      integer*4, dimension(:), allocatable :: ind_locX, ind_locY, ind_locZ
      integer*4, dimension(:), allocatable :: ind_gloX, ind_gloY, ind_gloZ
      integer*4, dimension(:), allocatable :: n_tra, n_tra_glo
      integer*4, dimension(nnod_loc) :: node_index_seq                   
      integer*4, dimension(nl_expl) :: num_ne                                 
      integer*4, dimension(6) :: num_sit                                        
      integer*4, dimension(nl_sism) :: num_ns                                         
      integer*4, dimension(nnod_loc) :: local_n_num
      integer*4, dimension(ne_loc) :: local_el_num
      integer*4, dimension(0:cs_nnz_loc) :: cs_loc
      integer*4, dimension(nm) :: tag_mat,type_mat,sdeg_mat
      integer*4, dimension(0:cs_nnz_bc_loc) :: cs_bc_loc
      integer*4, dimension(nl_dirX) :: fun_dirX, tag_dirX
      integer*4, dimension(nl_dirY) :: fun_dirY, tag_dirY
      integer*4, dimension(nl_dirZ) :: fun_dirZ, tag_dirZ
      integer*4, dimension(nl_neuX) :: fun_neuX, tag_neuX
      integer*4, dimension(nl_neuY) :: fun_neuY, tag_neuY
      integer*4, dimension(nl_neuZ) :: fun_neuZ, tag_neuZ
      integer*4, dimension(nl_neuN) :: fun_neuN, tag_neuN
      integer*4, dimension(nl_poiX) :: fun_poiX
      integer*4, dimension(nl_poiY) :: fun_poiY
      integer*4, dimension(nl_poiZ) :: fun_poiZ
      integer*4, dimension(nl_traX) :: fun_traX
      integer*4, dimension(nl_traY) :: fun_traY
      integer*4, dimension(nl_traZ) :: fun_traZ
      integer*4, dimension(nl_plaX) :: fun_plaX                         
      integer*4, dimension(nl_plaY) :: fun_plaY                         
      integer*4, dimension(nl_plaZ) :: fun_plaZ                         
      integer*4, dimension(nl_plaX) :: tag_plaX                                
      integer*4, dimension(nl_plaY) :: tag_plaY                         
      integer*4, dimension(nl_plaZ) :: tag_plaZ                         
      integer*4, dimension(nl_sism) :: fun_sism                         
      integer*4, dimension(nl_sism) :: tag_sism                         
      integer*4, dimension(nl_expl) :: fun_expl                         
      integer*4, dimension(nl_expl) :: tag_expl                         
      integer*4, dimension(nl_forX) :: fun_forX
      integer*4, dimension(nl_forY) :: fun_forY
      integer*4, dimension(nl_forZ) :: fun_forZ
      integer*4, dimension(nl_forc) :: fun_forc
      integer*4, dimension(nl_pres) :: fun_pres
      integer*4, dimension(n_test) :: fun_test
      integer*4, dimension(nl_shea) :: fun_shea
      integer*4, dimension(nfunc) :: tag_func
      integer*4, dimension(nl_poiX) :: node_poiX
      integer*4, dimension(nl_poiY) :: node_poiY
      integer*4, dimension(nl_poiZ) :: node_poiZ
      integer*4, dimension(nfunc) :: func_type
      integer*4, dimension(nfunc +1) :: func_indx
      

      integer*4, dimension(nhexa,9) :: con_hexa                                        
      integer*4, dimension(max_num_ns,nl_sism) :: sour_ns                        
      integer*4, dimension(max_num_ne,nl_expl) :: sour_ne                 

      integer*4, dimension(:), allocatable :: node_counter_poiX,node_counter_poiY,&
                                              node_counter_poiZ, recv_poiZ, &
                                              recv_poiY, recv_poiX
                                              
      real*8 :: dxdx,dxdy,dxdz,dydx,dydy,dydz,dzdx,dzdy,dzdz,det_j
      real*8 :: lambda,mu,alpha,tref,pi,lambda2, mu2, rho2
      real*8 :: rho,ellez                                                         
      real*8 :: l1x,l1y,l1z,l2x,l2y,l2z,area,v1,v2,v3,v4,v,term,rr
      real*8 :: x0,y0,z0,x1,x2,x3,r1,r2,r3,f1,f2,f3,phii,theta,psi
      real*8 :: dist, phi
      real*8 :: x,y,z
      real*8 :: csi, eta, zeta, normal_x,normal_y,normal_z
      real*8 :: slip1_sism,slip2_sism,slip3_sism                                 
      real*8 :: norm1_sism,norm2_sism,norm3_sism                                 
      real*8 :: amp_sism                                                         
      real*8 :: tau_sism                                                        
      real*8 :: slip1_expl,slip2_expl,slip3_expl                         
      real*8 :: norm1_expl,norm2_expl,norm3_expl                         
      real*8 :: amp_expl                                                 
      real*8 :: z_node_bottom,z_node_bottom1,z_node_bottom2
      real*8 :: a1, a2, a3, b1, b2, b3, c1, c2, c3, w1, cost1, cost2, cost3
      real*8 :: u1, u2, u3, u1_12, u1_13, u2_12, u2_23, u3_13, u3_23, k1, k2, k3
      real*8 :: vp,vs,omega, n1, n2, n3                        

      real*8, dimension(:), allocatable :: val_gloX, val_gloY, val_gloZ
      real*8, dimension(:), allocatable :: val_locX, val_locY, val_locZ
      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(:), allocatable :: normal_nx_el_neuN                
      real*8, dimension(:), allocatable :: normal_ny_el_neuN                
      real*8, dimension(:), allocatable :: normal_nz_el_neuN        
      real*8, dimension(:), allocatable :: dist_tra, dist_tra_glo
      real*8, dimension(:), allocatable :: x_tra,y_tra,z_tra,x_tra_glo,y_tra_glo,z_tra_glo
      real*8, dimension(:), allocatable :: xt_tra,yt_tra,zt_tra,xt_tra_glo,yt_tra_glo,zt_tra_glo
      real*8, dimension(:), allocatable :: dist_tra_real
      real*8, dimension(nnod_loc) :: xs_loc,ys_loc,zs_loc
      real*8, dimension(nnod_loc) :: lambda_nhe, rho_nhe, mu_nhe
      real*8, dimension(:,:,:), allocatable :: lambda_el, rho_el, mu_el, gamma_el
      real*8, dimension(nm) :: tref_mat
      real*8, dimension(ne_loc) :: alfa11,alfa12,alfa13, alfa21,alfa22,alfa23, alfa31,alfa32,alfa33
      real*8, dimension(ne_loc) :: beta11,beta12,beta13, beta21,beta22,beta23, beta31,beta32,beta33
      real*8, dimension(ne_loc) :: gamma1,gamma2,gamma3, delta1,delta2,delta3
      real*8, dimension(6) :: sum_facs

      real*8, dimension(:,:), allocatable :: dd
      real*8, dimension(nm,4) :: prop_mat
      real*8, dimension(nl_dirX,4) :: val_dirX
      real*8, dimension(nl_dirY,4) :: val_dirY
      real*8, dimension(nl_dirZ,4) :: val_dirZ
      real*8, dimension(nl_neuX,4) :: val_neuX
      real*8, dimension(nl_neuY,4) :: val_neuY
      real*8, dimension(nl_neuZ,4) :: val_neuZ
      real*8, dimension(nl_neuN,4) :: val_neuN                 
      real*8, dimension(nl_poiX,4) :: val_poiX
      real*8, dimension(nl_poiY,4) :: val_poiY
      real*8, dimension(nl_poiZ,4) :: val_poiZ
      real*8, dimension(nl_traX,4) :: val_traX
      real*8, dimension(nl_traY,4) :: val_traY
      real*8, dimension(nl_traZ,4) :: val_traZ
      real*8, dimension(nl_plaX,1) :: val_plaX                                 
      real*8, dimension(nl_plaY,1) :: val_plaY                                 
      real*8, dimension(nl_plaZ,1) :: val_plaZ                                 
      real*8, dimension(nl_sism,szsism) :: val_sism                         
      real*8, dimension(nl_expl,20) :: val_expl                         
      real*8, dimension(nl_forX,4) :: val_forX
      real*8, dimension(nl_forY,4) :: val_forY
      real*8, dimension(nl_forZ,4) :: val_forZ
      real*8, dimension(nl_forc,10) :: val_forc
      real*8, dimension(nl_pres,10) :: val_pres
      real*8, dimension(nl_shea,10) :: val_shea
      !real*8, dimension(n_test,10) :: val_fun_test
      real*8, dimension(nl_expl,6) :: facsexpl                                 
      real*8, dimension(nfunc,3*nnod_loc) :: fmat
      real*8, dimension(nl_sism,6) :: facsmom                                         
      real*8, dimension(nl_sism,1) :: tausmom                                         
      real*8, dimension(3,3) :: rot
      real*8, dimension(nfunc_data) :: func_data      


!***********************************************************************************
      pi = 4.d0*datan(1.d0);

!**********************************************************************************

! Initialization "seismic moment" variables - begin
      length_cns = 0                                
      if (nl_sism.gt.0) then                        
        do isism = 1,nl_sism                        
            facsmom(isism,1) = 0.0d0                 
            facsmom(isism,2) = 0.0d0                
            facsmom(isism,3) = 0.0d0                
            facsmom(isism,4) = 0.0d0                
            facsmom(isism,5) = 0.0d0                
            facsmom(isism,6) = 0.0d0                
            tausmom(isism,1) = 0.0d0                
        enddo                                        
      endif                                        
! Initialazation "seismic moment" variables - end

! Initialization "explosive source" variables - begin
      length_cne = 0                                
      if (nl_expl.gt.0) then                        
        do iexpl = 1,nl_expl                        
          facsexpl(iexpl,1) = 0.0d0                 
          facsexpl(iexpl,2) = 0.0d0                
          facsexpl(iexpl,3) = 0.0d0                
          facsexpl(iexpl,4) = 0.0d0                
          facsexpl(iexpl,5) = 0.0d0                
          facsexpl(iexpl,6) = 0.0d0                
        enddo                                        
      endif                                        
! Initialazation "explosive source" variables - end



!*************************************************************************************
! Calculation of the normal vector with respect to the surface element in which the
! Neumann load is enforced
! normal_nx_el_neuN, normal_ny_el_neuN, normal_nz_el_neuN =  are the 3 components
! of the normal vector
!*************************************************************************************

      if (nl_neuN.gt.0) then


        nelem_neuN = 0
        allocate(i4count(nnod_loc))
        i4count = 0

        call GET_NODE_FROM_FACE(nnod_loc,cs_nnz_bc_loc,cs_bc_loc,nl_neuN,tag_neuN,&
              nnode_neuN, i4count, local_n_num)



        ! nnode_neuN = node number in which the Neumann load is applied
        ! i4count(i) different from 0 if on the i-th the Neumann load is applied
        ! ne_loc = number of local elements
        ! ielem = global index for the element 

        do im = 1,nm
          nn = sdeg_mat(im) +1

          do ie = 1,ne_loc

            if (cs_loc(cs_loc(ie -1) +0).eq.tag_mat(im)) then
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)
             
            if ((i4count(ne1).ne.0).and.(i4count(ne2).ne.0) .and.(i4count(ne3).ne.0) .and.(i4count(ne4).ne.0)) then
               nelem_neuN = nelem_neuN +1
            endif

               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)
            
            if ((i4count(ne1).ne.0).and.(i4count(ne2).ne.0) .and.(i4count(ne3).ne.0) .and.(i4count(ne4).ne.0)) then
               nelem_neuN = nelem_neuN +1
            endif
            
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
            
            if ((i4count(ne1).ne.0).and.(i4count(ne2).ne.0) .and.(i4count(ne3).ne.0) .and.(i4count(ne4).ne.0)) then
             nelem_neuN = nelem_neuN +1
            endif
            
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
            
            if ((i4count(ne1).ne.0).and.(i4count(ne2).ne.0) .and.(i4count(ne3).ne.0) .and.(i4count(ne4).ne.0)) then
             nelem_neuN = nelem_neuN +1
            endif
            
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)
            
            if ((i4count(ne1).ne.0).and.(i4count(ne2).ne.0) .and.(i4count(ne3).ne.0) .and.(i4count(ne4).ne.0)) then
             nelem_neuN = nelem_neuN +1
            endif
            
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)
            
            if ((i4count(ne1).ne.0).and.(i4count(ne2).ne.0) .and.(i4count(ne3).ne.0) .and.(i4count(ne4).ne.0)) then
             nelem_neuN = nelem_neuN +1
            endif
            
            endif
          enddo
        enddo
            
            
            
         allocate(i4normal(nelem_neuN))
         allocate(normal_nx_el_neuN(nelem_neuN))
         allocate(normal_ny_el_neuN(nelem_neuN)) 
         allocate(normal_nz_el_neuN(nelem_neuN))
         
         normal_nx_el_neuN = 0.0d0
         normal_ny_el_neuN = 0.0d0
         normal_nz_el_neuN = 0.0d0
         nelem_neuN = 0
         
         do im = 1,nm
            nn = sdeg_mat(im) +1
            
            do ie = 1,ne_loc
            
            if (cs_loc(cs_loc(ie -1) +0).eq.tag_mat(im)) then
              ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
              ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
              ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)
              ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)
              
              if ((i4count(ne1).ne.0).and.(i4count(ne2).ne.0) .and.(i4count(ne3).ne.0) .and.(i4count(ne4).ne.0)) then
            
                  !face x = -1
                   
                   nelem_neuN = nelem_neuN +1
                  
                   call GET_ELEM_FROM_SURF(cs_nnz_bc_loc, cs_bc_loc, ne1, ne2, ne3, ne4, ie_surf)
                   if(ie_surf .eq. 0) then
                      write(*,*) 'surf not found'
                      call EXIT(EXIT_SURF_NOTFOUND)
                   endif   
                   i4normal(nelem_neuN) = ie_surf
                  
                   call MAKE_NORMAL(1,xs_loc(ne1), xs_loc(ne2), xs_loc(ne3), xs_loc(ne4), &
                                      ys_loc(ne1), ys_loc(ne2), ys_loc(ne3), ys_loc(ne4), &
                                      zs_loc(ne1), zs_loc(ne2), zs_loc(ne3), zs_loc(ne4) ,&
                                      normal_x,normal_y,normal_z, 0, 0 )
                  
                  
                   normal_nx_el_neuN(nelem_neuN) = normal_x
                   normal_ny_el_neuN(nelem_neuN) = normal_y
                   normal_nz_el_neuN(nelem_neuN) = normal_z
                  
                  
               endif
            
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)
               
               if ((i4count(ne1).ne.0).and.(i4count(ne2).ne.0) .and.(i4count(ne3).ne.0) .and.(i4count(ne4).ne.0)) then
               
               !  face y = -1                
                nelem_neuN = nelem_neuN +1
               
                call GET_ELEM_FROM_SURF(cs_nnz_bc_loc,cs_bc_loc, ne1, ne2, ne3, ne4, ie_surf)
                if(ie_surf .eq. 0) then
                   write(*,*) 'surf not found'
                   call EXIT(EXIT_SURF_NOTFOUND)
                endif   
                i4normal(nelem_neuN) = ie_surf
               
                call MAKE_NORMAL(2,xs_loc(ne1), xs_loc(ne2), xs_loc(ne3), xs_loc(ne4), &
                                   ys_loc(ne1), ys_loc(ne2), ys_loc(ne3), ys_loc(ne4), &
                                   zs_loc(ne1), zs_loc(ne2), zs_loc(ne3), zs_loc(ne4) ,&
                                   normal_x,normal_y,normal_z, 0, 0 )
               
               
                normal_nx_el_neuN(nelem_neuN) = normal_x
                normal_ny_el_neuN(nelem_neuN) = normal_y
                normal_nz_el_neuN(nelem_neuN) = normal_z
               
               endif
            
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
               
               if ((i4count(ne1).ne.0).and.(i4count(ne2).ne.0) .and.(i4count(ne3).ne.0) .and.(i4count(ne4).ne.0)) then
               
               !face z = -1
                nelem_neuN = nelem_neuN +1
               
                call GET_ELEM_FROM_SURF(cs_nnz_bc_loc,cs_bc_loc, ne1, ne2, ne3, ne4, ie_surf)
                if(ie_surf .eq. 0) then
                   write(*,*) 'surf not found'
                   call EXIT(EXIT_SURF_NOTFOUND)
                endif   
                i4normal(nelem_neuN) = ie_surf
               
                call MAKE_NORMAL(3,xs_loc(ne1), xs_loc(ne2), xs_loc(ne3), xs_loc(ne4), &
                                   ys_loc(ne1), ys_loc(ne2), ys_loc(ne3), ys_loc(ne4), &
                                   zs_loc(ne1), zs_loc(ne2), zs_loc(ne3), zs_loc(ne4) ,&
                                   normal_x,normal_y,normal_z, 0, 0 )
               
               
                normal_nx_el_neuN(nelem_neuN) = normal_x
                normal_ny_el_neuN(nelem_neuN) = normal_y
                normal_nz_el_neuN(nelem_neuN) = normal_z
               
               endif
               
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(1 -1) +(nn -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
               
               if ((i4count(ne1).ne.0).and.(i4count(ne2).ne.0) .and.(i4count(ne3).ne.0) .and.(i4count(ne4).ne.0)) then
               
               ! face x = 1
                
                nelem_neuN = nelem_neuN +1
               
                call GET_ELEM_FROM_SURF(cs_nnz_bc_loc,cs_bc_loc, ne1, ne2, ne3, ne4, ie_surf)
                if(ie_surf .eq. 0) then
                   write(*,*) 'surf not found'
                   call EXIT(EXIT_SURF_NOTFOUND)
                endif   
                
                i4normal(nelem_neuN) = ie_surf
               
                call MAKE_NORMAL(4,xs_loc(ne1), xs_loc(ne2), xs_loc(ne3), xs_loc(ne4), &
                                   ys_loc(ne1), ys_loc(ne2), ys_loc(ne3), ys_loc(ne4), &
                                   zs_loc(ne1), zs_loc(ne2), zs_loc(ne3), zs_loc(ne4) ,&
                                   normal_x,normal_y,normal_z, 0, 0 )
               
               
                normal_nx_el_neuN(nelem_neuN) = normal_x
                normal_ny_el_neuN(nelem_neuN) = normal_y
                normal_nz_el_neuN(nelem_neuN) = normal_z
               
               endif
               
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(1 -1) +nn*(nn -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)
               
               if ((i4count(ne1).ne.0).and.(i4count(ne2).ne.0) .and.(i4count(ne3).ne.0) .and.(i4count(ne4).ne.0)) then
               
               ! face y = 1
               
                nelem_neuN = nelem_neuN +1
               
                call GET_ELEM_FROM_SURF(cs_nnz_bc_loc,cs_bc_loc, ne1, ne2, ne3, ne4, ie_surf)
                if(ie_surf .eq. 0) then
                   write(*,*) 'surf not found'
                   call EXIT(EXIT_SURF_NOTFOUND)
                endif   
               
               
                i4normal(nelem_neuN) = ie_surf
               
                call MAKE_NORMAL(5,xs_loc(ne1), xs_loc(ne2), xs_loc(ne3), xs_loc(ne4), &
                                   ys_loc(ne1), ys_loc(ne2), ys_loc(ne3), ys_loc(ne4), &
                                   zs_loc(ne1), zs_loc(ne2), zs_loc(ne3), zs_loc(ne4) ,&
                                   normal_x,normal_y,normal_z, 0, 0 )
               
               
                normal_nx_el_neuN(nelem_neuN) = normal_x
                normal_ny_el_neuN(nelem_neuN) = normal_y
                normal_nz_el_neuN(nelem_neuN) = normal_z
               
               endif
               
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(1 -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(1 -1) +(nn -1) +1)
               ne3 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(nn -1) +1)
               ne4 = cs_loc(cs_loc(ie -1) +nn*nn*(nn -1) +nn*(nn -1) +(1 -1) +1)
               
               if ((i4count(ne1).ne.0).and.(i4count(ne2).ne.0) .and.(i4count(ne3).ne.0) .and.(i4count(ne4).ne.0)) then
               
               ! face z = 1
                nelem_neuN = nelem_neuN +1
               
                call GET_ELEM_FROM_SURF(cs_nnz_bc_loc,cs_bc_loc, ne1, ne2, ne3, ne4, ie_surf)
                if(ie_surf .eq. 0) then
                   write(*,*) 'surf not found'
                   call EXIT(EXIT_SURF_NOTFOUND)
                endif   
                i4normal(nelem_neuN) = ie_surf
               
                call MAKE_NORMAL(6,xs_loc(ne1), xs_loc(ne2), xs_loc(ne3), xs_loc(ne4), &
                                   ys_loc(ne1), ys_loc(ne2), ys_loc(ne3), ys_loc(ne4), &
                                   zs_loc(ne1), zs_loc(ne2), zs_loc(ne3), zs_loc(ne4) ,&
                                   normal_x,normal_y,normal_z, 0, 0 )
               
               
                normal_nx_el_neuN(nelem_neuN) = normal_x
                normal_ny_el_neuN(nelem_neuN) = normal_y
                normal_nz_el_neuN(nelem_neuN) = normal_z
               
               endif
            
            endif
            enddo
         enddo
            
         deallocate(i4count)
            
      endif   !if(nl_neuN.gt.0)
            
            ! *******************************************************************************************
            
            ! Fmat(nfunc, 3*nnode_dom) ---> Fmat(nfunc, nnod_loc)
            
      Fmat = 0.0d0
            
            !*********************************************************************************************
            ! If nl_poiX is different from 0 first compute the nearest local nodes to the one in which the
            ! load is applied and then compute the nearest global node by comparing the resutls on each 
            ! mpi process
            
      if (nl_poiX .gt. 0) then
            
            allocate(val_locX(nl_poiX), val_gloX(nl_poiX*mpi_np), ind_locX(nl_poiX), ind_gloX(nl_poiX*mpi_np))
            
            do i = 1,nl_poiX          
              call GET_NEAREST_NODE(nnod_loc,xs_loc,ys_loc,zs_loc, val_poiX(i,1),val_poiX(i,2),val_poiX(i,3),node_poiX(i), dist)
              val_locX(i) = dist
              ind_locX(i) = local_n_num(node_poiX(i))   !local_n_num(i)             
            enddo
            
            
            call MPI_BARRIER(mpi_comm,mpi_ierr)
            
            call MPI_ALLGATHER(val_locX, nl_poiX, SPEED_DOUBLE, val_gloX, nl_poiX, SPEED_DOUBLE, mpi_comm, mpi_ierr)
            call MPI_ALLGATHER(ind_locX, nl_poiX, SPEED_INTEGER, ind_gloX, nl_poiX, SPEED_INTEGER, mpi_comm, mpi_ierr)        
            
            ind_locX = 0
            
            call GET_MINVALUES(ind_gloX, val_gloX, nl_poiX*mpi_np, ind_locX, nl_poiX, mpi_np)
            
            do i = 1, nl_poiX
               node_poiX(i) = ind_gloX(ind_locX(i))
            enddo
            
            deallocate(ind_locX, val_locX, ind_gloX, val_gloX)        
            allocate(node_counter_poiX(nl_poiX),recv_poiX(nl_poiX))
            node_counter_poiX = 0;
            
            
            
            do ie = 1,ne_loc
               im = cs_loc(cs_loc(ie -1) +0);
               nn = sdeg_mat(im) +1
               do ip = 1, nl_poiX   
            
                  do k = 1,nn
                     do j = 1,nn
                        do i = 1,nn
                           is = nn*nn*(k -1) +nn*(j -1) +i
                           in = cs_loc(cs_loc(ie -1) +is)
            
                           if (local_n_num(in) .eq. node_poiX(ip)) then
                              node_counter_poiX(ip) = node_counter_poiX(ip) + 1  
                           endif
                        enddo
                     enddo
                  enddo
               enddo    
            
            enddo
            
            
            !      write(*,*) mpi_id, node_counter_poiX
            call MPI_BARRIER(mpi_comm,mpi_ierr)
            call MPI_ALLREDUCE(node_counter_poiX,recv_poiX, nl_poiX, SPEED_INTEGER, MPI_SUM, mpi_comm, mpi_ierr)
            node_counter_poiX = recv_poiX;
            deallocate(recv_poiX);
                 
      endif
            
             
      if (nl_poiY .gt. 0) then
            
            allocate(val_locY(nl_poiY), val_gloY(nl_poiY*mpi_np), ind_locY(nl_poiY), ind_gloY(nl_poiY*mpi_np))
            
            do i = 1,nl_poiY
               call GET_NEAREST_NODE(nnod_loc,xs_loc,ys_loc,zs_loc, val_poiY(i,1),val_poiY(i,2),val_poiY(i,3),node_poiY(i), dist)
               val_locY(i) = dist
               ind_locY(i) = local_n_num(node_poiY(i))   !local_n_num(i)             
            enddo
            
            call MPI_BARRIER(mpi_comm,mpi_ierr)
            
            call MPI_ALLGATHER(val_locY, nl_poiY, SPEED_DOUBLE, val_gloY, nl_poiY, SPEED_DOUBLE, mpi_comm, mpi_ierr)
            call MPI_ALLGATHER(ind_locY, nl_poiY, SPEED_INTEGER, ind_gloY, nl_poiY, SPEED_INTEGER, mpi_comm, mpi_ierr)                
            
            call GET_MINVALUES(ind_gloY,val_gloY, nl_poiY*mpi_np, ind_locY, nl_poiY, mpi_np)
            
            
            do i = 1, nl_poiY
               node_poiY(i) = ind_gloY(ind_locY(i))
            enddo
            
            deallocate(ind_locY, val_locY, ind_gloY, val_gloY)        
            allocate(node_counter_poiY(nl_poiY),recv_poiY(nl_poiY))
            node_counter_poiY = 0;
            
            
            
            do ie = 1,ne_loc
               im = cs_loc(cs_loc(ie -1) +0);
               nn = sdeg_mat(im) +1
               do ip = 1, nl_poiY   
            
                  do k = 1,nn
                     do j = 1,nn
                        do i = 1,nn
                           is = nn*nn*(k -1) +nn*(j -1) +i
                           in = cs_loc(cs_loc(ie -1) +is)
            
                           if (local_n_num(in) .eq. node_poiY(ip)) then
                              node_counter_poiY(ip) = node_counter_poiY(ip) + 1  
                           endif
                        enddo
                     enddo
                  enddo
               enddo    
            
            enddo
            
            
            
            call MPI_BARRIER(mpi_comm,mpi_ierr)
            call MPI_ALLREDUCE(node_counter_poiY,recv_poiY, nl_poiY, SPEED_INTEGER, MPI_SUM, mpi_comm, mpi_ierr)
            node_counter_poiY = recv_poiY;
            deallocate(recv_poiY);
            
            
      endif
            
      if (nl_poiZ .gt. 0) then
            
            allocate(val_locZ(nl_poiZ), val_gloZ(nl_poiZ*mpi_np), ind_locZ(nl_poiZ), ind_gloZ(nl_poiZ*mpi_np))
            
            do i = 1,nl_poiZ
               call GET_NEAREST_NODE(nnod_loc,xs_loc,ys_loc,zs_loc, val_poiZ(i,1),val_poiZ(i,2),val_poiZ(i,3),node_poiZ(i), dist)
               val_locZ(i) = dist
               ind_locZ(i) = local_n_num(node_poiZ(i))   !local_n_num(i)      
             enddo
            
            call MPI_BARRIER(mpi_comm,mpi_ierr)
            
            call MPI_ALLGATHER(val_locZ, nl_poiZ, SPEED_DOUBLE, val_gloZ, nl_poiZ, SPEED_DOUBLE, mpi_comm, mpi_ierr)
            call MPI_ALLGATHER(ind_locZ, nl_poiZ, SPEED_INTEGER, ind_gloZ, nl_poiZ, SPEED_INTEGER, mpi_comm, mpi_ierr)                
            call GET_MINVALUES(ind_gloZ, val_gloZ, nl_poiZ*mpi_np, ind_locZ, nl_poiZ, mpi_np)
 
            do i = 1, nl_poiZ
              node_poiZ(i) = ind_gloZ(ind_locZ(i))
            !    write(*,*) i, ind_locZ(i), ind_gloZ(ind_locZ(i)), node_poiZ(i)
            enddo
            ! write(*,*) node_poiZ
            
            deallocate(ind_locZ, val_locZ, ind_gloZ, val_gloZ)        
            allocate(node_counter_poiZ(nl_poiZ),recv_poiZ(nl_poiZ))
            node_counter_poiZ = 0;
            
            
            
            do ie = 1,ne_loc
               im = cs_loc(cs_loc(ie -1) +0);
               nn = sdeg_mat(im) +1
               do ip = 1, nl_poiZ   
            
                  do k = 1,nn
                     do j = 1,nn
                        do i = 1,nn
                           is = nn*nn*(k -1) +nn*(j -1) +i
                           in = cs_loc(cs_loc(ie -1) +is)

                           if (local_n_num(in) .eq. node_poiZ(ip)) then
                              node_counter_poiZ(ip) = node_counter_poiZ(ip) + 1  
                           endif
                        enddo
                     enddo
                  enddo
               enddo    
            
            enddo
            
            
            !write(*,*) mpi_id, node_counter_poiZ
            call MPI_BARRIER(mpi_comm,mpi_ierr)
            call MPI_ALLREDUCE(node_counter_poiZ, recv_poiZ, nl_poiZ, SPEED_INTEGER, MPI_SUM, mpi_comm, mpi_ierr)
            node_counter_poiZ = recv_poiZ;
            deallocate(recv_poiZ);
            !write(*,*) mpi_id, node_counter_poiZ
            call MPI_BARRIER(mpi_comm,mpi_ierr)
            
            !stop
      endif
            
      ne_loc = cs_loc(0) -1



      do ie = 1, ne_loc
               !write(*,*) ie, '/', ne_loc
         im = cs_loc(cs_loc(ie -1) +0);
               
         nn = sdeg_mat(im) +1
         allocate(ct(nn),ww(nn),dd(nn,nn))
         call MAKE_LGL_NW(nn,ct,ww,dd)
               
         rho = prop_mat(im,1) 
         lambda = prop_mat(im,2)
         mu = prop_mat(im,3)

         if (nmat_nhe.gt.0) then
                 allocate(rho_el(nn,nn,nn), lambda_el(nn,nn,nn), mu_el(nn,nn,nn), gamma_el(nn,nn,nn))
                 call GET_MECH_PROP_NH_ENHANCED(ie, nn, nnod_loc, cs_nnz_loc, cs_loc, &
                                                 rho_nhe, lambda_nhe, mu_nhe, 100.d0, 3.d0, &
                                                 rho_el, lambda_el, mu_el, gamma_el)
                  rho = sum(rho_el)/(nn*nn*nn)
                  lambda = sum(lambda_el)/(nn*nn*nn)
                  mu = sum(mu_el)/(nn*nn*nn)
                  deallocate(rho_el, lambda_el, mu_el, gamma_el)
         endif

         if (nl_poiX.gt.0) then    ! Point load X
                  do ip = 1,nl_poiX
                     fn = 0
                     do ifun = 1,nfunc
                        if (fun_poiX(ip).eq.tag_func(ifun)) fn = ifun
                     enddo
                     
                     if (fn.gt.0) then

                        call GET_INDLOC_FROM_INDGLO(local_n_num, nnod_loc, node_poiX(ip), ic1)

                        if (ic1 .ne. 0) then

                           do k = 1,nn
                              do j = 1,nn
                                 do i = 1,nn
                                    is = nn*nn*(k -1) +nn*(j -1) +i
                                    in = cs_loc(cs_loc(ie -1) +is)

                                    if (local_n_num(in) .eq. node_poiX(ip)) then
                                       fmat(fn,(3*(in -1) +1)) = fmat(fn,(3*(in -1) +1)) + val_poiX(ip,4)/node_counter_poiX(ip)  
                                       ! node_counter_poiX = node_counter_poiX + 1;                                 
                                       
                                       ! SDOF SS (Need to define, building ids corresponding to each node) 
                                    endif
                                 enddo
                              enddo
                           enddo
                        endif  

                     endif
                  enddo
         endif  
         
         if (nl_poiY.gt.0) then    ! Point load Y
                  do ip = 1,nl_poiY
                     fn = 0
                     do ifun = 1,nfunc
                        if (fun_poiY(ip).eq.tag_func(ifun)) fn = ifun
                     enddo
    
                     if (fn.gt.0) then

                        call GET_INDLOC_FROM_INDGLO(local_n_num, nnod_loc, node_poiY(ip), ic1)
                        if (ic1 .ne. 0) then
                           do k = 1,nn
                              do j = 1,nn
                                 do i = 1,nn
                                 
                                    is = nn*nn*(k -1) +nn*(j -1) +i
                                    in = cs_loc(cs_loc(ie -1) +is)
                                 
                                    if (local_n_num(in) .eq. node_poiY(ip)) then
                                        fmat(fn,(3*(in -1) +2)) = fmat(fn,(3*(in -1) +2)) + val_poiY(ip,4)/node_counter_poiY(ip)
                                    endif
                                 enddo
                              enddo
                           enddo
                        endif   
                     endif
                  enddo
         endif !(nl_poiY .gt. 0)
                            
         if (nl_poiZ.gt.0) then    ! Point load Z
            do ip = 1,nl_poiZ
                  fn = 0
               do ifun = 1,nfunc
                  if (fun_poiZ(ip).eq.tag_func(ifun)) fn = ifun
               enddo
    
               if (fn.gt.0) then

                  call GET_INDLOC_FROM_INDGLO(local_n_num, nnod_loc, node_poiZ(ip), ic1)
                  if (ic1 .ne. 0) then

                     do k = 1,nn
                        do j = 1,nn
                           do i = 1,nn

                              is = nn*nn*(k -1) +nn*(j -1) +i                        
                              in = cs_loc(cs_loc(ie -1) +is)                                        
                                               
                              if (local_n_num(in) .eq. node_poiZ(ip)) then
                      
                                 fmat(fn,(3*(in -1) +3)) = fmat(fn,(3*(in -1) +3)) &
                                           + val_poiZ(ip,4)/node_counter_poiZ(ip) 

                              endif

                           enddo
                        enddo
                     enddo
                  endif
               endif               
            enddo
         endif   !(nl_poiZ .gt. 0)


          !#############################################################################
          !##################                PLANE WAVE BEGIN             ##############
          !#############################################################################
 
         if (nl_plaX.gt.0) then   ! Plane Wave X
                              
            do ipl = 1,nl_plaX
                  
               if (tag_mat(im).eq.tag_plaX(ipl)) then  !Check on the Plane Wave Material - start
                                           
                  ! Recognizing bottom face - begin
                  sum_node_bottom1 = 0
                  sum_node_bottom2 = 0
                  last_node_bottom1 = 0
                  last_node_bottom2 = 0
                         
                  do i = 1,8
                     ie_glob = local_el_num(ie)
                             
                     call GET_INDLOC_FROM_INDGLO(local_n_num, nnod_loc, con_hexa(ie_glob,2), ic1)
                     call GET_INDLOC_FROM_INDGLO(local_n_num, nnod_loc, con_hexa(ie_glob,i+1), ic2)
                             
                     if(ic1 .ne. 0 .and. ic2 .ne. 0) then                   
                        if (dabs(zs_loc(ic1) - zs_loc(ic2)) .le. 1.d-4) then
                           last_node_bottom1 = i
                           sum_node_bottom1 = sum_node_bottom1 + i
                           z_node_bottom1 = zs_loc(ic1)
                  
                        else
                           last_node_bottom2 = i
                              sum_node_bottom2 = sum_node_bottom2 + i
                              z_node_bottom2 = zs_loc(ic2)
                                     
                        endif
                     endif   
                              
                  enddo
                          
                  if (z_node_bottom1.lt.z_node_bottom2) then
                     sum_node_bottom = sum_node_bottom1
                     sum_node_bottom_first3 = sum_node_bottom1 - last_node_bottom1
                   else
                     sum_node_bottom = sum_node_bottom2
                     sum_node_bottom_first3 = sum_node_bottom2 - last_node_bottom2
                  endif
                         
                  ! Recognizing bottom face - end
                  C=dsqrt(mu*rho)
                  ellez=dabs(z_node_bottom1 - z_node_bottom2)
                  
                  fn = 0
                  do ifun = 1,nfunc
                     if (fun_plaX(ipl).eq.tag_func(ifun)) fn = ifun
                  enddo
                  if (fn.gt.0) then
                  
                            !Table of the possible situation:
                            !hypothesis:
                            !- the plane wave load source is placed on a rectangular element, aligned 
                            !  along the principal axis (x,y,z). It does NOT work properly on element 
                            !  with different height (z co-ordinate MUST BE EQUAL)!!!
                            !- the node orientation is ALWAYS counter-clockwise!!!
                            !
                            !     Counter Clock wise means:     |  Looking from the top:
                            !                                   |  The number of nodes increses 
                            !         ^ Z                       |  going from the x axis versus y axis
                            !         |                         |  1 - 2 - 3 - 4
                            !         |                         |
                            !         5------8                  |
                            !        /|     /|                  |   Z
                            !       / |    / |                  |   o----------------->  Y
                            !      6------7  |                  |   |  
                            !      |  1---|--4 ------>  Y       |   |  1-------4
                            !      | /    | /                   |   |  |       |
                            !      |/     |/                    |   |  |       |
                            !      2------3                     |   |  |       |
                            !     /                             |   |  2-------3
                            !    /                              |   |
                            !   /                               |   v X
                            !  v  X                             |   
                            !
                            !
                            !  Possible position of the BOTTOM FACE
                            !
                            !         S1                    S2                   S3       
                            !                                                             
                            !      5------8              8------4             4------2    
                            !     /|     /|             /|     /|            /------/|    
                            !    / |    / |            / |    //|           /------/ |    
                            !   6------7  |           7------3//|          3------1  |    
                            !   |  1---|--4           |  5---|//2          |  8---|--5    
                            !   | /----|-/            | /    |//           | /    | /     
                            !   |/-----|/             |/     |/            |/     |/      
                            !   2------3              6------1             7------6     
                            !
                            !   sum_node_bottom
                            !   1+2+3+4=10           5+6+1+2=14            8+7+6+5=26
                            !
                            !         S4                    S5                   S6       
                            !                                                             
                            !      2------5              6------7             1------4    
                            !     /|     /|             /|     /|            /|-----/|    
                            !    //|    / |            / |    / |           / |----/-|    
                            !   1------6  |           2------3  |          5------8--|    
                            !   |//4---|--8           |- 5---|--8          |  2---|--3    
                            !   |//    | /            |------| /           | /    | /     
                            !   |/     |/             |------|/            |/     |/      
                            !   3------7              1------4             6------7     
                            !
                            !   
                            !   4+3+7+8=22           5+1+4+8=18            2+6+7+3=18
                            !
                            
                            !****** Situation S1 or S3 - begin ******
                            
                     if ((sum_node_bottom.eq.10).or.(sum_node_bottom.eq.26)) then
                              !
                               !*** Situation S1 or S3 ***
                               !
                               ! e.g.: spectral degree = 3, nn = 4
                               ! [#] = macro nodes
                               !  #  = micro nodes
                               !
                               ! i = 1,nn = 1,4
                               ! j = 1,nn = 1,4
                               ! z = ((nn-1)/2)+1 + num_sit(sit) = 2
                               !
                               !
                               ! *******************************************************
                               ! EXAMPLE DRAWING STILL ON WORKING!!!
                               !
                               !
                               !         ^ Z                  
                               !         |                    
                               !         |                    
                               !         5------8             
                               !        /|     /|             
                               !       / |    / |             
                               !      6------7  |             
                               !      |[1]
                               !         1---|--4 ------>  Y  
                               !      |/    | /              
                               !      2
                               !     /
                               !     3
                               !    /
                               !   4            
                               ! [2]-- 8--12--16[3]              
                               !      2------3                
                               !     /                        
                               !    /                         
                               !   /                          
                               !  v  X                        
                               !
                               !
                               !          S1
                               !   [2]4--3--2--1[1]
                               !    |            |
                               !    |8---7--6---5|   Applied plane wave load on this surface
                               !    |            |
                               !    |12--11-10--9|
                               !    |            |
                               !   [3]16-15-14-13[4]
                               !
                               ! e.g.: spectral degree = 4, nn = 5
                               !
                               !           S1
                               !   [2]5--4--3--2--1[1]
                               !    |               |
                               !    |10--9--8--7---6|
                               !    |               |
                               !  ->|15-14-13-12--11|<- Applied plane wave load
                               !    |               |
                               !    |20-19-18-17--16|
                               !    |               |
                               !   [3]25-24-23-22-21[4]
                               !
                               ! *******************************************************
                               ! EXAMPLE DRAWING STILL ON WORKING!!!
                               !
                               !
                               
                               !*** Situation S1 ***
                               
                        if (sum_node_bottom.eq.10) then
                               
                           sit = 1 
                            
                               !*** Situation S3 ***
                               
                         else
                               
                           sit = 2
                            
                        endif
                  
                        ! Even spectral polynomial degree
                        if (mod(nn-1,2).eq.0) then
                           k = ((nn-1)/2)+1
                                 
                         ! Odd spectral polynomial degree
                         else
                           k = (int((nn-1)/2))+num_sit(sit)
                        endif
                                 
                                 
                        do i = 1,nn
                           do j = 1,nn
                                       dxdx = alfa11(ie) +beta12(ie)*ct(k) &
                                              + beta13(ie)*ct(j) +gamma1(ie)*ct(j)*ct(k)
                                       dydx = alfa21(ie) +beta22(ie)*ct(k) &
                                              + beta23(ie)*ct(j) +gamma2(ie)*ct(j)*ct(k)
                                       dzdx = alfa31(ie) +beta32(ie)*ct(k) &
                                              + beta33(ie)*ct(j) +gamma3(ie)*ct(j)*ct(k)                 
                                       dxdy = alfa12(ie) +beta11(ie)*ct(k) &
                                              + beta13(ie)*ct(i) +gamma1(ie)*ct(k)*ct(i)
                                       dydy = alfa22(ie) +beta21(ie)*ct(k) &
                                              + beta23(ie)*ct(i) +gamma2(ie)*ct(k)*ct(i)
                                       dzdy = alfa32(ie) +beta31(ie)*ct(k) &
                                              + beta33(ie)*ct(i) +gamma3(ie)*ct(k)*ct(i)
                                       dxdz = alfa13(ie) +beta11(ie)*ct(j) &
                                              + beta12(ie)*ct(i) +gamma1(ie)*ct(i)*ct(j)
                                       dydz = alfa23(ie) +beta21(ie)*ct(j) &
                                              + beta22(ie)*ct(i) +gamma2(ie)*ct(i)*ct(j)
                                       dzdz = alfa33(ie) +beta31(ie)*ct(j) &
                                              + beta32(ie)*ct(i) +gamma3(ie)*ct(i)*ct(j)
                                       det_j = dxdz * (dydx*dzdy - dzdx*dydy) &
                                             - dydz * (dxdx*dzdy - dzdx*dxdy) &
                                             + dzdz * (dxdx*dydy - dydx*dxdy)
                                       is = nn*nn*(k -1) +nn*(j -1) +i
                                       in = cs_loc(cs_loc(ie -1) + is)
                             
                                       term = 4 * C * val_plaX(ipl,1) * ww(i) * ww (j) * dabs(det_j) / ellez
                                       fmat(fn,(3*(in -1) +1)) = fmat(fn,(3*(in -1) +1)) + term
                  
                                      prova = 3*(in -1) +1
                  
                           enddo
                        enddo
                               
                               
                              
                               
                            !****** Situation S2 or S4 - begin ******
                            
                      elseif ((sum_node_bottom.eq.14).or.(sum_node_bottom.eq.22)) then
                  
                             ! 
                             !*** Situation S2 ***
                               
                        if (sum_node_bottom.eq.14) then
                               
                                  sit = 1 
                            
                               !*** Situation S4 ***
                               
                         else
                               
                                  sit = 2
                            
                        endif
                  
                               ! Even spectral polynomial degree
                        if (mod(nn-1,2).eq.0) then
                           j = ((nn-1)/2)+1
                                 
                               ! Odd spectral polynomial degree
                         else
                           j = (int((nn-1)/2))+num_sit(sit)
                        endif
                                 
                        do i = 1,nn
                           do k = 1,nn
                                     dxdx = alfa11(ie) +beta12(ie)*ct(k) &
                                            + beta13(ie)*ct(j) +gamma1(ie)*ct(j)*ct(k)
                                     dydx = alfa21(ie) +beta22(ie)*ct(k) &
                                            + beta23(ie)*ct(j) +gamma2(ie)*ct(j)*ct(k)
                                     dzdx = alfa31(ie) +beta32(ie)*ct(k) &
                                            + beta33(ie)*ct(j) +gamma3(ie)*ct(j)*ct(k)
                                     dxdy = alfa12(ie) +beta11(ie)*ct(k) &
                                            + beta13(ie)*ct(i) +gamma1(ie)*ct(k)*ct(i)
                                     dydy = alfa22(ie) +beta21(ie)*ct(k) &
                                            + beta23(ie)*ct(i) +gamma2(ie)*ct(k)*ct(i)
                                     dzdy = alfa32(ie) +beta31(ie)*ct(k) &
                                            + beta33(ie)*ct(i) +gamma3(ie)*ct(k)*ct(i)
                                     dxdz = alfa13(ie) +beta11(ie)*ct(j) &
                                            + beta12(ie)*ct(i) +gamma1(ie)*ct(i)*ct(j)
                                     dydz = alfa23(ie) +beta21(ie)*ct(j) &
                                            + beta22(ie)*ct(i) +gamma2(ie)*ct(i)*ct(j)
                                     dzdz = alfa33(ie) +beta31(ie)*ct(j) &
                                            + beta32(ie)*ct(i) +gamma3(ie)*ct(i)*ct(j)
                                     det_j = dxdz * (dydx*dzdy - dzdx*dydy) &
                                            - dydz * (dxdx*dzdy - dzdx*dxdy) &
                                            + dzdz * (dxdx*dydy - dydx*dxdy)
                  
                                     is = nn*nn*(k -1) +nn*(j -1) +i
                                     in = cs_loc(cs_loc(ie -1) +is)
                  
                                     term = 4 * C * val_plaX(ipl,1) * ww(i) * ww (k) * dabs(det_j) / ellez
                  
                                     prova = 3*(in -1) +1
                                      
                  
                                     fmat(fn,(3*(in -1) +1)) = fmat(fn,(3*(in -1) +1)) + term
                           enddo
                        enddo
                              
                              
                              
                            !****** Situation S5 or S6 - begin ******
                            
                      elseif (sum_node_bottom.eq.18) then
                  
                             !*** Situation S5 ***
                               
                        if ((sum_node_bottom_first3.eq.10) .or.(sum_node_bottom_first3.eq.13) &
                                  .or.(sum_node_bottom_first3.eq.14)) then
                               
                                  sit = 1 
                            
                               !*** Situation S6 ***
                               
                               else
                               
                                  sit = 2
                            
                        endif
                  
                               ! Even spectral polynomial degree
                        if (mod(nn-1,2).eq.0) then
                                 i = ((nn-1)/2)+1
                                 
                               ! Odd spectral polynomial degree
                               else
                                 i = (int((nn-1)/2))+num_sit(sit)
                        endif
                                 
                              do j = 1,nn
                                 do k = 1,nn
                                      dxdx = alfa11(ie) +beta12(ie)*ct(k) &
                                             + beta13(ie)*ct(j) +gamma1(ie)*ct(j)*ct(k)
                                      dydx = alfa21(ie) +beta22(ie)*ct(k) &
                                             + beta23(ie)*ct(j) +gamma2(ie)*ct(j)*ct(k)
                                      dzdx = alfa31(ie) +beta32(ie)*ct(k) &
                                             + beta33(ie)*ct(j) +gamma3(ie)*ct(j)*ct(k)
                                      dxdy = alfa12(ie) +beta11(ie)*ct(k) &
                                             + beta13(ie)*ct(i) +gamma1(ie)*ct(k)*ct(i)
                                      dydy = alfa22(ie) +beta21(ie)*ct(k) &
                                             + beta23(ie)*ct(i) +gamma2(ie)*ct(k)*ct(i)
                                      dzdy = alfa32(ie) +beta31(ie)*ct(k) &
                                             + beta33(ie)*ct(i) +gamma3(ie)*ct(k)*ct(i)
                                      dxdz = alfa13(ie) +beta11(ie)*ct(j) &
                                             + beta12(ie)*ct(i) +gamma1(ie)*ct(i)*ct(j)
                                      dydz = alfa23(ie) +beta21(ie)*ct(j) &
                                             + beta22(ie)*ct(i) +gamma2(ie)*ct(i)*ct(j)
                                      dzdz = alfa33(ie) +beta31(ie)*ct(j) &
                                             + beta32(ie)*ct(i) +gamma3(ie)*ct(i)*ct(j)
                                      det_j = dxdz * (dydx*dzdy - dzdx*dydy) &
                                                    - dydz * (dxdx*dzdy - dzdx*dxdy) &
                                                    + dzdz * (dxdx*dydy - dydx*dxdy)
                  
                                      is = nn*nn*(k -1) +nn*(j -1) +i
                                      in = cs_loc(cs_loc(ie -1) +is)
                  
                                      term = 4 * C * val_plaX(ipl,1) * ww(j) * ww (k) * dabs(det_j) / ellez
                                      fmat(fn,(3*(in -1) +1)) = fmat(fn,(3*(in -1) +1)) + term
                                      
                                      
                                 enddo
                        enddo
                              
                     endif 
                               
                  
                  endif !if (fn.gt.0)
               endif  !Check on the Plane Wave Material - end
                  
            enddo !ipl = 1,nl_plaX
         endif  ! Plane Wave X - end                
                  
         if (nl_plaY.gt.0) then       ! Plane Wave Y
                   do ipl = 1,nl_plaY
                  
                      if (tag_mat(im).eq.tag_plaY(ipl)) then  !Check on the Plane Wave Material - start
                      
                         ! Recognizing bottom face - begin
                         sum_node_bottom1 = 0
                         sum_node_bottom2 = 0
                         last_node_bottom1 = 0
                         last_node_bottom2 = 0
                  
                         do i = 1,8
                             ie_glob = local_el_num(ie)
                             call GET_INDLOC_FROM_INDGLO(local_n_num, nnod_loc, con_hexa(ie_glob,2), ic1)
                             call GET_INDLOC_FROM_INDGLO(local_n_num, nnod_loc, con_hexa(ie_glob,i+1), ic2)
                             
                             if(ic1 .ne. 0 .and. ic2 .ne. 0) then                           
                                 if (dabs(zs_loc(ic1) - zs_loc(ic2)) .le. 1.d-4) then
                                 !if (zs_loc(ic1) .eq. zs_loc(ic2)) then
                                     last_node_bottom1 = i
                                     sum_node_bottom1 = sum_node_bottom1 + i
                                     z_node_bottom1 = zs_loc(ic1)
                                 else
                                     last_node_bottom2 = i
                                     sum_node_bottom2 = sum_node_bottom2 + i
                                     z_node_bottom2 = zs_loc(ic2)
                                 endif
                              endif   
                          enddo
                          
                          if (z_node_bottom1.lt.z_node_bottom2) then
                                 sum_node_bottom = sum_node_bottom1
                                 sum_node_bottom_first3 = sum_node_bottom1 - last_node_bottom1
                          else
                                 sum_node_bottom = sum_node_bottom2
                                 sum_node_bottom_first3 = sum_node_bottom2 - last_node_bottom2
                          endif
                         
                         
                         ! Recognizing bottom face - end
                         
                      
                         C=dsqrt(mu*rho)
                         ellez=dabs(z_node_bottom1 - z_node_bottom2)
                  
                         fn = 0
                         do ifun = 1,nfunc
                            if (fun_plaY(ipl).eq.tag_func(ifun)) fn = ifun
                         enddo
                         if (fn.gt.0) then
                            
                            !****** Situation S1 or S3 - begin ******
                            
                            if ((sum_node_bottom.eq.10).or.(sum_node_bottom.eq.26)) then
                             
                               
                               !*** Situation S1 ***
                               
                               if (sum_node_bottom.eq.10) then
                               
                                  sit = 1 
                            
                               !*** Situation S3 ***
                               
                               else
                               
                                  sit = 2
                            
                               endif
                  
                               ! Even spectral polynomial degree
                               if (mod(nn-1,2).eq.0) then
                                 k = ((nn-1)/2)+1
                                 
                               ! Odd spectral polynomial degree
                               else
                                 k = (int((nn-1)/2))+num_sit(sit)
                               endif
                                 
                               do i = 1,nn
                                 do j = 1,nn
                                      dxdx = alfa11(ie) +beta12(ie)*ct(k) &
                                           + beta13(ie)*ct(j) +gamma1(ie)*ct(j)*ct(k)
                                      dydx = alfa21(ie) +beta22(ie)*ct(k) &
                                           + beta23(ie)*ct(j) +gamma2(ie)*ct(j)*ct(k)
                                      dzdx = alfa31(ie) +beta32(ie)*ct(k) &
                                           + beta33(ie)*ct(j) +gamma3(ie)*ct(j)*ct(k)
                                      dxdy = alfa12(ie) +beta11(ie)*ct(k) &
                                           + beta13(ie)*ct(i) +gamma1(ie)*ct(k)*ct(i)
                                      dydy = alfa22(ie) +beta21(ie)*ct(k) &
                                           + beta23(ie)*ct(i) +gamma2(ie)*ct(k)*ct(i)
                                      dzdy = alfa32(ie) +beta31(ie)*ct(k) &
                                           + beta33(ie)*ct(i) +gamma3(ie)*ct(k)*ct(i)
                                      dxdz = alfa13(ie) +beta11(ie)*ct(j) &
                                           + beta12(ie)*ct(i) +gamma1(ie)*ct(i)*ct(j)
                                      dydz = alfa23(ie) +beta21(ie)*ct(j) &
                                           + beta22(ie)*ct(i) +gamma2(ie)*ct(i)*ct(j)
                                      dzdz = alfa33(ie) +beta31(ie)*ct(j) &
                                           + beta32(ie)*ct(i) +gamma3(ie)*ct(i)*ct(j)
                                  
                                      det_j = dxdz * (dydx*dzdy - dzdx*dydy) &
                                            - dydz * (dxdx*dzdy - dzdx*dxdy) &
                                            + dzdz * (dxdx*dydy - dydx*dxdy)
                  
                                      is = nn*nn*(k -1) +nn*(j -1) +i
                                      in = cs_loc(cs_loc(ie -1) +is)
                  
                                      term = 4 * C * val_plaY(ipl,1) * ww(i) * ww (j) * dabs(det_j) / ellez
                                      fmat(fn,(3*(in -1) +2)) = fmat(fn,(3*(in -1) +2)) + term
                                 enddo
                             enddo
                               
                               
                              
                               
                            !****** Situation S2 or S4 - begin ******
                            
                            elseif ((sum_node_bottom.eq.14).or.(sum_node_bottom.eq.22)) then
                  
                             !*** Situation S2 ***
                               
                               if (sum_node_bottom.eq.14) then
                               
                                  sit = 1 
                            
                               !*** Situation S4 ***
                               
                               else
                               
                                  sit = 2
                            
                               endif
                  
                               ! Even spectral polynomial degree
                               if (mod(nn-1,2).eq.0) then
                                 j = ((nn-1)/2)+1
                                 
                               ! Odd spectral polynomial degree
                               else
                                 j = (int((nn-1)/2))+num_sit(sit)
                               endif
                                 
                               do i = 1,nn
                                 do k = 1,nn 
                  
                                    dxdx = alfa11(ie) +beta12(ie)*ct(k) &
                                         + beta13(ie)*ct(j) +gamma1(ie)*ct(j)*ct(k)
                                    dydx = alfa21(ie) +beta22(ie)*ct(k) &
                                         + beta23(ie)*ct(j) +gamma2(ie)*ct(j)*ct(k)
                                    dzdx = alfa31(ie) +beta32(ie)*ct(k) &
                                         + beta33(ie)*ct(j) +gamma3(ie)*ct(j)*ct(k)
                                    dxdy = alfa12(ie) +beta11(ie)*ct(k) &
                                         + beta13(ie)*ct(i) +gamma1(ie)*ct(k)*ct(i)
                                    dydy = alfa22(ie) +beta21(ie)*ct(k) &
                                         + beta23(ie)*ct(i) +gamma2(ie)*ct(k)*ct(i)
                                    dzdy = alfa32(ie) +beta31(ie)*ct(k) &
                                         + beta33(ie)*ct(i) +gamma3(ie)*ct(k)*ct(i)
                                    dxdz = alfa13(ie) +beta11(ie)*ct(j) &
                                         + beta12(ie)*ct(i) +gamma1(ie)*ct(i)*ct(j)
                                    dydz = alfa23(ie) +beta21(ie)*ct(j) &
                                         + beta22(ie)*ct(i) +gamma2(ie)*ct(i)*ct(j)
                                    dzdz = alfa33(ie) +beta31(ie)*ct(j) &
                                         + beta32(ie)*ct(i) +gamma3(ie)*ct(i)*ct(j)
                                    det_j = dxdz * (dydx*dzdy - dzdx*dydy) &
                                          - dydz * (dxdx*dzdy - dzdx*dxdy) &
                                          + dzdz * (dxdx*dydy - dydx*dxdy)
                  
                                    is = nn*nn*(k -1) +nn*(j -1) +i
                                    in = cs_loc(cs_loc(ie -1) +is)
                               
                                    term = 4 * C * val_plaY(ipl,1) * ww(i) * ww (k) * dabs(det_j) / ellez
                                    fmat(fn,(3*(in -1) +2)) = fmat(fn,(3*(in -1) +2)) + term
                                 enddo
                               enddo
                              
                              
                              
                           !****** Situation S5 or S6 - begin ******
                            
                            elseif (sum_node_bottom.eq.18) then
                  
                             !*** Situation S5 ***
                               
                               if ((sum_node_bottom_first3.eq.10) &
                                  .or.(sum_node_bottom_first3.eq.13) &
                                  .or.(sum_node_bottom_first3.eq.14)) then
                               
                                  sit = 1 
                            
                               !*** Situation S6 ***
                               
                               else
                               
                                  sit = 2
                            
                               endif
                  
                               ! Even spectral polynomial degree
                               if (mod(nn-1,2).eq.0) then
                                 i = ((nn-1)/2)+1
                                 
                               ! Odd spectral polynomial degree
                               else
                                 i = (int((nn-1)/2))+num_sit(sit)
                               endif
                                 
                               do j = 1,nn
                                 do k = 1,nn
                  
                                    dxdx = alfa11(ie) +beta12(ie)*ct(k) &
                                         + beta13(ie)*ct(j) +gamma1(ie)*ct(j)*ct(k)
                                    dydx = alfa21(ie) +beta22(ie)*ct(k) &
                                         + beta23(ie)*ct(j) +gamma2(ie)*ct(j)*ct(k)
                                    dzdx = alfa31(ie) +beta32(ie)*ct(k) &
                                         + beta33(ie)*ct(j) +gamma3(ie)*ct(j)*ct(k)
                                    dxdy = alfa12(ie) +beta11(ie)*ct(k) &
                                         + beta13(ie)*ct(i) +gamma1(ie)*ct(k)*ct(i)
                                    dydy = alfa22(ie) +beta21(ie)*ct(k) &
                                         + beta23(ie)*ct(i) +gamma2(ie)*ct(k)*ct(i)
                                    dzdy = alfa32(ie) +beta31(ie)*ct(k) &
                                         + beta33(ie)*ct(i) +gamma3(ie)*ct(k)*ct(i)
                                    dxdz = alfa13(ie) +beta11(ie)*ct(j) &
                                         + beta12(ie)*ct(i) +gamma1(ie)*ct(i)*ct(j)
                                    dydz = alfa23(ie) +beta21(ie)*ct(j) &
                                         + beta22(ie)*ct(i) +gamma2(ie)*ct(i)*ct(j)
                                    dzdz = alfa33(ie) +beta31(ie)*ct(j) &
                                         + beta32(ie)*ct(i) +gamma3(ie)*ct(i)*ct(j)
                                    det_j = dxdz * (dydx*dzdy - dzdx*dydy) &
                                          - dydz * (dxdx*dzdy - dzdx*dxdy) &
                                          + dzdz * (dxdx*dydy - dydx*dxdy)
                  
                                    is = nn*nn*(k -1) +nn*(j -1) +i
                                    in = cs_loc(cs_loc(ie -1) +is)
                               
                                    term = 4 * C * val_plaY(ipl,1) * ww(j) * ww (k) * dabs(det_j) / ellez
                                    fmat(fn,(3*(in -1) +2)) = fmat(fn,(3*(in -1) +2)) + term
                  
                                 enddo
                               enddo
                              
                            endif 
                               
                  
                         endif
                      endif  !Check on the Plane Wave Material - end
                  
                   enddo
         endif  ! Plane Wave Y - end
                               
         if (nl_plaZ.gt.0) then       ! Plane Wave Z
                   do ipl = 1,nl_plaZ
                  
                      if (tag_mat(im).eq.tag_plaZ(ipl)) then  !Check on the Plane Wave Material - start
                      
                         ! Recognizing bottom face - begin
                         sum_node_bottom1 = 0
                         sum_node_bottom2 = 0
                         last_node_bottom1 = 0
                         last_node_bottom2 = 0
                         
                         do i = 1,8
                             ie_glob = local_el_num(ie)
                             call GET_INDLOC_FROM_INDGLO(local_n_num, nnod_loc, con_hexa(ie_glob,2), ic1)
                             call GET_INDLOC_FROM_INDGLO(local_n_num, nnod_loc, con_hexa(ie_glob,i+1), ic2)
                             
                             if(ic1 .ne. 0 .and. ic2 .ne. 0) then                           
                                 if (dabs(zs_loc(ic1) - zs_loc(ic2)) .le. 1.d-4) then
                                     last_node_bottom1 = i
                                     sum_node_bottom1 = sum_node_bottom1 + i
                                     z_node_bottom1 = zs_loc(ic1)
                                 else
                                     last_node_bottom2 = i
                                     sum_node_bottom2 = sum_node_bottom2 + i
                                     z_node_bottom2 = zs_loc(ic2)
                                 endif     
                             endif
                          enddo
                          
                          if (z_node_bottom1.lt.z_node_bottom2) then
                                 sum_node_bottom = sum_node_bottom1
                                 sum_node_bottom_first3 = sum_node_bottom1 - last_node_bottom1
                          else
                                 sum_node_bottom = sum_node_bottom2
                                 sum_node_bottom_first3 = sum_node_bottom2 - last_node_bottom2
                          endif
                         
                         
                         ! Recognizing bottom face - end
                         
                      
                         C=dsqrt((lambda+2*mu)*rho)
                         ellez=dabs(z_node_bottom1 - z_node_bottom2)
                  
                         fn = 0
                         do ifun = 1,nfunc
                            if (fun_plaZ(ipl).eq.tag_func(ifun)) fn = ifun
                         enddo
                         if (fn.gt.0) then
                            
                            !****** Situation S1 or S3 - begin ******
                            
                            if ((sum_node_bottom.eq.10).or.(sum_node_bottom.eq.26)) then
                             
                               
                               !*** Situation S1 ***
                               
                               if (sum_node_bottom.eq.10) then
                               
                                  sit = 1 
                            
                               !*** Situation S3 ***
                               
                               else
                               
                                  sit = 2
                            
                               endif
                  
                               ! Even spectral polynomial degree
                               if (mod(nn-1,2).eq.0) then
                                 k = ((nn-1)/2)+1
                                 
                               ! Odd spectral polynomial degree
                               else
                                 k = (int((nn-1)/2))+num_sit(sit)
                               endif
                                 
                               do i = 1,nn
                                 do j = 1,nn
                                        dxdx = alfa11(ie) +beta12(ie)*ct(k) &
                                             + beta13(ie)*ct(j) +gamma1(ie)*ct(j)*ct(k)
                                        dydx = alfa21(ie) +beta22(ie)*ct(k) &
                                             + beta23(ie)*ct(j) +gamma2(ie)*ct(j)*ct(k)
                                        dzdx = alfa31(ie) +beta32(ie)*ct(k) &
                                             + beta33(ie)*ct(j) +gamma3(ie)*ct(j)*ct(k)
                                        dxdy = alfa12(ie) +beta11(ie)*ct(k) &
                                             + beta13(ie)*ct(i) +gamma1(ie)*ct(k)*ct(i)
                                        dydy = alfa22(ie) +beta21(ie)*ct(k) &
                                             + beta23(ie)*ct(i) +gamma2(ie)*ct(k)*ct(i)
                                        dzdy = alfa32(ie) +beta31(ie)*ct(k) &
                                             + beta33(ie)*ct(i) +gamma3(ie)*ct(k)*ct(i)
                                        dxdz = alfa13(ie) +beta11(ie)*ct(j) &
                                             + beta12(ie)*ct(i) +gamma1(ie)*ct(i)*ct(j)
                                        dydz = alfa23(ie) +beta21(ie)*ct(j) &
                                             + beta22(ie)*ct(i) +gamma2(ie)*ct(i)*ct(j)
                                        dzdz = alfa33(ie) +beta31(ie)*ct(j) &
                                             + beta32(ie)*ct(i) +gamma3(ie)*ct(i)*ct(j)
                                        det_j = dxdz * (dydx*dzdy - dzdx*dydy) &
                                              - dydz * (dxdx*dzdy - dzdx*dxdy) &
                                              + dzdz * (dxdx*dydy - dydx*dxdy)
                  
                                         is = nn*nn*(k -1) +nn*(j -1) +i
                                         in = cs_loc(cs_loc(ie -1) +is)
                               
                                         term = 4 * C * val_plaZ(ipl,1)  * ww(i) * ww (j) * dabs(det_j) / ellez
                                         fmat(fn,(3*(in -1) +3)) = fmat(fn,(3*(in -1) +3)) + term
                                 enddo
                               enddo
                               
                               
                              
                               
                            !****** Situation S2 or S4 - begin ******
                            
                            elseif ((sum_node_bottom.eq.14).or.(sum_node_bottom.eq.22)) then
                  
                             !*** Situation S2 ***
                               
                               if (sum_node_bottom.eq.14) then
                               
                                  sit = 1 
                            
                               !*** Situation S4 ***
                               
                               else
                               
                                  sit = 2
                            
                               endif
                  
                               ! Even spectral polynomial degree
                               if (mod(nn-1,2).eq.0) then
                                 j = ((nn-1)/2)+1
                                 
                               ! Odd spectral polynomial degree
                               else
                                 j = (int((nn-1)/2))+num_sit(sit)
                               endif
                                 
                               do i = 1,nn
                                 do k = 1,nn
                                        dxdx = alfa11(ie) +beta12(ie)*ct(k) &
                                             + beta13(ie)*ct(j) +gamma1(ie)*ct(j)*ct(k)
                                        dydx = alfa21(ie) +beta22(ie)*ct(k) &
                                             + beta23(ie)*ct(j) +gamma2(ie)*ct(j)*ct(k)
                                        dzdx = alfa31(ie) +beta32(ie)*ct(k) &
                                             + beta33(ie)*ct(j) +gamma3(ie)*ct(j)*ct(k)
                                 
                                        dxdy = alfa12(ie) +beta11(ie)*ct(k) &
                                             + beta13(ie)*ct(i) +gamma1(ie)*ct(k)*ct(i)
                                        dydy = alfa22(ie) +beta21(ie)*ct(k) &
                                             + beta23(ie)*ct(i) +gamma2(ie)*ct(k)*ct(i)
                                        dzdy = alfa32(ie) +beta31(ie)*ct(k) &
                                             + beta33(ie)*ct(i) +gamma3(ie)*ct(k)*ct(i)
                                 
                                        dxdz = alfa13(ie) +beta11(ie)*ct(j) &
                                             + beta12(ie)*ct(i) +gamma1(ie)*ct(i)*ct(j)
                                        dydz = alfa23(ie) +beta21(ie)*ct(j) &
                                             + beta22(ie)*ct(i) +gamma2(ie)*ct(i)*ct(j)
                                        dzdz = alfa33(ie) +beta31(ie)*ct(j) &
                                             + beta32(ie)*ct(i) +gamma3(ie)*ct(i)*ct(j)
                                 
                                        det_j = dxdz * (dydx*dzdy - dzdx*dydy) &
                                              - dydz * (dxdx*dzdy - dzdx*dxdy) &
                                              + dzdz * (dxdx*dydy - dydx*dxdy)
                  
                                         is = nn*nn*(k -1) +nn*(j -1) +i
                                         in = cs_loc(cs_loc(ie -1) +is)
                  
                                         term = 4 * C * val_plaZ(ipl,1) * ww(i) * ww (k) * dabs(det_j) / ellez
                                         fmat(fn,(3*(in -1) +3)) = fmat(fn,(3*(in -1) +3)) + term
                                 enddo
                               enddo
                              
                              
                              
                           !****** Situation S5 or S6 - begin ******
                            
                            elseif (sum_node_bottom.eq.18) then
                  
                             !*** Situation S5 ***
                               
                               if ((sum_node_bottom_first3.eq.10) &
                                  .or.(sum_node_bottom_first3.eq.13) &
                                  .or.(sum_node_bottom_first3.eq.14)) then
                               
                                  sit = 1 
                            
                               !*** Situation S6 ***
                               
                               else
                               
                                  sit = 2
                            
                               endif
                  
                               ! Even spectral polynomial degree
                               if (mod(nn-1,2).eq.0) then
                                 i = ((nn-1)/2)+1
                                 
                               ! Odd spectral polynomial degree
                               else
                                 i = (int((nn-1)/2))+num_sit(sit)
                               endif
                                 
                               do j = 1,nn
                                 do k = 1,nn 
                  
                                     dxdx = alfa11(ie) +beta12(ie)*ct(k) &
                                          + beta13(ie)*ct(j) +gamma1(ie)*ct(j)*ct(k)
                                     dydx = alfa21(ie) +beta22(ie)*ct(k) &
                                          + beta23(ie)*ct(j) +gamma2(ie)*ct(j)*ct(k)
                                     dzdx = alfa31(ie) +beta32(ie)*ct(k) &
                                          + beta33(ie)*ct(j) +gamma3(ie)*ct(j)*ct(k)
                                 
                                     dxdy = alfa12(ie) +beta11(ie)*ct(k) &
                                          + beta13(ie)*ct(i) +gamma1(ie)*ct(k)*ct(i)
                                     dydy = alfa22(ie) +beta21(ie)*ct(k) &
                                          + beta23(ie)*ct(i) +gamma2(ie)*ct(k)*ct(i)
                                     dzdy = alfa32(ie) +beta31(ie)*ct(k) &
                                          + beta33(ie)*ct(i) +gamma3(ie)*ct(k)*ct(i)
                                  
                                     dxdz = alfa13(ie) +beta11(ie)*ct(j) &
                                          + beta12(ie)*ct(i) +gamma1(ie)*ct(i)*ct(j)
                                     dydz = alfa23(ie) +beta21(ie)*ct(j) &
                                          + beta22(ie)*ct(i) +gamma2(ie)*ct(i)*ct(j)
                                     dzdz = alfa33(ie) +beta31(ie)*ct(j) &
                                          + beta32(ie)*ct(i) +gamma3(ie)*ct(i)*ct(j)
                                  
                                     det_j = dxdz * (dydx*dzdy - dzdx*dydy) &
                                           - dydz * (dxdx*dzdy - dzdx*dxdy) &
                                           + dzdz * (dxdx*dydy - dydx*dxdy)
                  
                                         is = nn*nn*(k -1) +nn*(j -1) +i
                                         in = cs_loc(cs_loc(ie -1) +is)
                               
                                         term = 4 * C * val_plaZ(ipl,1) * ww(j) * ww (k) * dabs(det_j) / ellez
                                         fmat(fn,(3*(in -1) +3)) = fmat(fn,(3*(in -1) +3)) + term
                                 enddo
                               enddo
                              
                            endif 
                               
                  
                         endif
                      endif  !Check on the Plane Wave Material - end
                  
                   enddo
         endif  ! Plane Wave Z - end


          !############################################################################
          !##################              PLANE WAVE END                ##############
          !############################################################################


          !############################################################################
          !##################             SEISMIC MOMENT  BEGIN          ##############
          !############################################################################


         if (nl_sism.gt.0) then
                  
                   do isism = 1, nl_sism
                      do ipsism = 1, num_ns(isism)
                                                            
                      fn = 0
                                                                                      
                      do ifun = 1,nfunc
                         if (fun_sism(isism).eq.tag_func(ifun)) fn = ifun                        
                      enddo                                                                
                      
                      if (fn.gt.0) then                                                        
                         do k = 1,nn                                                        
                            do j = 1,nn                                                        
                               do i = 1,nn                                                
                                  dxdx = alfa11(ie) +beta12(ie)*ct(k) &                        
                                       + beta13(ie)*ct(j) +gamma1(ie)*ct(j)*ct(k)        
                                  dydx = alfa21(ie) +beta22(ie)*ct(k) &                        
                                       + beta23(ie)*ct(j) +gamma2(ie)*ct(j)*ct(k)        
                                  dzdx = alfa31(ie) +beta32(ie)*ct(k) &                        
                                       + beta33(ie)*ct(j) +gamma3(ie)*ct(j)*ct(k)        
                                  
                                  dxdy = alfa12(ie) +beta11(ie)*ct(k) &                        
                                       + beta13(ie)*ct(i) +gamma1(ie)*ct(k)*ct(i)        
                                  dydy = alfa22(ie) +beta21(ie)*ct(k) &                        
                                       + beta23(ie)*ct(i) +gamma2(ie)*ct(k)*ct(i)        
                                  dzdy = alfa32(ie) +beta31(ie)*ct(k) &                        
                                       + beta33(ie)*ct(i) +gamma3(ie)*ct(k)*ct(i)        
                                  
                                  dxdz = alfa13(ie) +beta11(ie)*ct(j) &                        
                                       + beta12(ie)*ct(i) +gamma1(ie)*ct(i)*ct(j)        
                                  dydz = alfa23(ie) +beta21(ie)*ct(j) &                        
                                       + beta22(ie)*ct(i) +gamma2(ie)*ct(i)*ct(j)        
                                  dzdz = alfa33(ie) +beta31(ie)*ct(j) &                        
                                       + beta32(ie)*ct(i) +gamma3(ie)*ct(i)*ct(j)        
                                  
                                  det_j = dxdz * (dydx*dzdy - dzdx*dydy) &                
                                        - dydz * (dxdx*dzdy - dzdx*dxdy) &                
                                        + dzdz * (dxdx*dydy - dydx*dxdy)                        
                                  
                                  is = nn*nn*(k -1) +nn*(j -1) +i                        
                                  in = cs_loc(cs_loc(ie -1) +is)                                        
                  
                                     
                                     if (local_n_num(in) .eq. sour_ns(ipsism,isism)) then                           
                                        facsmom(isism,1) = facsmom(isism,1) + det_j * ww(i) * ww(j) * ww(k)        
                                        facsmom(isism,2) = facsmom(isism,2) + det_j * ww(i) * ww(j) * ww(k)        
                                        facsmom(isism,3) = facsmom(isism,3) + det_j * ww(i) * ww(j) * ww(k)        
                                        facsmom(isism,4) = facsmom(isism,4) + det_j * ww(i) * ww(j) * ww(k)        
                                        facsmom(isism,5) = facsmom(isism,5) + det_j * ww(i) * ww(j) * ww(k)        
                                        facsmom(isism,6) = facsmom(isism,6) + det_j * ww(i) * ww(j) * ww(k)         
                  
                                        length_cns = length_cns + 1                        
                  
                                     endif                                                
                               enddo !i                                                        
                            enddo !j                                                        
                         enddo !k                                                        
                      endif !if (fn.gt.0) then                                                
                                              
                   enddo !ip                                                                
                  
                   enddo !isism                                                                
                  
         endif !if (nl_sism.gt.0) then                                                
                               
                                                    
          !############################################################################
          !##################             SEISMIC MOMENT END             ##############
          !############################################################################

       
          !############################################################################
          !##################          EXPLOSIVE SOURCE BEGIN            ##############
          !############################################################################

         if (nl_expl.gt.0) then
                   do iexpl = 1,nl_expl
                      do ipexpl = 1,num_ne(iexpl)
                                                            
                      fn = 0                                                                
                      do ifun = 1,nfunc                                                        
                         if (fun_expl(iexpl).eq.tag_func(ifun)) fn = ifun                        
                      enddo                                                                
                      
                      if (fn.gt.0) then                                                        
                         do k = 1,nn                                                        
                            do j = 1,nn                                                        
                               do i = 1,nn                                                
                                  dxdx = alfa11(ie) +beta12(ie)*ct(k) &                        
                                       + beta13(ie)*ct(j) +gamma1(ie)*ct(j)*ct(k)        
                                  dydx = alfa21(ie) +beta22(ie)*ct(k) &                        
                                       + beta23(ie)*ct(j) +gamma2(ie)*ct(j)*ct(k)        
                                  dzdx = alfa31(ie) +beta32(ie)*ct(k) &                        
                                       + beta33(ie)*ct(j) +gamma3(ie)*ct(j)*ct(k)        
                                  
                                  dxdy = alfa12(ie) +beta11(ie)*ct(k) &                        
                                       + beta13(ie)*ct(i) +gamma1(ie)*ct(k)*ct(i)        
                                  dydy = alfa22(ie) +beta21(ie)*ct(k) &                        
                                       + beta23(ie)*ct(i) +gamma2(ie)*ct(k)*ct(i)        
                                  dzdy = alfa32(ie) +beta31(ie)*ct(k) &                        
                                       + beta33(ie)*ct(i) +gamma3(ie)*ct(k)*ct(i)        
                                  
                                  dxdz = alfa13(ie) +beta11(ie)*ct(j) &                        
                                       + beta12(ie)*ct(i) +gamma1(ie)*ct(i)*ct(j)        
                                  dydz = alfa23(ie) +beta21(ie)*ct(j) &                        
                                       + beta22(ie)*ct(i) +gamma2(ie)*ct(i)*ct(j)        
                                  dzdz = alfa33(ie) +beta31(ie)*ct(j) &                        
                                       + beta32(ie)*ct(i) +gamma3(ie)*ct(i)*ct(j)        
                                  
                                  det_j = dxdz * (dydx*dzdy - dzdx*dydy) &                
                                        - dydz * (dxdx*dzdy - dzdx*dxdy) &                
                                        + dzdz * (dxdx*dydy - dydx*dxdy)                        
                                  
                                  is = nn*nn*(k -1) +nn*(j -1) +i                        
                                  in = cs_loc(cs_loc(ie -1) +is)                                        
                                  
                                     if (local_n_num(in) .eq. sour_ne(ipexpl,iexpl)) then                                        
                                        facsexpl(iexpl,1) = facsexpl(iexpl,1) + det_j * ww(i) * ww(j) * ww(k)        
                                        facsexpl(iexpl,2) = facsexpl(iexpl,2) + det_j * ww(i) * ww(j) * ww(k)        
                                        facsexpl(iexpl,3) = facsexpl(iexpl,3) + det_j * ww(i) * ww(j) * ww(k)        
                                        facsexpl(iexpl,4) = facsexpl(iexpl,4) + det_j * ww(i) * ww(j) * ww(k)        
                                        facsexpl(iexpl,5) = facsexpl(iexpl,5) + det_j * ww(i) * ww(j) * ww(k)        
                                        facsexpl(iexpl,6) = facsexpl(iexpl,6) + det_j * ww(i) * ww(j) * ww(k)         
                  
                                        length_cne = length_cne + 1                        
                  
                                     endif                                                
                               enddo !i                                                        
                            enddo !j                                                        
                         enddo !k                                                        
                      endif !if (fn.gt.0) then                                                
                                              
                   enddo !ip                                                                
                  
                   enddo !iexpl                                                                
                  
                  
         endif !if (nl_expl.gt.0) then                                                

                                       
                                       
          !############################################################################
          !##################          EXPLOSIVE SOURCE END              ##############
          !############################################################################
         deallocate(ct,dd,ww)
      enddo  




!############################################################################
!##################          SEISMIC MOMENT BEGIN              ##############
!############################################################################
      if (nl_sism.gt.0) then                                                                        
            do isism = 1,nl_sism                                                                
               call MPI_BARRIER(mpi_comm,mpi_ierr)
               call MPI_ALLREDUCE(facsmom(isism,:),sum_facs, 6, SPEED_DOUBLE, MPI_SUM, mpi_comm, mpi_ierr)
               facsmom(isism,:) = sum_facs                    
            enddo                        
         
      endif        

      if (nl_sism.gt.0) then                                                                        
                 if (srcmodflag.eq.0) then
                    val_st = 12
                  elseif (srcmodflag.eq.1) then
                    val_st = 6
                  endif
         do isism = 1,nl_sism                                                                
                     slip1_sism = val_sism(isism,val_st + 1)                                                
                     slip2_sism = val_sism(isism,val_st + 2)                                                
                     slip3_sism = val_sism(isism,val_st + 3)                                                
                     
                     norm1_sism = val_sism(isism,val_st + 4)                                                
                     norm2_sism = val_sism(isism,val_st + 5)                                                
                     norm3_sism = val_sism(isism,val_st + 6)                                                
                     
                     amp_sism = val_sism(isism,val_st + 8)
                                         
                     tau_sism = val_sism(isism,val_st + 9)                                                        

                    ! Seismic Tensor is given is terms of norma and slip vector:
                                        ! nx, ny, nz and sx, sy, sz
                                        ! i.e.: SULMONA Case
                    if(facsmom(isism,1) .ne. 0) &
                                       facsmom(isism,1) = 1/facsmom(isism,1) &                                        
                                       * (slip1_sism*norm1_sism+slip1_sism*norm1_sism) &        
                                       * amp_sism                                                
                    if(facsmom(isism,2) .ne. 0) &
                                       facsmom(isism,2) = 1/facsmom(isism,2) &                                        
                                       * (slip2_sism*norm2_sism+slip2_sism*norm2_sism) &        
                                       * amp_sism                                                
                    if(facsmom(isism,3) .ne. 0) &
                                       facsmom(isism,3) = 1/facsmom(isism,3) &                                        
                                       * (slip3_sism*norm3_sism+slip3_sism*norm3_sism) &        
                                       * amp_sism                                                
                    if(facsmom(isism,4) .ne. 0) &
                                       facsmom(isism,4) = 1/facsmom(isism,4) &                                        
                                       * (slip2_sism*norm3_sism+slip3_sism*norm2_sism) &        
                                       * amp_sism                                                
                    if(facsmom(isism,5) .ne. 0) &
                                       facsmom(isism,5) = 1/facsmom(isism,5) &                                        
                                       * (slip1_sism*norm3_sism+slip3_sism*norm1_sism) &        
                                       * amp_sism                                                
                    if(facsmom(isism,6) .ne. 0) &
                                       facsmom(isism,6) = 1/facsmom(isism,6) &                                        
                                       * (slip1_sism*norm2_sism+slip2_sism*norm1_sism) &        
                                       * amp_sism                          
                                                                           
                     tausmom(isism,1) = tau_sism                        

                    ! Seismic Tensor is given is terms of norma and slip vector 
                                        ! mxx,myy,mzz,mxy,mxz,myz
                                        !  i.e.: CASHIMA (ONLY)
                    !facsmom(isism,1) = 1/facsmom(isism,1) &                                        
                    !                   * slip1_sism                                                
                    !facsmom(isism,2) = 1/facsmom(isism,2) &                                        
                    !                   * slip2_sism                                                
                    !facsmom(isism,3) = 1/facsmom(isism,3) &                                        
                    !                   * slip3_sism                                        
                    !facsmom(isism,4) = 1/facsmom(isism,4) &                                        
                    !                   * norm1_sism                                        
                    !facsmom(isism,5) = 1/facsmom(isism,5) &                                        
                    !                   * norm2_sism                                                        
                    !facsmom(isism,6) = 1/facsmom(isism,6) &                                        
                    !                   * norm3_sism                                                                                     
                    !facsmom(isism,1) = 1/facsmom(isism,1) &                                        
                    !                   * 1 &        
                    !                   * amp_sism                                                
                    !facsmom(isism,2) = 1/facsmom(isism,2) &                                        
                    !                   * 1 &        
                    !                   * amp_sism                                                
                    !facsmom(isism,3) = 1/facsmom(isism,3) &                                        
                    !                   * 1 &        
                    !                   * amp_sism                                                
                    !facsmom(isism,4) = 0.0d0                                                
                    !facsmom(isism,5) = 0.0d0                                                
                    !facsmom(isism,6) = 0.0d0                                                  
         enddo                        
              
      endif        
          

!############################################################################
!##################          SEISMIC MOMENT END                ##############
!############################################################################


!############################################################################
!##################         EXPLOSIVE SOURCE BEGIN             ##############
!############################################################################

      if (nl_expl.gt.0) then                                                                        
         do iexpl = 1,nl_expl                                                                

            call MPI_BARRIER(mpi_comm,mpi_ierr)
            call MPI_ALLREDUCE(facsexpl(iexpl,:),sum_facs, 6, SPEED_DOUBLE, MPI_SUM, mpi_comm, mpi_ierr)
            facsexpl(iexpl,:) = sum_facs                    
         enddo
      endif      



      if (nl_expl.gt.0) then                                                                        
             
         do iexpl = 1,nl_expl                                                                
            slip1_expl = val_expl(iexpl,13)                                                
            slip2_expl = val_expl(iexpl,14)                                                
            slip3_expl = val_expl(iexpl,15)                                                
             
            norm1_expl = val_expl(iexpl,16)                                                
            norm2_expl = val_expl(iexpl,17)                                                
            norm3_expl = val_expl(iexpl,18)                                                
             
            amp_expl = val_expl(iexpl,20)                                                
         
            facsexpl(iexpl,1) = 1/facsexpl(iexpl,1) &                                        
                               * slip1_expl &        
                               * amp_expl                                                
            facsexpl(iexpl,2) = 1/facsexpl(iexpl,2) &                                        
                               * slip2_expl &        
                               * amp_expl                                                
            facsexpl(iexpl,3) = 1/facsexpl(iexpl,3) &                                        
                               * slip3_expl &        
                               * amp_expl                                                
            facsexpl(iexpl,4) = 1/facsexpl(iexpl,4) &                                        
                               * norm1_expl &        
                               * amp_expl                                                
            facsexpl(iexpl,5) = 1/facsexpl(iexpl,5) &                                        
                               * norm2_expl &        
                               * amp_expl                                                
            facsexpl(iexpl,6) = 1/facsexpl(iexpl,6) &                                        
                               * norm3_expl &        
                               * amp_expl                                         
         enddo                                                                                
      endif        

!############################################################################
!##################         EXPLOSIVE SOURCE END               ##############
!############################################################################

      ne_loc = cs_loc(0) -1

      do ie = 1,ne_loc
         im = cs_loc(cs_loc(ie -1) +0);
         nn = sdeg_mat(im) +1
   
         allocate(ct(nn),ww(nn),dd(nn,nn))
         call MAKE_LGL_NW(nn,ct,ww,dd)
         rho = prop_mat(im,1) 
         lambda = prop_mat(im,2)
         mu = prop_mat(im,3)

         if (n_test.gt.0) then    ! Test load

            do ip = 1, n_test

               fn = 0
               do ifun = 1,nfunc
                  if (fun_test(ip).eq.tag_func(ifun)) fn = ifun
               enddo                  
                    
               if (fn.gt.0) then
   
                 rho = prop_mat(im,1); lambda = prop_mat(im,2); mu = prop_mat(im,3)

                         
                  do k = 1,nn
                     do j = 1,nn
                        do i = 1,nn
                           dxdx = alfa11(ie) +beta12(ie)*ct(k) &
                            + beta13(ie)*ct(j) +gamma1(ie)*ct(j)*ct(k)
                          dydx = alfa21(ie) +beta22(ie)*ct(k) &
                              + beta23(ie)*ct(j) +gamma2(ie)*ct(j)*ct(k)
                           dzdx = alfa31(ie) +beta32(ie)*ct(k) &
                              + beta33(ie)*ct(j) +gamma3(ie)*ct(j)*ct(k)
                
                          dxdy = alfa12(ie) +beta11(ie)*ct(k) &
                              + beta13(ie)*ct(i) +gamma1(ie)*ct(k)*ct(i)
                           dydy = alfa22(ie) +beta21(ie)*ct(k) &
                              + beta23(ie)*ct(i) +gamma2(ie)*ct(k)*ct(i)
                           dzdy = alfa32(ie) +beta31(ie)*ct(k) &
                              + beta33(ie)*ct(i) +gamma3(ie)*ct(k)*ct(i)
                
                           dxdz = alfa13(ie) +beta11(ie)*ct(j) &
                              + beta12(ie)*ct(i) +gamma1(ie)*ct(i)*ct(j)
                           dydz = alfa23(ie) +beta21(ie)*ct(j) &
                              + beta22(ie)*ct(i) +gamma2(ie)*ct(i)*ct(j)
                           dzdz = alfa33(ie) +beta31(ie)*ct(j) &
                              + beta32(ie)*ct(i) +gamma3(ie)*ct(i)*ct(j)
                
                           det_j = dxdz * (dydx*dzdy - dzdx*dydy) &
                                 - dydz * (dxdx*dzdy - dzdx*dxdy) &
                                 + dzdz * (dxdx*dydy - dydx*dxdy)
                
                           is = nn*nn*(k -1) +nn*(j -1) +i
                           in = cs_loc(cs_loc(ie -1) +is)

               
                            x =  xs_loc(in); y = ys_loc(in); z = zs_loc(in); 

                           !EXTERNAL LOAD PAPER WITH BLANCA
                           fmat(fn,(3*(in -1) +1)) = fmat(fn,(3*(in -1) +1)) + det_j * ww(i) * ww(j) * ww(k) * ( &
                                 4.d0*pi**2.d0*dcos(pi*y)*dcos(pi*z)*dsin(pi*y)*dsin(pi*z) * &
                                 (2.d0*lambda - 8.d0*mu + 9.d0*rho - 4.d0*lambda*dcos(pi*x)**2.d0 &
                                 + 8.d0*mu*dcos(pi*x)**2.d0 - 9.d0*rho*dcos(pi*x)**2.d0) )
                    
                           fmat(fn,(3*(in -1) +2)) = fmat(fn,(3*(in -1) +2)) + det_j * ww(i) * ww(j) * ww(k) * ( &
                                 4.d0*pi**2.d0*dcos(pi*x)*dcos(pi*z)*dsin(pi*x)*dsin(pi*z) * &
                                 (2.d0*lambda + 12.d0*mu - 9.d0*rho - 4.d0*lambda*dcos(pi*y)**2.d0 &
                                 - 16.d0*mu*dcos(pi*y)**2.d0 + 9.d0*rho*dcos(pi*y)**2.d0) )

                           fmat(fn,(3*(in -1) +3)) = fmat(fn,(3*(in -1) +3)) + det_j * ww(i) * ww(j) * ww(k) * (&
                                 4.d0*pi**2.d0*dcos(pi*x)*dcos(pi*y)*dsin(pi*x)*dsin(pi*y) * & 
                                 (2.d0*lambda + 12.d0*mu - 9.d0*rho - 4.d0*lambda*dcos(pi*z)**2.d0 &
                                 - 16.d0*mu*dcos(pi*z)**2.d0 + 9.d0*rho*dcos(pi*z)**2.d0) ) 


                        enddo            
                     enddo
                  enddo
               endif
            enddo                    
         endif
         deallocate(ct,ww,dd)
      enddo


if (cs_nnz_bc_loc.gt.0) then
nface_loc = cs_bc_loc(0) -1

do ie = 1,nface_loc
!if ((cs_bc_loc(ie) - cs_bc_loc(ie -1) -1) .ne. nn*nn) then
!deallocate(ct,ww,dd)
nn = int(sqrt(real(cs_bc_loc(ie) - cs_bc_loc(ie -1) -1)))

allocate(ct(nn),ww(nn),dd(nn,nn))
call MAKE_LGL_NW(nn,ct,ww,dd)
!endif

ic1 = cs_bc_loc(cs_bc_loc(ie -1) +nn*nn)
ic2 = cs_bc_loc(cs_bc_loc(ie -1) +1)
ic3 = cs_bc_loc(cs_bc_loc(ie -1) +nn*(nn -1) +1)
ic4 = cs_bc_loc(cs_bc_loc(ie -1) +nn)

l1x = 0.d0
l2x = 0.d0
l1y = 0.d0
l2y = 0.d0 
l1z = 0.d0
l2z = 0.d0
area = 0.d0

if (ic1 .ne. 0 .and. ic2 .ne. 0 .and. ic3 .ne. 0 .and. ic4 .ne. 0) then

l1x = xs_loc(ic1) - xs_loc(ic2)
l1y = ys_loc(ic1) - ys_loc(ic2)
l1z = zs_loc(ic1) - zs_loc(ic2)

l2x = xs_loc(ic3) - xs_loc(ic4)
l2y = ys_loc(ic3) - ys_loc(ic4)
l2z = zs_loc(ic3) - zs_loc(ic4)

area = l1y*l2z -l1z*l2y +l1z*l2x -l1x*l2z +l1x*l2y -l1y*l2x          
area = dabs(0.5d0*area)

endif


if (nl_neuX.gt.0) then   ! Neumann load X
do il = 1,nl_neuX
 if (tag_neuX(il).eq.cs_bc_loc(cs_bc_loc(ie -1) +0)) then
    fn = 0
    do ifun = 1,nfunc
       if (fun_neuX(il).eq.tag_func(ifun)) fn = ifun
    enddo
    
    if (fn.gt.0) then
    
       v1 = val_neuX(il,1)
       v2 = val_neuX(il,2)
       v3 = val_neuX(il,3)
       v4 = val_neuX(il,4)
       
       do j = 1,nn
          do i = 1,nn
             is = nn*(j -1) +i
             in = cs_bc_loc(cs_bc_loc(ie -1) +is)
                v = 0.25*(1.0 -ct(i))*(1.0 -ct(j))*v1 &
                  + 0.25*(1.0 +ct(i))*(1.0 -ct(j))*v2 &
                  + 0.25*(1.0 +ct(i))*(1.0 +ct(j))*v3 &
                  + 0.25*(1.0 -ct(i))*(1.0 +ct(j))*v4
                
                term = 0.25*area*v * ww(i)*ww(j)
                                            
                fmat(fn,(3*(in -1) +1)) = fmat(fn,(3*(in -1) +1)) + term
          enddo
       enddo
    endif
 endif
enddo
endif


if (nl_neuY.gt.0) then   ! Neumann load Y
do il = 1,nl_neuY
 if (tag_neuY(il).eq.cs_bc_loc(cs_bc_loc(ie -1) +0)) then
    fn = 0
    do ifun = 1,nfunc
       if (fun_neuY(il).eq.tag_func(ifun)) fn = ifun
    enddo
    
    if (fn.gt.0) then
       v1 = val_neuY(il,1)
       v2 = val_neuY(il,2)
       v3 = val_neuY(il,3)
       v4 = val_neuY(il,4)
       
       do j = 1,nn
          do i = 1,nn
             is = nn*(j -1) +i
             in = cs_bc_loc(cs_bc_loc(ie -1) +is)
                v = 0.25*(1.0 -ct(i))*(1.0 -ct(j))*v1 &
                  + 0.25*(1.0 +ct(i))*(1.0 -ct(j))*v2 &
                  + 0.25*(1.0 +ct(i))*(1.0 +ct(j))*v3 &
                  + 0.25*(1.0 -ct(i))*(1.0 +ct(j))*v4
                
                term = 0.25*area*v * ww(i)*ww(j)
                
                fmat(fn,(3*(in -1) +2)) = fmat(fn,(3*(in -1) +2)) + term
          enddo
       enddo
    endif
 endif
enddo
endif


if (nl_neuZ.gt.0) then   ! Neumann load Z
do il = 1,nl_neuZ
 if (tag_neuZ(il).eq.cs_bc_loc(cs_bc_loc(ie -1) +0)) then
    fn = 0
    do ifun = 1,nfunc
       if (fun_neuZ(il).eq.tag_func(ifun)) fn = ifun
    enddo
    
    if (fn.gt.0) then
       v1 = val_neuZ(il,1)
       v2 = val_neuZ(il,2)
       v3 = val_neuZ(il,3)
       v4 = val_neuZ(il,4)
       
       do j = 1,nn
          do i = 1,nn
             is = nn*(j -1) +i
             in = cs_bc_loc(cs_bc_loc(ie -1) +is)
                v = 0.25*(1.0 -ct(i))*(1.0 -ct(j))*v1 &
                  + 0.25*(1.0 +ct(i))*(1.0 -ct(j))*v2 &
                  + 0.25*(1.0 +ct(i))*(1.0 +ct(j))*v3 &
                  + 0.25*(1.0 -ct(i))*(1.0 +ct(j))*v4
                
                term = 0.25*area*v * ww(i)*ww(j)
                
                fmat(fn,(3*(in -1) +3)) = fmat(fn,(3*(in -1) +3)) + term
          enddo
       enddo
    endif
 endif
enddo
endif

! ------------------------------------------------------------------------

if (nl_neuN.gt.0) then   
do il = 1,nl_neuN   
 if (tag_neuN(il).eq.cs_bc_loc(cs_bc_loc(ie -1) +0)) then   
    fn = 0   
    do ifun = 1,nfunc  
       if (fun_neuN(il).eq.tag_func(ifun)) fn = ifun   
    enddo   
    
    if (fn.gt.0) then   
      
       rho = prop_mat(1,1); lambda = prop_mat(1,2); mu = prop_mat(1,3)
       
       k1 = val_neuN(il,1)
       k2 = val_neuN(il,2)
       k3 = val_neuN(il,3)

       omega = val_neuN(il,4)*2*PI

       vp = dsqrt((lambda+2*mu)/rho)
       k3 = k3*omega/vp;
       
       
       v1 = 0; !val_neuN(il,1)  
       v2 = 0 !val_neuN(il,2)  
       v3 = 0 !val_neuN(il,3)   

       
       do j = 1,nn  
          do i = 1,nn   
             is = nn*(j -1) +i  
             in = cs_bc_loc(cs_bc_loc(ie -1) +is)
             
             x =  xs_loc(in); y = ys_loc(in); z = zs_loc(in); 

             do index = 1,nelem_neuN 
                     if (i4normal(index) .eq. ie) index_vector = index
             enddo
             
             n1 = normal_nx_el_neuN(index_vector)
             n2 = normal_ny_el_neuN(index_vector)
             n3 = normal_nz_el_neuN(index_vector)
             
             !write(*,*) area
             !read(*,*)
             
             if (tag_func(fn) .eq. 202) then 
                v1 = k3*lambda*n1*cos(k1*x + k2*y + k3*z) + k1*mu*n3*cos(k1*x + k2*y + k3*z)
                v2 = k3*lambda*n2*cos(k1*x + k2*y + k3*z) + k2*mu*n3*cos(k1*x + k2*y + k3*z)
                v3 = n3*(k3*lambda*cos(k1*x + k2*y + k3*z) + 2*k3*mu*cos(k1*x + k2*y + k3*z)) &
                   + k1*mu*n1*cos(k1*x + k2*y + k3*z) + k2*mu*n2*cos(k1*x + k2*y + k3*z)
                 
                 v1 = -v1; v2 = -v2; v3 = -v3;
                  
             elseif(tag_func(fn) .eq. 203) then 
                v1 = - k3*lambda*n1*sin(k1*x + k2*y + k3*z) - k1*mu*n3*sin(k1*x + k2*y + k3*z)
                v2 = - k3*lambda*n2*sin(k1*x + k2*y + k3*z) - k2*mu*n3*sin(k1*x + k2*y + k3*z)
                v3 = - n3*(k3*lambda*sin(k1*x + k2*y + k3*z) + 2*k3*mu*sin(k1*x + k2*y + k3*z)) &
                     - k1*mu*n1*sin(k1*x + k2*y + k3*z) - k2*mu*n2*sin(k1*x + k2*y + k3*z)
                
                v1 = -v1; v2 = -v2; v3 = -v3;
              else
                write(*,*) 'Case not implemented!!! Check inputs.'
             
             endif
              ! v1 = 0; v2 = 0; v3 = 0;
             if (n_test.gt.0) then
               
              v1 = n1*(lambda*(2*pi*cos(pi*y)*sin(2*pi*x)*sin(pi*y)*sin(2*pi*z) &
                   - 2*pi*cos(pi*x)*sin(pi*x)*sin(2*pi*y)*sin(2*pi*z) &
                   + 2*pi*cos(pi*z)*sin(2*pi*x)*sin(2*pi*y)*sin(pi*z)) - &
                   4*mu*pi*cos(pi*x)*sin(pi*x)*sin(2*pi*y)*sin(2*pi*z)) &
                   + mu*n2*(2*pi*cos(2*pi*x)*sin(pi*y)**2*sin(2*pi*z)  &
                   - 2*pi*cos(2*pi*y)*sin(pi*x)**2*sin(2*pi*z)) &
                   + mu*n3*(2*pi*cos(2*pi*x)*sin(2*pi*y)*sin(pi*z)**2 &
                   - 2*pi*cos(2*pi*z)*sin(pi*x)**2*sin(2*pi*y))
             v2 = n2*(lambda*(2*pi*cos(pi*y)*sin(2*pi*x)*sin(pi*y)*sin(2*pi*z) & 
                - 2*pi*cos(pi*x)*sin(pi*x)*sin(2*pi*y)*sin(2*pi*z) &
                + 2*pi*cos(pi*z)*sin(2*pi*x)*sin(2*pi*y)*sin(pi*z)) &
                + 4*mu*pi*cos(pi*y)*sin(2*pi*x)*sin(pi*y)*sin(2*pi*z)) &
                + mu*n1*(2*pi*cos(2*pi*x)*sin(pi*y)**2*sin(2*pi*z) &
                - 2*pi*cos(2*pi*y)*sin(pi*x)**2*sin(2*pi*z)) &
                + mu*n3*(2*pi*cos(2*pi*y)*sin(2*pi*x)*sin(pi*z)**2 &
                + 2*pi*cos(2*pi*z)*sin(2*pi*x)*sin(pi*y)**2)
             v3 = n3*(lambda*(2*pi*cos(pi*y)*sin(2*pi*x)*sin(pi*y)*sin(2*pi*z) &
                - 2*pi*cos(pi*x)*sin(pi*x)*sin(2*pi*y)*sin(2*pi*z) &
                + 2*pi*cos(pi*z)*sin(2*pi*x)*sin(2*pi*y)*sin(pi*z)) &
                + 4*mu*pi*cos(pi*z)*sin(2*pi*x)*sin(2*pi*y)*sin(pi*z)) &
                + mu*n1*(2*pi*cos(2*pi*x)*sin(2*pi*y)*sin(pi*z)**2 &
                - 2*pi*cos(2*pi*z)*sin(pi*x)**2*sin(2*pi*y)) &
                + mu*n2*(2*pi*cos(2*pi*y)*sin(2*pi*x)*sin(pi*z)**2 &
                + 2*pi*cos(2*pi*z)*sin(2*pi*x)*sin(pi*y)**2)
               endif
                
                fmat(fn,(3*(in -1) +1)) = fmat(fn,(3*(in -1) +1)) +  0.25*area * v1 * ww(i)*ww(j)  
                fmat(fn,(3*(in -1) +2)) = fmat(fn,(3*(in -1) +2)) +  0.25*area * v2 * ww(i)*ww(j)  
                fmat(fn,(3*(in -1) +3)) = fmat(fn,(3*(in -1) +3)) +  0.25*area * v3 * ww(i)*ww(j)   
                                       
          enddo  
       enddo   
    endif  
 endif 
enddo  

endif   

! ------------------------------------------------------------------------

if (nl_dirX.gt.0) then    ! Dirichlet X
do il = 1,nl_dirX
 if (tag_dirX(il).eq.cs_bc_loc(cs_bc_loc(ie -1) +0)) then
    fn = 0
    do ifun = 1,nfunc
       if (fun_dirX(il).eq.tag_func(ifun)) fn = ifun
    enddo
    
    if (fn.gt.0) then
       v1 = val_dirX(il,1)
       v2 = val_dirX(il,2)
       v3 = val_dirX(il,3)
       v4 = val_dirX(il,4)
       
       do j = 1,nn
          do i = 1,nn
             is = nn*(j -1) +i
             in = cs_bc_loc(cs_bc_loc(ie -1) +is)
                v = 0.25*(1.0 -ct(i))*(1.0 -ct(j))*v1 &
                  + 0.25*(1.0 +ct(i))*(1.0 -ct(j))*v2 &
                  + 0.25*(1.0 +ct(i))*(1.0 +ct(j))*v3 &
                  + 0.25*(1.0 -ct(i))*(1.0 +ct(j))*v4
                
                fmat(fn,(3*(in -1) +1)) = v
          enddo
       enddo
    endif
 endif
enddo
endif

if (nl_dirY.gt.0) then    ! Dirichlet Y
do il = 1,nl_dirY
 if (tag_dirY(il).eq.cs_bc_loc(cs_bc_loc(ie -1) +0)) then
    fn = 0
    do ifun = 1,nfunc
       if (fun_dirY(il).eq.tag_func(ifun)) fn = ifun
    enddo
    
    if (fn.gt.0) then
       v1 = val_dirY(il,1)
       v2 = val_dirY(il,2)
       v3 = val_dirY(il,3)
       v4 = val_dirY(il,4)
       
       do j = 1,nn
          do i = 1,nn
             is = nn*(j -1) +i
             in = cs_bc_loc(cs_bc_loc(ie -1) +is)
                v = 0.25*(1.0 -ct(i))*(1.0 -ct(j))*v1 &
                  + 0.25*(1.0 +ct(i))*(1.0 -ct(j))*v2 &
                  + 0.25*(1.0 +ct(i))*(1.0 +ct(j))*v3 &
                  + 0.25*(1.0 -ct(i))*(1.0 +ct(j))*v4
                
                fmat(fn,(3*(in -1) +2)) = v
          enddo
       enddo
    endif
 endif
enddo
endif


if (nl_dirZ.gt.0) then    ! Dirichlet Z
do il = 1,nl_dirZ
 if (tag_dirZ(il).eq.cs_bc_loc(cs_bc_loc(ie -1) +0)) then
    fn = 0
    do ifun = 1,nfunc
       if (fun_dirZ(il).eq.tag_func(ifun)) fn = ifun
    enddo
    
    if (fn.gt.0) then
       v1 = val_dirZ(il,1)
       v2 = val_dirZ(il,2)
       v3 = val_dirZ(il,3)
       v4 = val_dirZ(il,4)
       
       do j = 1,nn
          do i = 1,nn
             is = nn*(j -1) +i
             in = cs_bc_loc(cs_bc_loc(ie -1) +is)
                v = 0.25*(1.0 -ct(i))*(1.0 -ct(j))*v1 &
                  + 0.25*(1.0 +ct(i))*(1.0 -ct(j))*v2 &
                  + 0.25*(1.0 +ct(i))*(1.0 +ct(j))*v3 &
                  + 0.25*(1.0 -ct(i))*(1.0 +ct(j))*v4
                
                fmat(fn,(3*(in -1) +3)) = v
          enddo
       enddo
    endif
 endif
enddo
endif
deallocate(ct,dd,ww)
enddo
endif

if (nl_neuN.gt.0) deallocate(i4normal, normal_nx_el_neuN, normal_ny_el_neuN, normal_nz_el_neuN)
if (nl_poiZ.gt.0) deallocate(node_counter_poiZ)
if (nl_poiY.gt.0) deallocate(node_counter_poiY)
if (nl_poiX.gt.0) deallocate(node_counter_poiX)

!      deallocate(ct,ww,dd)

!      do i = 1, nfunc 
!         do j = 1, nnod_loc
!           write(*,*) fmat(i,3*(j-1)+1), fmat(i,3*(j-1)+2),fmat(i,3*(j-1)+3)
!         enddo
!         read(*,*)
!      enddo 
!stop


return
end subroutine MAKE_EXTINT_FORCES



