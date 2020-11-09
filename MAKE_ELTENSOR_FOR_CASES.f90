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

!> @brief Assignes material properties node by node. 
!! @author Ilario Mazzieri
!> @date September, 2014 
!> @version 1.0
!> @param[in] tcase label for not honoring case  
!> @param[in] vcase value case for non linear block 
!> @param[in] nn  polynomial degree
!> @param[in] nn_loc number of local nodes
!> @param[in] zs_elev elevation of the nodes from topography
!> @param[in] zs_all elevation of the nodes from alluvial
!> @param[in] vs_nodes vs30 of the nodes 
!> @param[in] thick_nodes thickness of sediments
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc connectivity vector
!> @param[in] ielem  element index
!> @param[in] nf number of functions
!> @param[in] sub_tag_all tags for multi-not honoring
!> @param[in] xs -coordinate of the nodes 
!> @param[in] ys y-coordinate of the nodes 
!> @param[in] zs z-coordinate of the nodes 
!> @param[in] mpi_id MPI process id 
!> @param[in] local_n_num local node numeration
!> @param[in] damping_tpe 1-Kosloff&Kosloff, 2-Standard Linear Solid
!> @param[in] check_case 1-debug mode, 0-standard mode
!> @param[out] rho_el material density 
!> @param[out] lambda_el Lame coefficient lambda
!> @param[out] mu_el Lame coefficient mu
!> @param[out] gamma_el damping coefficient gamma (Kosloff&Kosloff)
!> @param[out] qs quality factor for s-waves (Standard Linear Solid)
!> @param[out] qp quality factor for p-waves (Standard Linear Solid)

      subroutine MAKE_ELTENSOR_FOR_CASES(tcase,vcase,&
                                 nn,rho_el,lambda_el,mu_el,gamma_el,&
                                 nn_loc, zs_elev, zs_all, vs_nodes, thick_nodes, &
                                 cs_nnz_loc, cs_loc, ielem, &
                                 sub_tag_all, zs, mpi_id, local_n_num, &
                                 damping_type, qs, qp, &
                                 xs, ys, check_case, label_case)
 
      
      implicit none
                                                      
      integer*4 :: tcase, check_case, label_case               
      integer*4 :: vcase, mpi_id        
      integer*4 :: nn
      integer*4 :: p, q, r, ic
      integer*4 :: nn_loc
      integer*4 :: cs_nnz_loc                                
      integer*4 :: is,in,ielem , damping_type                               

      integer*4, dimension(nn_loc) :: local_n_num
      integer*4, dimension(nn_loc) :: sub_tag_all
      integer*4, dimension(0:cs_nnz_loc) :: cs_loc                

      real*8 :: Depth, Depth_real, vs_all, vp_all, thickness, vs30, pig
      real*8 :: VS,VP, rho,lambda,mu,gamma,ni, qs,qp
      real*8 :: x1,y1,x2,y2,coef_a, coef_b, coef_c, numer, den, distance, f_distance

      real*8, dimension(nn_loc) :: zs_elev
      real*8, dimension(nn_loc) :: zs_all
      real*8, dimension(nn_loc) :: vs_nodes, thick_nodes

      real*8, dimension(nn_loc) :: zs, xs, ys                

      real*8, dimension(nn,nn,nn) :: rho_el,lambda_el,mu_el,gamma_el

      real*8 :: stat_id1_x, stat_id1_y, stat_id2_x, stat_id2_y
      
      character*70 :: filename
      character*5 :: filesuffix
         
      
      pig = 4.d0*datan(1.d0);
      
!     STRESS CALCULATION
      
      if (check_case .eq. 1) then
      
          filesuffix = '.dat'
          write(filename, '(A,I5.5,A5)') 'NHCheck', mpi_id, filesuffix
          !write(*,*) mpi_id, ielem, filename
          !read(*,*)
          open(1000 + mpi_id,file=filename,position='APPEND')         
      endif
      
      do r = 1,nn
          do q = 1,nn
              do p = 1,nn
                  is = nn*nn*(r -1) +nn*(q -1) +p
                  in = cs_loc(cs_loc(ielem -1) +is)
                  ic = in

                  if (ic .eq. 0 ) write(*,*) 'Error in MAKE_ELTENSOR_FOR_CASES '

                  if (tcase.eq.1) then
                  ! CASE 1: GRENOBLE 1
                      call MAKE_MECH_PROP_CASE_001(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                   xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                   vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))
                                               
                  elseif (tcase.eq.2) then
                  ! CASE 2: GRENOBLE 2
                      call MAKE_MECH_PROP_CASE_002(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                   xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                   vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))

                  elseif (tcase.eq.3) then
                  ! CASE 3: GUBBIO
                      call MAKE_MECH_PROP_CASE_003(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                   xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                   vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))
                                
                  elseif (tcase.eq.4) then
                  ! CASE 4: SULMONA NOT honoring
                      call MAKE_MECH_PROP_CASE_004(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                   xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                   vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))
                                
                  elseif (tcase.eq.5) then
                  ! CASE 5: VOLVI CASHIMA benchmark -  NOT honoring
                      call MAKE_MECH_PROP_CASE_005(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                   xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                   vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))
                
                                        
                  elseif (tcase.eq.6) then !NOT FOUND
                  ! CASE 6: FRIULI (Tagliamento river valley) 
                      call MAKE_MECH_PROP_CASE_006(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                   xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                   vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))

                  elseif (tcase.eq.7) then
             	  ! CASE 7: AQUILA (Evangelista et al. 2017)
                       call MAKE_MECH_PROP_CASE_007(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                    xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                    vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))

    		       elseif (tcase.eq.8) then
                   ! CASE 8: SANTIAGO
                       call MAKE_MECH_PROP_CASE_008(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                    xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                    vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))
			
                   elseif (tcase.eq.10) then
                   ! CASE 10: CHRISTCHURCH INGV - Staircase
                        call MAKE_MECH_PROP_CASE_010(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                     xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                     vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))
                                
                   elseif (tcase.eq.11) then
                   ! CASE 11: CHRISTCHURCH (last model)
                        call MAKE_MECH_PROP_CASE_011(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                      xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))

                   elseif (tcase.eq.12) then
                   ! CASE 12: PO PLAIN (new model)  

                        !straight line passing by (x1,y1,0), (x2,y2,0)
                        x1 =  654957.002352;  y1 = 4974060.299450;
                        x2 =  688420.525202;  y2 = 4957613.600935;
                        !coefficient of the line ax + by + c = 0
                        coef_a = 1.d0/(x2-x1);
                        coef_b = 1.d0/(y1-y2);
                        coef_c = - y1/(y1-y2) + x1/(x1-x2);
                        !distance between (x,y,0) and the line ax + by + c = 0
                        numer = coef_a*xs(ic) + coef_b*ys(ic) + coef_c
                        den = dsqrt(coef_a**2 + coef_b**2)
                        distance = dabs(numer/den)
                        f_distance = 150.d0 + 1850.d0/(1.d0 + dexp(-0.0012*(distance-5000.d0)));    

                        call MAKE_MECH_PROP_CASE_012(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                      xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic), f_distance)


                   elseif (tcase.eq.13) then
                   ! CASE 13: PO PLAIN-BEDROCK  
                        call MAKE_MECH_PROP_CASE_013(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                      xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))
                        
                   elseif (tcase.eq.14) then
                   ! CASE 14: WELLINGTON (simplified model, Benites et al. 2005)
                        call MAKE_MECH_PROP_CASE_014(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                      xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))
                        
                   elseif (tcase.eq.15) then
                    ! CASE 15: MARSICA  
                        call MAKE_MECH_PROP_CASE_015(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                      xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))
                                                                
                   elseif (tcase.eq.16) then
                    ! CASE 16: ISTANBUL 
                         call MAKE_MECH_PROP_CASE_016(rho, lambda, mu, gamma, qs, qp, & !outputs
                                                      xs(ic), ys(ic), zs(ic),&
                                                      zs_elev(ic), zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))
                                                 
                   elseif (tcase.eq.18) then
                   ! CASE 18: BEIJING -- SIMPLIFIED MODEL
                        call MAKE_MECH_PROP_CASE_018(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                      xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))

                   elseif (tcase.eq.19) then
                   ! CASE 19: SALONICCO
                        call MAKE_MECH_PROP_CASE_019(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                      xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))

                   elseif (tcase.eq.20) then
                   ! CASE 20: ATENE  
                        call MAKE_MECH_PROP_CASE_020(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                      xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))

                   elseif (tcase.eq.21) then
                   ! CASE 21: BEIJING
                        call MAKE_MECH_PROP_CASE_021(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                      xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))

                   elseif (tcase.eq.22) then
                   !-------------------------------------------------------------------
                   ! CASE 22: NORCIA
                        call MAKE_MECH_PROP_CASE_022(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                      xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))

                   elseif (tcase.eq.30) then
                   !-------------------------------------------------------------------
                   ! CASE 30: ATHENS - PARTHENON
                        call MAKE_MECH_PROP_CASE_030(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                      xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))

                   elseif (tcase.eq.31) then
                   !------------------------------------------------------------------
                   ! CASE 31: GRONINGEN
                        call MAKE_MECH_PROP_CASE_031(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                      xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))
					 
                   elseif (tcase.eq.32) then
                   !------------------------------------------------------------------
                   ! CASE 32: GRONINGEN - LAYERED MODEL
                        call MAKE_MECH_PROP_CASE_032(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                      xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))

                   elseif (tcase.eq.33) then
                   !------------------------------------------------------------------
                   ! CASE 33: GRONINGEN - 2ND LAYERED MODEL
                        call MAKE_MECH_PROP_CASE_033(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                      xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))

                   elseif (tcase.eq.40) then
                   !-------------------------------------------------------------------
                   ! CASE 40: KUTCH  
                        call MAKE_MECH_PROP_CASE_040(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                      xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))

                   elseif (tcase.eq.50) then
                   !-------------------------------------------------------------------
                   ! CASE 50: PLANE-WAVE benchmark -  MULTI NOT honoring
                        call MAKE_MECH_PROP_CASE_040(rho,lambda,mu,gamma,qs,qp, & !outputs
                                                      xs(ic),ys(ic),zs(ic),zs_elev(ic),zs_all(ic), &
                                                      vs_nodes(ic), thick_nodes(ic), sub_tag_all(ic))
                    
                   elseif (tcase.eq.98) then
                   !-------------------------------------------------------------------
                   ! CASE 98: TEST honoring (only TOPOGRAPHY)
                        VS  = 100        !VS: S velocity in m/s
                        VP  = 200        !VS: S velocity in m/s
                        rho = 2000        !RHO: MASS DENSITY in kg/m^3
                        lambda = rho * (VP**2 - 2*VS**2)
                        mu = rho * VS**2
                        gamma = 0.0d0
                                
                   elseif (tcase.eq.99) then
                   !-------------------------------------------------------------------
                   ! CASE 99: TEST honoring (TOPOGRAPHY&ALLUVIAL)
                        Depth = zs_elev(ic)    !D: depth in m
                        if ((Depth .ge. 0.0d0) .and. (zs_all(ic) .ge. 0.0d0)) then
                            VS  = 1000.d0         !VS: S velocity in m/s
	                        VP  = 2081.9942d0     !VP: P velocity in m/s 
	                        rho = 2000.d0         !rho_el: MASS DENSITY in kg/m^3
	                        lambda = rho * (VP**2 - 2*VS**2)
	                        mu = rho * VS**2
	                        qp = 200;
	                        qs = 100;
                            gamma = 4.d0*datan(1.d0)/qs;
                        else    
                            VS  = 2000.d0         !VS: S velocity in m/s
	                        VP  = 3463.9976     !VP: P velocity in m/s 
	                        rho = 2500.d0         !rho_el: MASS DENSITY in kg/m^3
	                        lambda = rho * (VP**2 - 2*VS**2)
	                        mu = rho * VS**2
	                        qp = 200;
	                        qs = 100;
                            gamma = 4.d0*datan(1.d0)/qs;
                        endif
                                
                   elseif (tcase.eq.100) then
                   !-------------------------------------------------------------------
                   ! CASE 99: TEST plane wave (TOPOGRAPHY&ALLUVIAL)
                       Depth = zs_elev(ic)    !D: depth in m
                       if ((Depth .ge. 0.0d0) .and. (zs_all(ic) .ge. 0.0d0)) then
                           VS  = 300.d0         !VS: S velocity in m/s
	                       VP  = 600.d0     !VP: P velocity in m/s 
	                       rho = 1800.d0         !rho_el: MASS DENSITY in kg/m^3
	                       lambda = rho * (VP**2 - 2*VS**2)
	                       mu = rho * VS**2
	                       qp = 60;
	                       qs = 30;
                           gamma = 4.d0*datan(1.d0)/qs;
                       else    
                           VS  = 2000.d0         !VS: S velocity in m/s
	                       VP  = 4000.d0     !VP: P velocity in m/s 
	                       rho = 2200.d0         !rho_el: MASS DENSITY in kg/m^3
	                       lambda = rho * (VP**2 - 2*VS**2)
	                       mu = rho * VS**2
                           qp = 400;
	                       qs = 200;
                           gamma = 4.d0*datan(1.d0)/qs;
                       endif
                   endif 
                        
                   if (check_case .eq. 1)  &
                       write(1000+mpi_id,*) xs(ic),ys(ic),zs(ic), &
                                            dsqrt(mu/rho), &
                                            dsqrt((lambda + 2.d0*mu)/rho), &
                                            rho, lambda, mu, &
                                            qp, qs, gamma, zs_elev(ic), zs_all(ic)                              

                        
                   rho_el(p,q,r) = rho
                   lambda_el(p,q,r) = lambda
                   mu_el(p,q,r) = mu
                   gamma_el(p,q,r) = gamma
               enddo
           enddo
       enddo
     
       if (check_case .eq. 1) close (1000+mpi_id)

       if (damping_type .eq. 2) then
           qs = 0; qp = 0;
           vs_all = 0.d0; vp_all=0.d0;
           do r = 1,nn
               do q = 1,nn
                   do p = 1,nn
                       vs_all = vs_all + dsqrt(mu_el(p,q,r)/rho_el(p,q,r))
                       vp_all = vp_all + dsqrt( (lambda_el(p,q,r) + 2.d0*mu_el(p,q,r))/rho_el(p,q,r) );
                   enddo
               enddo
           enddo
       
           qs = 0.1d0*vs_all/nn**3;     
           qp = 0.1d0*vp_all/nn**3;

        endif

        return
      
       end subroutine MAKE_ELTENSOR_FOR_CASES

