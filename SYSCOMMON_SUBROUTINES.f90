!    Copyright (c) 2012 The SPEED FOUNDATION
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

!---------------------Central Difference Scheme---------------------------
!> @brief Computes the displacement of the oscillator through central difference scheme
!> @author Aline Herlin, Srihari
!> @date November, 2020
!> @version 1.0
!> @param[in] m           mass of the system
!> @param[in] c           damping of system
!> @param[in] dT          time step for system
!> @param[in] U1, U0      displacement at time instants n and n-1
!> @param[in] P           equivalent force [-Mu''-Fint]
!> @param[out] U          displacement at time instant n+1

subroutine CENTRAL_DIFFERENCE(NDOF, M, M_inv, C, dT, U, U1, U0, P, flag_Minv)
    
    implicit none

    integer, intent(in) :: NDOF
    integer, intent(inout) :: flag_Minv
    integer :: idof, j

    real*8, dimension(NDOF, NDOF), intent(in) :: M, C
    real*8, dimension(NDOF), intent(in) :: U1, U0, P
    real*8, intent(in) :: dT

    real*8, dimension(NDOF, NDOF), intent(inout) :: M_inv
    real*8, dimension(NDOF), intent(inout) :: U
    real*8, dimension(NDOF, NDOF) :: m1, m2, m3, m4
    real*8, dimension(NDOF) :: T
    real*8 :: x, y
    
    ! Check if for some cases we need to change damping matrix (C); then we can pre define these variables (m1, m2, m3, m4); may be we dont need M,M1_inv,C matrices later
    U=0; 

    m1=M/dT/dT + C/2./dT    !!! coeff for u(n+1) -> 1/M_inv
    m2=-2.*M/dT/dT          !!! coeff for u(n)
    m3=M/dT/dT-C/2./dT      !!! coeff for u(n-1)

    if (flag_Minv.eq.0) then
        if (NDOF.gt.1) then
            call matinv(m1, m4, NDOF)
        elseif (NDOF.eq.1) then
            m4 = 1/m1
        endif
        M_inv = m4
        flag_Minv = 1
    else
        m4 = M_inv
    endif

    do idof= 1,NDOF
        x=0;    y=0;
        do j=1, NDOF;
            x = x + m2(idof,j)*U1(j);  
            y = y + m3(idof,j)*U0(j);
        enddo
        T(idof) = P(idof) -x -y
    enddo

    !For SDOF
    !T=P-m2*U1-m3*U0
    !U=T/m1

    do idof=1,NDOF
        do j=1,NDOF
            U(idof) = U(idof) + M_inv(idof,j)*T(j)
        enddo
    enddo

    return
end subroutine CENTRAL_DIFFERENCE

!------------------Inverse of a Matrix------------------------------------
subroutine matinv(A, B, N) 

    implicit none

    integer(4), intent(in) :: N
    integer(4) :: IS(N),JS(N)
    integer(4) :: I, J, K

    real(8), intent (in) :: A(N,N)
    real(8) :: B(N,N)
    real(8) :: D, T

    B=A
    do K=1,N
        D=0.0D0

        do I=K,N
            do J=K,N
                if(abs(B(I,J))>D) then
                    D=abs(B(I,J)); IS(K)=I; JS(K)=J
                end if
            enddo
        enddo

        do J=1,N
            T=B(K,J); B(K,J)=B(IS(K),J); B(IS(K),J)=T
        enddo

        do I=1,N
            T=B(I,K);   B(I,K)=B(I,JS(K));  B(I,JS(K))=T
        enddo

        B(K,K)=1.d0/B(K,K);
        do J=1,N; if(J.NE.K) B(K,J)=B(K,J)*B(K,K);  enddo
        
        do I=1,N
            if(I.NE.K) then
                do J=1,N; if(J.NE.K) B(I,J)=B(I,J)-B(I,K)*B(K,J); enddo
            endif
        enddo

        do I=1,N; if(I.NE.K) B(I,K)=-B(I,K)*B(K,K); enddo
    enddo

    do K=N,1,-1
        do J=1,N
            T=B(K,J); B(K,J)=B(JS(K),J); B(JS(K),J)=T
        enddo
        do I=1,N
            T=B(I,K); B(I,K)=B(I,IS(K)); B(I,IS(K))=T
        enddo
    enddo

    return
end subroutine matinv


!------------------Const-Law: Non-linear Shear Spring----------------------------
subroutine ksteel02(props,s,e,de,Et,statev,spd, yield, IDeath) !M, ndof
!
    implicit none
    real*8 E0,sy0,eta,mu,gama,esoft,alpha,beta,a_k,Omega
    real*8 emax,emin,ert,srt,erc,src,Ehc,Eh1,dt,dc,eu
    real*8 de,s,e,s0,Et,e_unload,sign,sy,evs,eve,epeak,smax,max
    real*8 sres,eres,x,e_slip,s_slip,e_close,s_close,srel,ET1
    real*8 smin,spd,strain_end
    
    !real*8 mu
    integer kon, yield, IDeath  !ndof
    !real*8 M(ndof,ndof)!
    real*8 props(10), statev(11)

    E0  = props(1) 
    sy0 = props(2) 
    eta = props(3) 
    mu  = props(4) 
    gama= props(5) 
    esoft=props(6) 
    alpha= props(7) 
    beta = props(8)
    a_k= props(9) 
    Omega= props(10) 
    emax  = statev(1) !maximum strain
    emin  = statev(2) !minimum strain
    ert   = statev(3) !strain at load reversal toward tension
    srt   = statev(4)
    erc   = statev(5) !strain at load reversal toward compression
    src   = statev(6)
    kon   = nint(statev(7)) !
    Ehc   = statev(8) !effective cummulative hysteresis energy
    Eh1   = statev(9) !hysteresis energy in a half cycle
    dt    = statev(10) !damage index for tension
    dc    = statev(11) !damage index for compression
      
    eu    = mu * sy0/E0 !characteristic ultimate strain


    if(Omega<=0.) Omega=0.5  
    if(a_k<0.) a_k=0. 
    if(eta<=0.) eta=1.d-6;      
    if(esoft>=0.) esoft=-1.d-6  
      
    if (kon.eq.0) then
        emax =  sy0/E0
        emin = -beta*sy0/E0
        if (de.ge.0.0) then
            kon = 1
        else
            kon = 2
        end if
    else if ((kon.eq.1).and.(de.lt.0.0)) then !Load reversal
            kon = 2
            if (s.gt.0.0) then
                erc = e
                src = s
            end if
            Ehc = Ehc + Eh1 * (erc / eu ) ** 2.0
            Eh1 = 0.0 !a new half cycle is to begin
            if (e.gt.emax) emax = e
    else if ((kon.eq.2).and.(de.gt.0.0)) then !Load reversal
            kon = 1
            if (s.lt.0.0) then
                ert = e
                srt = s
            endif
            Ehc = Ehc + Eh1 * (ert / eu ) ** 2.0
            Eh1 = 0.0 !a new half cycle is to begin
            if (e.lt.emin) emin = e
    end if
!
    s0=s
    s = s + E0 * de
    Et = E0
    
    if(a_k>0.) then
        if(s0>0.) E_unload=E0*(abs(emax/(sy0/E0)))**(-a_k)
        if(s0<0.) E_unload=E0*(abs(emin/(sy0/E0)))**(-a_k) 
        if(E_unload<0.1*E0) E_unload=0.1*E0
    else
        E_unload=E0
    end if

    if(s0*de<0.) then 
        s = s0 + E_unload * de
        Et = E_unload
        if(s*s0<0.) then 
            de=de-s0/E_unload
            s0=1.D-6*sy0*sign(1.d0,s) 
            Et=E0
        end if
    end if

        
    if ( de .ge. 0.0 .and. s0>=0.) then
        sy = (1.0 - dt) * sy0
        !loading envelope
        ! Hardening
        if(e+de>sy/E0) then
            evs = max( sy + ( e + de - sy/E0) * eta * E0, 0.)
            evE = eta * E0
           if (s .ge. evs) then
              s = evs
              Et = evE
              yield=1
           end if
        end if
        ! Softening
        epeak=sy/E0+(alpha-1.)*sy/E0/eta
        if(e+0.5*de>epeak) then
            evs=max(sy*alpha+esoft*E0*(e+de-epeak),0.0*sy)
            if(sy*alpha+esoft*E0*(e+de-epeak)<=0.) IDeath=1 ! Complete damage
            evE=esoft*E0
            if (s .ge. evs) then
               s = evs
               Et = evE
               yield=1
            end if
        end if

        !reloading envelope
        smax = max(sy, sy + (emax - sy/E0) * eta * E0)  
        if(emax>epeak) then
            smax=max(sy*alpha+esoft*E0*(emax-epeak),0.0*sy)
        end if
        sres = 0.02 * smax                         
        eres = ert - (srt - sres) / E_unload                 

        x=emax-smax/E0  
        e_slip=gama*emax+(1.-gama)*x
        s_slip=smax*gama 
        e_close=e_slip*Omega 
        s_close=(e_close-eres)/(e_slip-eres) * (s_slip-sres) + sres

        if (eres .le. emax - smax / E0) then   
            if(e+0.5*de<e_close)  then  
                srel = (e+de-eres)/(e_slip-eres) * (s_slip-sres) + sres
                Et1=(s_slip-sres)/(e_slip-eres)
            else                 
                srel=(e+de-e_close)/(emax-e_close)*(smax-s_close)+s_close
                Et1 = (smax - s_close) / (emax - e_close)
            end if
            if (s .gt. srel) then
               s = max( srel, 0.)
               Et = Et1
            end if
        end if

    elseif ( de .lt. 0.0 .and. s0<0. ) then
        sy = (1.0 - dc) * sy0 *beta
        !loading envelope
        ! Hardening
        if(e+de<-sy/E0) then
            evs =  min(-sy + ( e + de + sy/E0) * eta * E0,-0.0*sy)
            evE = eta * E0
            if (s .le. evs) then
                s = evs
                Et = evE
                yield=1
            end if
        end if
        ! Softening
        epeak=-sy/E0-(alpha-1.)*sy/E0/eta
        if(e+0.5*de<epeak) then
            evs=min(-sy*alpha+esoft*E0*(e+de-epeak),-0.*sy)
            if(-sy*alpha+esoft*E0*(e+de-epeak)>=0.) IDeath=1 
            evE=esoft*E0
            if (s .le. evs) then
                s = evs
                Et = evE
                yield=1
          end if
        end if

        !reloading envelope
        smin = min(-sy, -sy + (emin + sy/E0) * eta * E0)
        if(emin<epeak)	then 
            smin=min(-sy*alpha+esoft*E0*(emin-epeak),0.)
        end if	
        sres = 0.02 * smin 
        eres = erc - (src - sres) /  E_unload

        x=emin-smin/E0 
        e_slip=gama*emin+(1.-gama)*x 
        s_slip=smin*gama 
        e_close=e_slip*Omega 
        s_close=(e_close-eres)/(e_slip-eres) * (s_slip-sres) + sres

        if (eres .ge. emin - smin / E0) then    
            if(e+0.5*de>e_close) then 
                srel = (e+de-eres)/(e_slip-eres) * (s_slip-sres) + sres
                Et1=(s_slip-sres)/(e_slip-eres)
            else         
                srel=(e+de-e_close)/(emin-e_close)*(smin-s_close)+s_close
                Et1 = (smin - s_close) / (emin - e_close)
            end if
            if (s .lt. srel) then
                s = min (srel, 0.)
                Et = Et1
            end if
        end if
    end if


    if (Et.ne.E0 .and. Et.ne. E_unload) then 
        spd = spd + s * de
        Eh1 = Eh1 + s * de
        if ( s .ge. 0.0 ) then
            dc = min(Ehc /(3.0 * beta* sy0 * eu), 0.7)
        else
            dt = min(Ehc /(3.0 * sy0 * eu), 0.7)
        end if
    end if

    x=max(sy0, beta*sy0)
    epeak=x/E0+(alpha-1.)*x/E0/eta
    Strain_End=epeak+abs(x/(esoft*E0))
    x=max(abs(emax),abs(emin),abs(e+de))

    if (IDeath==1) then
        s=0; 
        Et=1E-6*E0
    end if
    !
    statev(1)   = emax
    statev(2)   = emin
    statev(3)   = ert
    statev(4)   = srt
    statev(5)   = erc
    statev(6)   = src
    statev(7)   = kon
    statev(8)   = Ehc
    statev(9)   = Eh1
    statev(10)  = dt
    statev(11)  = dc
    return
end subroutine ksteel02