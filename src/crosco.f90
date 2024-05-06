!!
!! crosco(par%nx,par%nz,npml,nx_pml,nz_pml,pad_top,xxb,zzb,xzb,xx,zz,xz,grdla,grdmu_a,grdmu_b,grdmu_c)

   subroutine crosco(nx,nz,npml,nx_pml,nz_pml,pad_top,txx,tzz,txz,rectxx,rectzz,rectxz,grdla,grdmu_a,grdmu_b,grdmu_c)
      implicit none
      integer, intent(in) :: nx,nz,npml,nx_pml,nz_pml,pad_top
      ! CROSCO - cross correlation

      integer iz,ix
      real(kind=kind(1.e0)) bacfld,forfld
      real(kind=kind(1.e0)) terma,termb,termc

      real,dimension(1:nz_pml,1:nx_pml)::grdla,grdmu_a,grdmu_b,grdmu_c
      real,dimension(1:nz_pml,1:nx_pml)::txx,tzz,txz,rectxx,rectzz,rectxz

!      print *,'inside crosco'
!      print *,'zsta',zsta, 'zend',zend,'xsta',xsta,'xend',xend

      do ix = npml+1,nx_pml-npml 
         do iz = npml+1,nz_pml-npml

            bacfld =    txx(iz,ix) +    tzz(iz,ix) ! Backward wavefield
            forfld = rectxx(iz,ix) + rectzz(iz,ix) ! Forward wavefield

            grdla(iz,ix)=grdla(iz,ix)+(forfld*bacfld)  ! Cross-correlation


            terma = (rectxx(iz,ix)*tzz(iz,ix))+(rectzz(iz,ix)*txx(iz,ix))
            grdmu_a(iz,ix)=grdmu_a(iz,ix)+terma

            termb = (rectxx(iz,ix)*txx(iz,ix))+(rectzz(iz,ix)*tzz(iz,ix))
            grdmu_b(iz,ix)=grdmu_b(iz,ix)+termb

            termc = rectxz(iz,ix)*txz(iz,ix)
            grdmu_c(iz,ix)=grdmu_c(iz,ix)+termc

         enddo
      enddo

    end subroutine crosco

    subroutine grconv_crosco(nx,nz,npml,nx_pml,nz_pml,pad_top,lambda,mu,grdla,grdmu,grdmu_a,grdmu_b,grdmu_c)

      implicit none
      integer, intent(in) :: nx,nz,npml,nx_pml,nz_pml,pad_top

      integer iz,ix
      real::denom,terma,termb,termc

      real,dimension(1:nz_pml,1:nx_pml)::grdla,grdmu,grdmu_a,grdmu_b,grdmu_c,&
                                         lambda,mu,lam2mu

      lam2mu = 0.0 

      lam2mu = lambda + 2.0 * mu

       do ix = npml+1,nx_pml-npml
          do iz = npml+1,nz_pml-npml

             grdla(iz,ix) = - grdla(iz,ix) / (2.0*(lambda(iz,ix)+mu(iz,ix)))
             grdla(iz,ix) = grdla(iz,ix) / (2.0*(lambda(iz,ix)+mu(iz,ix)))

             if (mu(iz,ix) > 0) then
                denom=(lambda(iz,ix)+mu(iz,ix))**2.0
                denom=8.0*(mu(iz,ix)**2.0)*denom

                terma=2.0*(lambda(iz,ix)*lam2mu(iz,ix))/denom
                terma=terma*grdmu_a(iz,ix)

                termb=((lambda(iz,ix)**2.0)+(lam2mu(iz,ix)**2))/denom
                termb=termb*grdmu_b(iz,ix)

                termc=1.0/(mu(iz,ix)**2.0)
                termc=termc*grdmu_c(iz,ix)
                grdmu(iz,ix) = terma - termb - termc
             else
                grdmu(iz,ix) = 0.0
             endif
          enddo
       enddo            
   
    end subroutine grconv_crosco
    

    subroutine grconv(nx,nz,npml,nx_pml,nz_pml,pad_top,dt,rho,winvp,winvs,grdvp,grdvs,grdla,grdmu)

     implicit none
     integer, intent(in) :: nx,nz,npml,nx_pml,nz_pml,pad_top       

     integer iz,ix
     real::terma,termb,dt

     real,dimension(1:nz_pml,1:nx_pml)::grdvp,grdvs,grdla,grdmu,&
                                        rho,rhoo,winvp,winvs

      rhoo = rho

      do ix = npml+1,nx_pml-npml
         do iz = npml+1,nz_pml-npml

            rhoo(iz,ix) = rhoo(iz,ix)
            grdvp(iz,ix) = 2.0*rhoo(iz,ix) * winvp(iz,ix) * &
                           grdla(iz,ix)

            ! Change polarity of gradient to allow all notation to agree with Heiner's

            grdvp(iz,ix) = - grdvp(iz,ix)

            terma=-4.0*rhoo(iz,ix)*winvs(iz,ix)*grdla(iz,ix)

            termb=2.0*rhoo(iz,ix)*winvs(iz,ix)*grdmu(iz,ix)

            grdvs(iz,ix) = terma+termb
            grdvs(iz,ix) = -grdvs(iz,ix)

!            rhoo(iz,ix) = dt/rhoo(iz,ix)

         enddo
      enddo            


    end subroutine grconv

