subroutine msd( dr_dt_3dmsd, dx_dt_3dmsd, dy_dt_3dmsd, dz_dt_3dmsd, n_dt_3dmsd, &
                   & ID, fname, &
                   & frame, aid, tot_frame, tot_atom, &
                   & maxdtsteps, n_ID)

        implicit none

        character(LEN=200), intent(in) :: fname

        integer, intent(in) :: frame, tot_frame, &
                   & tot_atom, n_ID, maxdtsteps, aid
        
        integer, intent(in), dimension(0:n_ID-1) :: ID
        
        real(4), intent(out), dimension(0:maxdtsteps-1) :: dr_dt_3dmsd, &
                   & dx_dt_3dmsd, dy_dt_3dmsd, dz_dt_3dmsd, n_dt_3dmsd
        
        real(4), dimension(0:tot_atom-1,0:2) :: r0, rt
        real(8), dimension(0:2) :: box0 = 0, binv0 = 0, hbox0 = 0
        real(8), dimension(0:2) :: boxt = 0, binvt = 0, hboxt = 0
        real(8), dimension(0:5) :: box5 = 0

        !real(4) :: dxt = 0.0, dyt = 0.0, dzt = 0.0, binctobin=0.0, &
        !           & dr2=0.0, zi0 = 0, zit = 0, fnskip=1.
        !integer :: nxt = 0, nyt = 0, nzt = 0, zref = 0, n = 0, &
        !           & m = 0, m2 = 0, i = 0, j = 0, k = 0, i0 = 0, &
        !           & it = 0, nskip = 1, &
        !           & nbroad = 0, nreplica = 0, &
        !           & lownrep = 0, highnrep = 0, &
        !           & ireplica = 0, ibroad = 0, lagt = 0
        integer :: i = 0, m = 0, n = 0, lagt = 0
        real(4) :: dr2
        real(4), dimension(0:2) :: dr0 = 0, drt = 0, dr0t = 0
        real :: T1, T2
        !z0 = 0
        !zt = 0
        !valz = 0
        !flag = 0
        !tflag = 0
        r0 = 0
        rt = 0
        
        !highnrep = int((nreplica-1)/2)+1
        !lownrep = -1*highnrep
        !nbroad = int(nbins/nreplica)
        print *, "Inside Alldyn"
        !binctobin = int(0.5*dbinwidth/dbincwidth)
        !zref = (nbins+1)/2-1
        !valz(:) =   
        dr_dt_3dmsd = 0
        dx_dt_3dmsd = 0
        dy_dt_3dmsd = 0
        dz_dt_3dmsd = 0
        n_dt_3dmsd = 0
        open(unit = 435, file = trim(fname), form = "UNFORMATTED", status = "OLD")
        read(435) 
        read(435)
        read(435)
        do i = 1, frame-1
                read(435)
                read(435)
                read(435)
                read(435)
        end do
        !!!!!!!!!!!!!!!!!!!!!!! Frames skipped !!!!!!!!!!!
        call cpu_time(T1)
        
         r0 = 0.0         
         read(435) box5
         read(435) r0(:,0)
         read(435) r0(:,1)
         read(435) r0(:,2)
         
         box0(0) = box5(0)
         box0(1) = box5(2)
         box0(2) = box5(5)
         binv0 = 1.0/box0
         hbox0 = box0*0.5
         !call pbcwrapxy(r0,r0(aid-1,:),tot_atom,box0,binv0)
        
         !flag = 1.0
         n = 0
        
        do m = frame+1, min(frame+1 + maxdtsteps-3,tot_frame) !-3*nskip, 1 !nskip
         !print *, m
         rt = 0.0         
         n = n+1             !!!! Keep count of frames read
         lagt = (m-frame-1)
         read(435) box5
         read(435) rt(:,0)
         read(435) rt(:,1)
         read(435) rt(:,2)
         !print *, "read frame"
         boxt(0) = box5(0)
         boxt(1) = box5(2)
         boxt(2) = box5(5)
         binvt = 1.0/boxt
         hboxt = boxt*0.5
         !call pbcwrapxy(rt,rt(aid-1,:),tot_atom,boxt,binvt)
         !print *, "pbc wrapped"
         do i = 0,n_ID-1
           !it = rt(ID(i)-1,2)
           !i0 = r0(ID(i)-1,2)
           !z0 = 0
           !zt = 0

           !do ireplica = lownrep,highnrep
           !     zit = nint((it - (ireplica-1)*dbincwidth)*1./dbinwidth)*nreplica + (ireplica-1) + zref !nint((it)*1./dbincwidth) + zref
           !     zt(zit) = 1.0
           !     zi0 = nint((i0 - (ireplica-1)*dbincwidth)*1./dbinwidth)*nreplica + (ireplica-1) + zref !nint((i0)*1./dbincwidth) + zref
           !     z0(zi0) = 1.0
           !end do
        
           !flag(:,i) = flag(:,i)*z0(:)*zt(:)
           !tflag(:,i) = tflag(:,i) + flag(:,i)
        
           !!!!!!!!! Two-dimensional MSD !!!!!!!!!!!!!!!
           dr0t(:) = rt(ID(i)-1,:) - r0(ID(i)-1,:)

           dr2 = dr0t(0)*dr0t(0)
           dx_dt_3dmsd(lagt) = dx_dt_3dmsd(lagt)+dr2 !-K_dt_2dmsd(:,lagt))*flag(:,i)
           !dx2_dt_3dmsd(lagt) = dx2_dt_2dmsd(:,lagt)+(dr2-K_dt_2dmsd(:,lagt))*(dr2-K_dt_2dmsd(:,lagt))*flag(:,i)

           dr2 = dr0t(1)*dr0t(1)
           dy_dt_3dmsd(lagt) = dy_dt_3dmsd(lagt)+(dr2) !-K_dt_2dmsd(:,lagt))*flag(:,i)
           !dy2_dt_3dmsd(lagt) = dy2_dt_3dmsd(lagt)+(dr2) !-K_dt_2dmsd(:,lagt))*(dr2-K_dt_2dmsd(:,lagt))*flag(:,i)

           dr2 = dr0t(2)*dr0t(2)
           dz_dt_3dmsd(lagt) = dz_dt_3dmsd(lagt)+(dr2) !-K_dt_2dmsd(:,lagt))*flag(:,i)
           !dz2_dt_3dmsd(lagt) = dz2_dt_3dmsd(lagt)+(dr2) !-K_dt_2dmsd(:,lagt))*(dr2-K_dt_2dmsd(:,lagt))*flag(:,i)

           dr2 = dot_product(dr0t,dr0t)
           dr_dt_3dmsd(lagt) = dr_dt_3dmsd(lagt)+(dr2) !-K_dt_2dmsd(:,lagt))*flag(:,i)
           !dr2_dt_3dmsd(lagt) = dr2_dt_3dmsd(lagt)+(dr2) !-K_dt_2dmsd(:,lagt))*(dr2-K_dt_2dmsd(:,lagt))*flag(:,i)

           n_dt_3dmsd(lagt) = n_dt_3dmsd(lagt)+1.0 !*flag(:,i)
        
           !!!!!!!!!!  Rotational correlation !!!!!!!!!!
           
           !!!!!!!!!!  Dipole vector contribution 1
           !dr0(:) = r0(ID(i)-1,:)-r0(ID(i),:)
           !dr0(0) = dr0(0) - box0(0)*anint(dr0(0)*binv0(0))
           !dr0(1) = dr0(1) - box0(1)*anint(dr0(1)*binv0(1))
           !dr0(2) = dr0(2) - box0(2)*anint(dr0(2)*binv0(2))
           !!!!!!!!!!  Dipole vector contribution 2
           !dr0t(:) = r0(ID(i)-1,:)-r0(ID(i)+1,:)
           !dr0t(0) = dr0t(0) - box0(0)*anint(dr0t(0)*binv0(0))
           !dr0t(1) = dr0t(1) - box0(1)*anint(dr0t(1)*binv0(1))
           !dr0t(2) = dr0t(2) - box0(2)*anint(dr0t(2)*binv0(2))
           !!!!!!!!!! Dipole vector at 0
           !dr0(:) = 0.5*(dr0(:)+dr0t(:)) 
           !dr0(0) = dr0(0) - box0(0)*anint(dr0(0)*binv0(0))
           !dr0(1) = dr0(1) - box0(1)*anint(dr0(1)*binv0(1))
           !dr0(2) = dr0(2) - box0(2)*anint(dr0(2)*binv0(2))
           !!!!!!!!!! Dipole vector contribution 1
           !drt(:) = rt(ID(i)-1,:)-rt(ID(i),:)
           !drt(0) = drt(0) - boxt(0)*anint(drt(0)*binvt(0))
           !drt(1) = drt(1) - boxt(1)*anint(drt(1)*binvt(1))
           !drt(2) = drt(2) - boxt(2)*anint(drt(2)*binvt(2))
           !!!!!!!!!! Dipole vector contribution 2
           !dr0t(:) = rt(ID(i)-1,:)-rt(ID(i)+1,:)
           !dr0t(0) = dr0t(0) - boxt(0)*anint(dr0t(0)*binvt(0))
           !dr0t(1) = dr0t(1) - boxt(1)*anint(dr0t(1)*binvt(1))
           !dr0t(2) = dr0t(2) - boxt(2)*anint(dr0t(2)*binvt(2))
           !!!!!!!!!! Dipole vector at t
           !drt(:) = 0.5*(drt(:)+dr0t(:)) 
           !drt(0) = drt(0) - boxt(0)*anint(drt(0)*binvt(0))
           !drt(1) = drt(1) - boxt(1)*anint(drt(1)*binvt(1))
           !drt(2) = drt(2) - boxt(2)*anint(drt(2)*binvt(2))
           !!!!!!!!!! Correlation calculation
           !dr2 = dot_product(dr0,drt)/NORM2(dr0)/NORM2(drt)

           !call pbcdr(dr0,r0(ID(i)-1,:),r0(ID(i),:),box0,binv0)
           !call pbcdr(dr0t,r0(ID(i)-1,:),r0(ID(i)+1,:),box0,binv0)
           !!!!!!!!!! Dipole vector at 0
           !dr0 = 0.5*(dr0+dr0t)
           !dr0 = dr0 - box0*anint(dr0*binv0)
        
           !!!!!!!!!! Dipole vector contribution 1
           !call pbcdr(drt,rt(ID(i)-1,:),rt(ID(i),:),boxt,binvt)
           !!!!!!!!!! Dipole vector contribution 2
           !call pbcdr(dr0t,rt(ID(i)-1,:),rt(ID(i)+1,:),boxt,binvt)
           !!!!!!!!!! Dipole vector at t
           !drt = 0.5*(drt+dr0t)
           !drt = drt - boxt*anint(drt*binvt) 
   
           !dr2 = dot_product(dr0,drt)/NORM2(dr0)/NORM2(drt)
           !dr_dt_rc(:,lagt) = dr_dt_rc(:,lagt) + (dr2-K_dt_rc(:,lagt))*flag(:,i)
           !dr2_dt_rc(:,lagt) = dr2_dt_rc(:,lagt) + (dr2-K_dt_rc(:,lagt))*(dr2-K_dt_rc(:,lagt))*flag(:,i)
           !n_dt_rc(:,lagt) = n_dt_rc(:,lagt)+1.0*flag(:,i)
        
         end do
        
        call cpu_time(T2)


        !do m2 = 1, nskip-1
        ! !print *, m
        ! !n = n+1             !!!! Keep count of frames read
        ! read(435) box5
        ! read(435) rt(:,0)
        ! read(435) rt(:,1)
        ! read(435) rt(:,2)

        !end do

        end do
        close(435)
        
        !!!!!!!!! Residence Time !!!!!!!!!!!!!!!!!!!!
        !do i=0,nbins-1
        !do j=0,n_ID-1
        !n_dt_RT(i,tflag(i,j)-1) = n_dt_RT(i,tflag(i,j)-1) + 1.0
        !end do
        !end do
end subroutine msd