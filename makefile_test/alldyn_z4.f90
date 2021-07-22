subroutine alldyn_z( dr_dt_rc, n_dt_rc, &
                   & dr2_dt_rc, K_dt_rc, &
                   & dr_dt_2dmsd, dx_dt_2dmsd, dy_dt_2dmsd,n_dt_2dmsd, &
                   & dr2_dt_2dmsd, dx2_dt_2dmsd, dy2_dt_2dmsd, K_dt_2dmsd, &
                   & n_dt_RT, &
                   & fname, outname, &
                   & ID, dbinwidth, dbincwidth, &
                   & frame, aid, tot_frame, tot_atom, &
                   & maxdtsteps, nbins, n_ID)

        implicit none
        character(LEN=200), intent(in) :: fname
        character(LEN=200), intent(in) :: outname
        integer, intent(in) :: frame, tot_frame, &
                   & tot_atom, n_ID, maxdtsteps, nbins, aid
        integer, intent(in), dimension(0:n_ID-1) :: ID
        real(4), intent(in) :: dbinwidth, dbincwidth
        real(4), intent(inout), dimension(0:nbins-1,0:maxdtsteps-1) :: dr_dt_rc, &
                   & n_dt_rc, dr_dt_2dmsd, dx_dt_2dmsd, dy_dt_2dmsd, n_dt_2dmsd, n_dt_RT, &
                   & dr2_dt_rc, dr2_dt_2dmsd, dx2_dt_2dmsd, dy2_dt_2dmsd
        real(4), intent(in), dimension(0:nbins-1,0:maxdtsteps-1) :: K_dt_rc, K_dt_2dmsd
        real(4), dimension(0:nbins-1) :: z0, zt
        real(4), dimension(0:nbins-1,0:n_ID-1) :: flag, tflag
        real(4), dimension(0:tot_atom-1,0:2) :: r0, rt
        real(8), dimension(0:2) :: box0, binv0, hbox0
        real(8), dimension(0:2) :: boxt, binvt, hboxt
        real(8), dimension(0:5) :: box5
        real(4) :: dxt, dyt, dzt, binctobin=0.0, dr2=0.0,zi0,zit, fnskip=1.
        integer :: nxt, nyt, nzt, zref, n, m, m2, i, j, k, i0, it, nskip=1
        real(4), dimension(0:2) :: dr0, drt, dr0t
        real :: T1, T2
        !print *, "Inside Alldyn"
        r0(:,:) = 0.0
        flag(:,:) = 0.0
        tflag(:,:) = 0.0
        z0(:) = 0.0
        zt(:) = 0.0
        rt(:,:) = 0.0
        box0(:) = 0.0
        binv0(:) = 0.0
        hbox0(:) = 0.0
        boxt(:) = 0.0
        hboxt(:) = 0.0
        binvt(:) = 0.0
        !binctobin = int(0.5*dbinwidth/dbincwidth)
        zref = (nbins+1)/2-1

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
        
         r0(:,:) = 0.0         
         read(435) box5
         read(435) r0(:,0)
         read(435) r0(:,1)
         read(435) r0(:,2)
         
         box0(0) = box5(0)
         box0(1) = box5(2)
         box0(2) = box5(5)
         binv0(:) = 1.0/box0(:)
         hbox0(:) = box0(:)*0.5
         call pbcwrapxy(r0,r0(aid-1,:),tot_atom,box0,binv0)
        
         flag(:,:) = 1.0
         n = 0
        
        do m = frame+1, frame+1+maxdtsteps*nskip-3*nskip, nskip
         !print *, m
         rt(:,:) = 0.0         
         n = n+1             !!!! Keep count of frames read
         read(435) box5
         read(435) rt(:,0)
         read(435) rt(:,1)
         read(435) rt(:,2)
         !print *, "read frame"
         boxt(0) = box5(0)
         boxt(1) = box5(2)
         boxt(2) = box5(5)
         binvt(:) = 1.0/boxt(:)
         hboxt(:) = boxt(:)*0.5
         call pbcwrapxy(rt,rt(aid-1,:),tot_atom,boxt,binvt)
         !print *, "pbc wrapped"
         do i = 0,n_ID-1
           it = rt(ID(i)-1,2)
           i0 = r0(ID(i)-1,2)
           !print *, "Item id",i,ID(i)
           zit = nint((it)*1./dbincwidth) + zref
           zi0 = nint((i0)*1./dbincwidth) + zref
           z0(:) = 0.0
           z0(zi0) = 1.0
           !z0(zi0-binctobin:zi0+binctobin) = 1.0
           zt(:) = 0.0
           zt(zit) = 1.0
           !zt(zit-binctobin:zit+binctobin) = 1.0
        
           flag(:,i) = flag(:,i)*z0(:)*zt(:)
           tflag(:,i) = tflag(:,i) + flag(:,i)
        
           !!!!!!!!! Two-dimensional MSD !!!!!!!!!!!!!!!
           dr0t(:) = rt(ID(i)-1,:) - r0(ID(i)-1,:)

           dr2 = dr0t(0)*dr0t(0)
           dx_dt_2dmsd(:,int((m-frame-1))) = dx_dt_2dmsd(:,int((m-frame-1)))+(dr2-K_dt_2dmsd(:,int((m-frame-1))))*flag(:,i)
           dx2_dt_2dmsd(:,int((m-frame-1))) = dx2_dt_2dmsd(:,int((m-frame-1)))+(dr2-K_dt_2dmsd(:,int((m-frame-1))))*(dr2-K_dt_2dmsd(:,int((m-frame-1))))*flag(:,i)

           dr2 = dr0t(1)*dr0t(1)
           dy_dt_2dmsd(:,int((m-frame-1))) = dy_dt_2dmsd(:,int((m-frame-1)))+(dr2-K_dt_2dmsd(:,int((m-frame-1))))*flag(:,i)
           dy2_dt_2dmsd(:,int((m-frame-1))) = dy2_dt_2dmsd(:,int((m-frame-1)))+(dr2-K_dt_2dmsd(:,int((m-frame-1))))*(dr2-K_dt_2dmsd(:,int((m-frame-1))))*flag(:,i)


           dr2 = dot_product(dr0t(0:1),dr0t(0:1))
           dr_dt_2dmsd(:,int((m-frame-1))) = dr_dt_2dmsd(:,int((m-frame-1)))+(dr2-K_dt_2dmsd(:,int((m-frame-1))))*flag(:,i)
           dr2_dt_2dmsd(:,int((m-frame-1))) = dr2_dt_2dmsd(:,int((m-frame-1)))+(dr2-K_dt_2dmsd(:,int((m-frame-1))))*(dr2-K_dt_2dmsd(:,int((m-frame-1))))*flag(:,i)


           n_dt_2dmsd(:,int((m-frame-1))) = n_dt_2dmsd(:,int((m-frame-1)))+1.0*flag(:,i)
        
           !!!!!!!!!!  Rotational correlation !!!!!!!!!!
           !!!!!!!!!!  Dipole vector contribution 1
           dr0(:) = r0(ID(i)-1,:)-r0(ID(i),:)
           dr0(0) = dr0(0) - box0(0)*anint(dr0(0)*binv0(0))
           dr0(1) = dr0(1) - box0(1)*anint(dr0(1)*binv0(1))
           dr0(2) = dr0(2) - box0(2)*anint(dr0(2)*binv0(2))
           !!!!!!!!!!  Dipole vector contribution 2
           dr0t(:) = r0(ID(i)-1,:)-r0(ID(i)+1,:)
           dr0t(0) = dr0t(0) - box0(0)*anint(dr0t(0)*binv0(0))
           dr0t(1) = dr0t(1) - box0(1)*anint(dr0t(1)*binv0(1))
           dr0t(2) = dr0t(2) - box0(2)*anint(dr0t(2)*binv0(2))
           !!!!!!!!!! Dipole vector at 0
           dr0(:) = 0.5*(dr0(:)+dr0t(:)) 
           dr0(0) = dr0(0) - box0(0)*anint(dr0(0)*binv0(0))
           dr0(1) = dr0(1) - box0(1)*anint(dr0(1)*binv0(1))
           dr0(2) = dr0(2) - box0(2)*anint(dr0(2)*binv0(2))
           
           !!!!!!!!!! Dipole vector contribution 1
           drt(:) = rt(ID(i)-1,:)-rt(ID(i),:)
           drt(0) = drt(0) - boxt(0)*anint(drt(0)*binvt(0))
           drt(1) = drt(1) - boxt(1)*anint(drt(1)*binvt(1))
           drt(2) = drt(2) - boxt(2)*anint(drt(2)*binvt(2))
           !!!!!!!!!! Dipole vector contribution 2
           dr0t(:) = rt(ID(i)-1,:)-rt(ID(i)+1,:)
           dr0t(0) = dr0t(0) - boxt(0)*anint(dr0t(0)*binvt(0))
           dr0t(1) = dr0t(1) - boxt(1)*anint(dr0t(1)*binvt(1))
           dr0t(2) = dr0t(2) - boxt(2)*anint(dr0t(2)*binvt(2))
           !!!!!!!!!! Dipole vector at t
           drt(:) = 0.5*(drt(:)+dr0t(:)) 
           drt(0) = drt(0) - boxt(0)*anint(drt(0)*binvt(0))
           drt(1) = drt(1) - boxt(1)*anint(drt(1)*binvt(1))
           drt(2) = drt(2) - boxt(2)*anint(drt(2)*binvt(2))
           !!!!!!!!!! Correlation calculation
           dr2 = dot_product(dr0,drt)/NORM2(dr0)/NORM2(drt)
           dr_dt_rc(:,int((m-frame-1))) = dr_dt_rc(:,int((m-frame-1))) + (dr2-K_dt_rc(:,int((m-frame-1))))*flag(:,i)

           dr2_dt_rc(:,int((m-frame-1))) = dr2_dt_rc(:,int((m-frame-1))) + (dr2-K_dt_rc(:,int((m-frame-1))))*(dr2-K_dt_rc(:,int((m-frame-1))))*flag(:,i)

           n_dt_rc(:,int((m-frame-1))) = n_dt_rc(:,int((m-frame-1)))+1.0*flag(:,i)
        
         end do
        
        call cpu_time(T2)
        open(199, file=trim(outname)//"rc_progress.txt", &
                        & STATUS="REPLACE",ACTION="WRITE")
        write(199,FMT=*) "Progress of frames: ", &
                        & int(float(m-frame)/float(maxdtsteps*nskip)*100), &
                        & "% done ", &
                        & int((T2-T1)/float(m-frame)*float(maxdtsteps*nskip-m)/60.0), &
                        & " min left", &
                        & int((T2-T1)/float(m-frame)*float(maxdtsteps*nskip-m)/60.0/60.0), &
                        & " hrs left"
        close(199)


        do m2 = 1, nskip-1
         !print *, m
         !n = n+1             !!!! Keep count of frames read
         read(435) box5
         read(435) rt(:,0)
         read(435) rt(:,1)
         read(435) rt(:,2)

        end do

        end do
        close(435)
        
        !!!!!!!!! Residence Time !!!!!!!!!!!!!!!!!!!!
        do i=0,nbins-1
        do j=0,n_ID-1
        n_dt_RT(i,tflag(i,j)-1) = n_dt_RT(i,tflag(i,j)-1) + 1.0
        end do
        end do
end subroutine alldyn_z
