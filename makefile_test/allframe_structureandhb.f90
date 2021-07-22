subroutine convert_c_string_f_string(str,strlen,f_string)
        implicit none
        integer, intent(in):: strlen
        character(len=1), intent(in), dimension(1:strlen) :: str
        character(len=strlen), intent(out) :: f_string
        integer :: i
        f_string = ""
        do i = 1,strlen
        f_string = trim(f_string) // str(i)
        end do
        f_string = trim(f_string)
end subroutine

subroutine structure_z(aid, t_t_outname, O_ID, HO_ID1, HO_ID2, OH_ID, &
        & HOH_ID1, HOH_ID2, HbOH_ID1, HbOH_ID2, &
        & Os_ID, subs_ID, bO_ID, bOH_ID, filenames, dt, frame, &
        & n_O_ID, n_OH_ID, n_Os_ID, n_subs_ID, n_bO_ID, n_bOH_ID, &
        & nfiles, tot_frame, tot_atom, &
        & nx0, ny0, nz0, nzdd, &
        & nbincrdf, lensolv2drdfarr, &
        & lensolv3drdfarr, maxdtdynsteps, nbinsdyn, &
        & firstsolv2drdf, lastsolv2drdf, firstsolv3drdf, &
        & lastsolv3drdf, dbincrdf, dbinrdf, dr_3drdf, &
        & dr_rdf, dbindynwidth, dbincdynwidth, ldist, box0, &
        & dx, dy, dz)

        implicit none

        integer, intent(in) :: n_O_ID, n_OH_ID, n_Os_ID, &
            & n_subs_ID, n_bOH_ID, n_bO_ID

        integer, intent(in) :: aid, nfiles, frame

        character(LEN=1), intent(in), &
         &  dimension(1:200) :: filenames,t_t_outname

        character(LEN=1), dimension(1:200) :: t_filename,t_outname

        integer, intent(in), dimension(0:n_O_ID-1) :: O_ID

        integer, intent(in), &
        & dimension(0:n_O_ID-1) :: HO_ID1,HO_ID2

        integer, intent(in), dimension(0:n_OH_ID-1) :: OH_ID

        integer, intent(in), &
        & dimension(0:n_OH_ID-1) :: HOH_ID1,HOH_ID2


        integer, intent(in), dimension(0:n_bOH_ID-1) :: bOH_ID

        integer, intent(in), &
        & dimension(0:n_bOH_ID-1) :: HbOH_ID1,HbOH_ID2


        integer, intent(in), &
        & dimension(0:n_subs_ID-1) :: subs_ID

        integer, intent(in), &
        & dimension(0:n_Os_ID-1) :: Os_ID

        integer, intent(in), &
        & dimension(0:n_bO_ID-1) :: bO_ID


        real(4), intent(in) :: dt

        character(LEN=200) :: f_filename

        character(LEN=200) :: outname

        !real(4), &
        !& dimension(0:2*(n_OH_ID+n_Os_ID)-1) :: occ,occkappa

        integer, intent(in) :: tot_frame, &
        & tot_atom, &
        & nx0, ny0, nz0, nzdd, &
        & nbincrdf, lensolv2drdfarr, &
        & lensolv3drdfarr, maxdtdynsteps, nbinsdyn, &
        & ldist
        real(4), intent(in) :: firstsolv2drdf, &
        & lastsolv2drdf, firstsolv3drdf, lastsolv3drdf, &
        & dbincrdf, dbinrdf, dr_3drdf, &
        & dr_rdf, dbindynwidth, dbincdynwidth

        real(8), intent(in), dimension(0:2) :: box0

        real(4), intent(in) :: dx, dy, dz

        real(4) :: dxt, dyt, dzt

        real(4) :: M_H2O=29910.8567, &
        & costheta=0.0, &
        & qt=0.0, dr2=0.0, &
        & const1=0.0, const2=0.0

        real(8), dimension(0:5) :: box5

        real(8), dimension(0:2) :: box,hbox,binv,hbox0,binv0

        real(4), dimension(0:2) :: dr,dr1,dra

        real(4), dimension(0:n_O_ID-1,0:3) :: neighbor

        integer, dimension(0:n_O_ID-1,0:3) :: neighborID

        integer :: i, j, k, n, n1, ja, m, o, m1, tnpairs, &
        & nx, ny, nz

        real(4) :: dix=0.0, diy=0.0, diz=0.0        

        real :: T1, T2

        real(4), dimension(0:tot_atom-1,0:2) :: r

        real(4), dimension(0:nbincrdf-1) :: nsolv2drdfarr

        real(4) :: nsolv3drdfarr,pi=4.0*atan(1.0)

        real(4), dimension(0:nz0-1) :: qz, costhz, Oz

        real(4), dimension(0:nz0-1,0:ldist-1) :: costhdistz, qdistz

        real(4), dimension(0:nbincrdf-1,0:lensolv2drdfarr-1) :: solv2drdfarr

        real(4), dimension(0:lensolv3drdfarr-1) :: solv3drdfarr

        real(4), dimension(0:nbinsdyn-1,0:maxdtdynsteps-1) :: dr_dt_rc, n_dt_rc, &
        & dr_dt_2dmsd, n_dt_2dmsd, n_dt_RT

        real(4), dimension(0:nz0-1) :: npairs_z_wdwa, &
        & npairs_z_wdsOa, npairs_z_wdsOHa, &
        & npairs_z_wdbOa, npairs_z_wdbOHa, &
        & npairs_z_sOHdsOHa, npairs_z_sOHdsOa, &
        & npairs_z_sOHdwa, &
        & npairs_z_sOHdbOa, npairs_z_sOHdbOHa, &
        & npairs_z_bOHdwa, npairs_z_bOHdsOHa, npairs_z_bOHdsOa 

        integer, dimension(0:nz0-1) :: nwd, nsd, nbd

        !real(4), dimension(0:nx0-1,0:ny0-1,0:nz0-1) :: Oxyz, costhxyz, tO, Okappa

        print *, "Testing inputs"
        print *, SIZE(O_ID)," == ",n_O_ID
        print *, SIZE(OH_ID)," == ",n_OH_ID
        print *, SIZE(HO_ID1)," == ",n_O_ID
        print *, SIZE(HO_ID2)," == ",n_O_ID
        print *, SIZE(HOH_ID1)," == ",n_OH_ID
        print *, SIZE(HOH_ID2)," == ",n_OH_ID
        print *, SIZE(Os_ID)," == ",n_Os_ID

        print *, filenames

        call convert_c_string_f_string(t_t_outname,200,outname)
        print *, filenames
        call convert_c_string_f_string(filenames,200,f_filename)
        print *, "f_filename is ",f_filename
        print *, "Size of 2drdf is ",size(solv2drdfarr), &
      & " = ",nbincrdf*lensolv2drdfarr
        print *, "Shape of 2drdf is ",shape(solv2drdfarr), &
      & " = ",nbincrdf,"x",lensolv2drdfarr

        !########### Initialize ###########
        solv2drdfarr(:,:) = 0.0
        nsolv2drdfarr(:) = 0.0
        solv3drdfarr(:) = 0.0
        nsolv3drdfarr = 0.0


        !########### Initialize ###########
        hbox0(:) = 0.5*box0(:)
        binv0(:) = 1.0/box0(:)

        !########### Initialize ###########
        dr_dt_rc(:,:) = 0
        n_dt_rc(:,:) = 0
        dr_dt_2dmsd(:,:) = 0
        n_dt_2dmsd(:,:) = 0
        n_dt_RT(:,:) = 0
        print *, "Size of dynamic parameters is ", &
        & shape(dr_dt_rc)," = ", &
        & nbinsdyn," X ",maxdtdynsteps

        !########### Initialize ############
        print *, "Grid created of cell width : ",dx,dy,dz
        print *, "Length of grid box : ",nx0,ny0,nz0
        print *, "Box length of grid : ",box0
        print *, "Half box : ",hbox0

        print *, "Allocated"
        costhz(:) = 0.0
        qz(:) = 0.0
        Oz(:) = 0.0
        !Oxyz(:,:,:) = 0.0
        !costhxyz(:,:,:) = 0.0
        !tO(:,:,:) = 0.0
        !Okappa(:,:,:) = 0.0
        r(:,:) = 0.0      
        qdistz(:,:) = 0.0
        costhdistz(:,:) = 0.0
        npairs_z_wdwa(:) = 0.0
        npairs_z_wdsOa(:) = 0.0
        npairs_z_wdsOHa(:) = 0.0
        npairs_z_sOHdsOa(:) = 0.0
        npairs_z_sOHdsOHa(:) = 0.0
        npairs_z_sOHdwa(:) = 0.0

        npairs_z_wdbOa(:) = 0.0
        npairs_z_wdbOHa(:) = 0.0
        npairs_z_sOHdbOa(:) = 0.0
        npairs_z_sOHdbOHa(:) = 0.0
        npairs_z_bOHdsOa(:) = 0.0
        npairs_z_bOHdsOHa(:) = 0.0
        npairs_z_bOHdwa(:) = 0.0

       
        nwd(:) = 0
        nsd(:) = 0
        nbd(:) = 0
        tnpairs = 0

        print *, "Initialized"
        
        neighborID(:,:) = 0
        neighbor(:,:) = 0.0
        box5(:) = 0.0
        box(:) = 0.0
        hbox(:) = 0.0
        binv(:) = 0.0
        dr(:) = 0.0
        dr1(:) = 0.0
        dra(:) = 0.0
        print *, "Initialized"
        n = 0
        do n1 = 1, nfiles
        
        call read_header(tot_atom, tot_frame, box, nx, ny, nz, dxt, dyt, dzt, f_filename)
        nx = nx+1
        ny = ny+1
        nz = nz+1
        dix = 1.0/dx
        diy = 1.0/dy
        diz = 1.0/dz
        print *, "Opening file",f_filename
        open(unit = 5200, file = trim(f_filename), form = "UNFORMATTED", status = "OLD")

        read(5200) 
        read(5200)
        read(5200)
        do i = 1, frame-1
        read(5200)
        read(5200)
        read(5200)
        read(5200)
        end do
        print *, "Skipped frames "
        call cpu_time(T1)
        do m = frame, tot_frame - maxdtdynsteps - 10
        !call alldyn_z(dr_dt_rc, n_dt_rc, dr_dt_2dmsd, n_dt_2dmsd, n_dt_RT, &
        !& f_filename, outname, O_ID, dbindynwidth, dbincdynwidth, &
        !& m, aid, tot_frame, tot_atom, maxdtdynsteps, nbinsdyn, n_O_ID)
        r(:,:) = 0.0 
        n = n+1     !!!! Keep count of frames read
        read(5200) box5
        read(5200) r(:,0)
        read(5200) r(:,1)
        read(5200) r(:,2)
        box(0) = box5(0)
        box(1) = box5(2)
        box(2) = box5(5)
        binv(:) = 1.0/box(:)
        hbox(:) = box(:)*0.5
        !print *, "Re-initialized"
        call pbcwrap(r, r(aid-1,:), tot_atom, box, binv)
        !print *, "Wrapped"
        !!!!!!!!!!!!!!!!!!!!!! Call functions that are frame property

        !call densfluc(occ, occkappa, [OH_ID,Os_ID], n_OH_ID+n_Os_ID, &
        !& O_ID, n_O_ID, r, tot_atom, box, binv)
        !print *, "outside Dens. fluc."

        !call findfournn(neighborID, O_ID, n_O_ID, [O_ID,OH_ID], &
        !& n_O_ID+n_OH_ID, r, tot_atom, box, binv)

        !print *, "Found NN"
        !call Odensity(Orz,r,O_ID,subs_ID,box,binv,dz,nzdd,n_O_ID,n_subs_ID,tot_atom)
        !!!!!! Convert code to work for 2D rdf, shell is a cylinder, 
        !!!!!! distance is 2D distance. nideal is Npairs in slab/Vol of slab
        !!!!!! Plot a heat map of 2d rdf (r2D(x-axis) vs z(y-axis)) as it changes away from interface
        !!!!! Write a separate code for hbonding coordination number, using distance criteria b/w Os, OHs etc
        !call twodrdfsolv_pf(nsolv2drdfarr, solv2drdfarr, &
        !   & O_ID, r, dbinrdf, dbincrdf, &
        !   & lastsolv2drdf, firstsolv2drdf, &
        !   & box, binv, hbox0, n_O_ID, tot_atom, &
        !   & lensolv2drdfarr, nbincrdf)
        !print *, "Outside 2D"
        !call threedrdfsolv_pf(nsolv3drdfarr, solv3drdfarr, O_ID, &
        !   & O_ID, r, lastsolv3drdf, &
        !   & firstsolv3drdf, box, binv, hbox0, &
        !   & n_O_ID, n_O_ID, tot_atom, &
        !   & lensolv3drdfarr)
        !print *, "Outside 3D"
        !if (m == tot_frame .AND. n1 == nfiles-1) then
        !call cpu_time(T2)
        !print *, "Before vol density calc: frames", &
        !       & m-frame, "of file", f_filename, &
        !       & "was computed in time", (T2-T1)
        !end if

        do j = 0, n_bOH_ID-1
        ja = bOH_ID(j)-1
        i = nint((r(ja,2))*1/dz) + (nz0+1)/2-1

        tnpairs=0
        call dynnc_bond_angle(tnpairs, [bOH_ID(j)], 1, &
        & Os_ID, n_Os_ID, r, "dist: angle:", &
        & tot_atom, box, binv, 12.25, &
        & COS(30.0*pi/180.0), &
        & reshape([HbOH_ID1(j),HbOH_ID2(j)],[1,2]))
        npairs_z_bOHdsOa(i) = npairs_z_bOHdsOa(i) + tnpairs

        tnpairs=0
        call dynnc_bond_angle(tnpairs, [bOH_ID(j)], 1, &
        & OH_ID, n_OH_ID, r, "dist: angle:", &
        & tot_atom, box, binv, 12.25, &
        & COS(30.0*pi/180.0), &
        & reshape([HbOH_ID1(j),HbOH_ID2(j)],[1,2]))
        npairs_z_bOHdsOHa(i) = npairs_z_bOHdsOHa(i) + 2.0*tnpairs

        tnpairs=0
        call dynnc_bond_angle(tnpairs, [bOH_ID(j)], 1, &
        & O_ID, n_O_ID, r, "dist: angle:", &
        & tot_atom, box, binv, 12.25, &
        & COS(30.0*pi/180.0), &
        & reshape([HbOH_ID1(j),HbOH_ID2(j)],[1,2]))
        npairs_z_bOHdwa(i) = npairs_z_bOHdwa(i) + tnpairs

        nbd(i) = nbd(i)+1
        end do


        do j = 0, n_OH_ID-1
        ja = OH_ID(j)-1
        i = nint((r(ja,2))*1/dz) + (nz0+1)/2-1

        tnpairs=0
        call dynnc_bond_angle(tnpairs, [OH_ID(j)], 1, &
        & Os_ID, n_Os_ID, r, "dist: angle:", &
        & tot_atom, box, binv, 12.25, &
        & COS(30.0*pi/180.0), &
        & reshape([HOH_ID1(j),HOH_ID2(j)],[1,2]))
        npairs_z_sOHdsOa(i) = npairs_z_sOHdsOa(i) + tnpairs

        tnpairs=0
        call dynnc_bond_angle(tnpairs, [OH_ID(j)], 1, &
        & OH_ID, n_OH_ID, r, "dist: angle:", &
        & tot_atom, box, binv, 12.25, &
        & COS(30.0*pi/180.0), &
        & reshape([HOH_ID1(j),HOH_ID2(j)],[1,2]))
        npairs_z_sOHdsOHa(i) = npairs_z_sOHdsOHa(i) + 2.0*tnpairs

        tnpairs=0
        call dynnc_bond_angle(tnpairs, [OH_ID(j)], 1, &
        & O_ID, n_O_ID, r, "dist: angle:", &
        & tot_atom, box, binv, 12.25, &
        & COS(30.0*pi/180.0), &
        & reshape([HOH_ID1(j),HOH_ID2(j)],[1,2]))
        npairs_z_sOHdwa(i) = npairs_z_sOHdwa(i) + tnpairs

        tnpairs=0
        call dynnc_bond_angle(tnpairs, [OH_ID(j)], 1, &
        & bO_ID, n_bO_ID, r, "dist: angle:", &
        & tot_atom, box, binv, 12.25, &
        & COS(30.0*pi/180.0), &
        & reshape([HOH_ID1(j),HOH_ID2(j)],[1,2]))
        npairs_z_sOHdbOa(i) = npairs_z_sOHdbOa(i) + tnpairs

        tnpairs=0
        call dynnc_bond_angle(tnpairs, [OH_ID(j)], 1, &
        & bOH_ID, n_bOH_ID, r, "dist: angle:", &
        & tot_atom, box, binv, 12.25, &
        & COS(30.0*pi/180.0), &
        & reshape([HOH_ID1(j),HOH_ID2(j)],[1,2]))
        npairs_z_sOHdbOHa(i) = npairs_z_sOHdbOHa(i) + tnpairs



        nsd(i) = nsd(i)+1
        end do

        !!!!!!!!!!!!!!!! Loop over water oxygens
        !tO(:,:,:) = 0.0
        do j = 0, n_O_ID-1
        !!!!!!!!!!!!!!!! Call functions that are atom property
        ja = O_ID(j)-1
        i = nint((r(ja,2))*1/dz) + (nz0+1)/2-1
        o = nint((r(ja,1)+hbox0(1))*1/dy)
        m1 = nint((r(ja,0)+hbox0(0))*1/dx)

        !if (i>=0 .AND. o>=0 .AND. m1>=0 .AND. i<=nz0-1 .AND. o<=ny0-1 .AND. m1<=nx0-1) then
        !!!!!!!!!!!!!!!! Hydrogen bond information !!!!!!!!!!!!!
        tnpairs=0
        call dynnc_bond_angle(tnpairs, [O_ID(j)], 1, &
        & O_ID, n_O_ID, r, "dist: angle:", &
        & tot_atom, box, binv, 12.25, &
        & COS(30.0*pi/180.0), &
        & reshape([HO_ID1(j),HO_ID2(j)],[1,2]))
        npairs_z_wdwa(i) = npairs_z_wdwa(i) + 2.0*tnpairs
       ! !print *, tnpairs
        tnpairs=0
        call dynnc_bond_angle(tnpairs, [O_ID(j)], 1, &
        & OH_ID, n_OH_ID, r, "dist: angle:", &
        & tot_atom, box, binv, 12.25, &
        & COS(30.0*pi/180.0), &
        & reshape([HO_ID1(j),HO_ID2(j)],[1,2]))
        npairs_z_wdsOHa(i) = npairs_z_wdsOHa(i) + tnpairs
       ! !print *, tnpairs
        tnpairs=0
        call dynnc_bond_angle(tnpairs, [O_ID(j)], 1, &
        & Os_ID, n_Os_ID, r, "dist: angle:", &
        & tot_atom, box, binv, 12.25, &
        & COS(30.0*pi/180.0), &
        & reshape([HO_ID1(j),HO_ID2(j)],[1,2]))
        npairs_z_wdsOa(i) = npairs_z_wdsOa(i) + tnpairs

        tnpairs=0
        call dynnc_bond_angle(tnpairs, [O_ID(j)], 1, &
        & bO_ID, n_bO_ID, r, "dist: angle:", &
        & tot_atom, box, binv, 12.25, &
        & COS(30.0*pi/180.0), &
        & reshape([HO_ID1(j),HO_ID2(j)],[1,2]))
        npairs_z_wdbOa(i) = npairs_z_wdbOa(i) + tnpairs

        tnpairs=0
        call dynnc_bond_angle(tnpairs, [O_ID(j)], 1, &
        & bOH_ID, n_bOH_ID, r, "dist: angle:", &
        & tot_atom, box, binv, 12.25, &
        & COS(30.0*pi/180.0), &
        & reshape([HO_ID1(j),HO_ID2(j)],[1,2]))
        npairs_z_wdbOHa(i) = npairs_z_wdbOHa(i) + tnpairs


        nwd(i) = nwd(i)+1
        !!!!!!!!!!!!!!!!  Dipole orientation wrt z

        costheta = 0.0
        call waterdipole(costheta, r(ja,:), &
        & r(ja+1,:), r(ja+2,:), box, binv)
        if (costheta < -1 .OR. costheta > 1) then
        print *, "Unusual dipole orientation &
        & wrt surface normal calculated for", &
        & O_ID(j), "=", costheta
        end if
        !costhxyz(m1,o,i) = costhxyz(m1,o,i) + costheta
        !!!!!!!!!!!!!!!! costheta varies between -1 and 1
        costhz(i) = costhz(i) + costheta
        costhdistz(i, nint((costheta - (-1))*(ldist-1)/2.0)) = &
        & costhdistz(i,nint((costheta - (-1))*(ldist-1)/2.0)) + 1.0

        !!!!!!!!!!!!!!!! qfour tetrahedral order parameter

        !qt = 0.0
        !call qfourcalc(qt, neighborID, j, ja, n_O_ID, &
        !& [O_ID,OH_ID], n_O_ID+n_OH_ID, r, &
        !& tot_atom, box, binv)
        !if (qt < -3 .OR. qt > 1) then
        !print *, "Unusual tetrahedral parameter &
        !& calculated for", O_ID(j), "=", qt
        !end if
        !qz(i) = qz(i) + qt

        !!!!!!!!!!!!!!!! Distribution of qfour. qfour varies between -3 and 1
        !qdistz(i, nint((qt - (-3))*(ldist-1.0)/4.0)) = &
        !& qdistz(i,nint((qt - (-3))*(ldist-1.0)/4.0)) + 1.0

        !!!!!!!!!!!!!!!! Number density

        Oz(i) = Oz(i)+1.0
        !tO(m1,o,i) = tO(m1,o,i) + 1.0
        !Oxyz(m1,o,i) = Oxyz(m1,o,i)+1.0
        !end if
        end do
        !!!!!!!!!!!!!!! End atom loop

        call cpu_time(T2)
        !print *, outname
        open(199, file=trim(outname)//"progress.txt", &
        & STATUS="REPLACE",ACTION="WRITE")
        write(199,FMT=*) "Progress of frames: ", &
        & int(float(m-frame)/float(tot_frame-frame)*100), &
        & "% done ", &
        & int((T2-T1)/float(m-frame)*float(tot_frame-m)/60.0), &
        & " min left", &
        & int((T2-T1)/float(m-frame)*float(tot_frame-m)/60.0/60.0), &
        & " hrs left"

        close(199)

        !Okappa(:,:,:) = Okappa(:,:,:) + tO(:,:,:)*tO(:,:,:)

        end do
        !!!!!!!!!!!!!!!! End frame loop

        close(5200)
        end do
        !!!!!!!!!!!!!!! End file loop
        print *, "Closing calculations; &
        & Starting normalization"

        do i = 0, nz0-1
        !do j = 0, ny0-1
        !do k = 0, nx0-1
        !if (Oxyz(k,j,i) .eq. 0.0) then
        !costhxyz(k,j,i) = -9999.0
        !else
        !costhxyz(k,j,i) = &
        !& costhxyz(k,j,i)/Oxyz(k,j,i)
        !end if
        !end do
        !end do

        if (Oz(i) .eq. 0.0) then
        !qz(i) = -9999
        costhz(i) = -9999
        else
        !qz(i) = qz(i)/Oz(i)
        costhz(i) = costhz(i)/Oz(i)
        end if

        if (nwd(i) .eq. 0) then
        npairs_z_wdwa(i) = -9999
        npairs_z_wdsOHa(i) = -9999
        npairs_z_wdsOa(i) = -9999
        npairs_z_wdbOa(i) = -9999
        npairs_z_wdbOHa(i) = -9999
        else
        npairs_z_wdwa(i) = &
           & npairs_z_wdwa(i)/nwd(i)
        npairs_z_wdsOHa(i) = &
           & npairs_z_wdsOHa(i)/nwd(i)
        npairs_z_wdsOa(i) = &
           & npairs_z_wdsOa(i)/nwd(i)
        npairs_z_wdbOa(i) = &
           & npairs_z_wdbOa(i)/nbd(i)
        npairs_z_wdbOHa(i) = &
           & npairs_z_wdbOHa(i)/nbd(i)
        end if

        if (nsd(i) .eq. 0) then
        npairs_z_sOHdsOa(i) = -9999
        npairs_z_sOHdsOHa(i) = -9999
        npairs_z_sOHdwa(i) = -9999
        npairs_z_sOHdbOa(i) = -9999
        npairs_z_sOHdbOHa(i) = -9999
        else
        npairs_z_sOHdsOa(i) = &
           & npairs_z_sOHdsOa(i)/nsd(i)
        npairs_z_sOHdsOHa(i) = &
           & npairs_z_sOHdsOHa(i)/nsd(i)
        npairs_z_sOHdwa(i) = &
           & npairs_z_sOHdwa(i)/nsd(i)
        npairs_z_sOHdbOa(i) = &
           & npairs_z_sOHdbOa(i)/nbd(i)
        npairs_z_sOHdbOHa(i) = &
           & npairs_z_sOHdbOHa(i)/nbd(i)
        end if

        if (nbd(i) .eq. 0) then
        npairs_z_bOHdsOa(i) = -9999
        npairs_z_bOHdsOHa(i) = -9999
        npairs_z_bOHdwa(i) = -9999
        else
        npairs_z_bOHdsOa(i) = &
           & npairs_z_bOHdsOa(i)/nbd(i)
        npairs_z_bOHdsOHa(i) = &
           & npairs_z_bOHdsOHa(i)/nbd(i)
        npairs_z_bOHdwa(i) = &
           & npairs_z_bOHdwa(i)/nbd(i)
        end if


        end do
       ! print *, "w",npairs_z_wdwa(:)+npairs_z_wdsOa(:)+npairs_z_wdsOHa(:)
        !print *, COS(30.0*pi/180.0)
        print *, "Normalization 1 done"

        const2 = 1.0/n
        !Okappa(:,:,:) = ( Okappa(:,:,:)/Oxyz(:,:,:) ) - (Oxyz(:,:,:)*const2)
        const1 = M_H2O*const2*binv(0)*binv(1)*diz
        Oz(:) = Oz(:)*const1
        const1 = M_H2O*const2*dix*diy*diz
        !Oxyz(:,:,:) = Oxyz(:,:,:)*const1
        print *, "Normalization 2 done"

        !occ(:) = occ(:)*const2
        !occkappa(:) = occkappa(:)*const2
        !occkappa(:) = (occkappa(:) - occ(:)*occ(:))/(occ(:))

!! Derivation of normalization for 2d rdf:
!! For j = 0,lensolv2drdfarr-1
!! Area element for j = 0 is pi*(( 0.5 )^2 - 0)*drdf^2
!! --> Area element        = pi*(0.5*0.5)*(drdf*drdf)
!! Area element for j > 0 is pi*(( j + 0.5 )^2 - ( j - 0.5 )^2)*drdf^2
!! --> Area element        = pi*((1)*(2j))*drdf^2
!! --> Area element        = pi*(2)*(drdf*drdf)*(j)
        const2 = 1.0/(pi*0.25*(dr_rdf*dr_rdf))
        const1 = 1.0/(pi*2.0*(dr_rdf*dr_rdf))

        !do i = 0,nbincrdf-1
        !if (nsolv2drdfarr(i)==0.0) then
        !solv2drdfarr(i,:) = -9999.0
        !else
        !do j = 0,lensolv2drdfarr-1
        !if (j == 0) then
        !solv2drdfarr(i,j) = &
        !   & (solv2drdfarr(i,j)*box(0)*box(1)*const2)/(nsolv2drdfarr(i))
        !else
        !solv2drdfarr(i,j) = &
        !   & (solv2drdfarr(i,j)*box(0)*box(1)*const1)/(nsolv2drdfarr(i)*(j))
        !end if
        !end do
        !end if
        !end do

!! Derivation of normalization for 3d rdf:
!! For j = 0,lensolv3drdfarr-1
!! Vol element for j = 0 is 4.0/3.0*pi*(( 0.5 )^3 - 0)*drdf^3
!! --> Vol element        = 4.0/3.0*pi*(0.5*0.5*0.5)*(drdf*drdf*drdf)
!! Vol element for j > 0 is 4.0/3.0*pi*(( j + 0.5 )^3 - ( j - 0.5 )^3)*drdf^3
!! --> Vol element        = 4.0/3.0*pi*((1)*(3*j^2 + 0.25))*drdf^3
!!
        const2 = 1.0/((4.0/3.0)*pi*0.5*0.5*0.5*&
        &(dr_3drdf*dr_3drdf*dr_3drdf))
        const1 = 1.0/((4.0/3.0)*pi*&
        &(dr_3drdf*dr_3drdf*dr_3drdf))

        !do j = 0,lensolv3drdfarr-1
        !if (j == 0) then
        !solv3drdfarr(j) = (solv3drdfarr(j)*box(0)*box(1)*&
        !&box(2)*const2)/(nsolv3drdfarr)
        !else
        !solv3drdfarr(j) = (solv3drdfarr(j)*box(0)*box(1)*&
        !&box(2)*const1)/(nsolv3drdfarr*(3*j*j + 0.25))
        !end if
        !end do

        !open(2001, file=trim(outname)//"Oxyz.dat", access="stream")
        !write(2001) Oxyz
        !close(2001)

        !open(200, file=trim(outname)//"Oxyz_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) nx0,ny0,nz0
        !close(200)

        !open(200, file=trim(outname)//"Oxyz_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) -hbox0(0),-hbox0(0)+(nx0-1)*dx,-hbox0(1),-hbox0(1)+(ny0-1)*dy,-hbox0(2),-hbox0(2)+(nz0-1)*dz
        !close(200)

        !open(2001, file=trim(outname)//"costhxyz.dat", access="stream")
        !write(2001) costhxyz
        !close(2001)

        !open(200, file=trim(outname)//"costhxyz_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) nx0,ny0,nz0
        !close(200)

        !open(200, file=trim(outname)//"costhxyz_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) -hbox0(0),-hbox0(0)+(nx0-1)*dx,-hbox0(1),-hbox0(1)+(ny0-1)*dy,-hbox0(2),-hbox0(2)+(nz0-1)*dz
        !close(200)

        open(2001, file=trim(outname)//"costhz.dat", access="stream")
        write(2001) costhz
        close(2001)

        open(200, file=trim(outname)//"costhz_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)

        open(200, file=trim(outname)//"costhz_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        !open(2001, file=trim(outname)//"qz.dat", access="stream")
        !write(2001) qz
        !close(2001)
        !open(200, file=trim(outname)//"qz_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) nz0
        !close(200)
        !open(200, file=trim(outname)//"qz_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        !close(200)

        !open(2001, file=trim(outname)//"qdistz.dat", access="stream")
        !write(2001) qdistz
        !close(2001)

        !open(200, file=trim(outname)//"qdistz_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) nz0,ldist
        !close(200)

        !open(200, file=trim(outname)//"qdistz_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz,-3.0,1.0
        !close(200)

        open(2001, file=trim(outname)//"Oz.dat", access="stream")
        write(2001) Oz
        close(2001)

        open(200, file=trim(outname)//"Oz_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)

        open(200, file=trim(outname)//"Oz_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        !open(2001, file=trim(outname)//"Okappa.dat", access="stream")
        !write(2001) Okappa
        !close(2001)

        !open(200, file=trim(outname)//"Okappa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) nx0,ny0,nz0
        !close(200)

        !open(200, file=trim(outname)//"Okappa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) -hbox0(0),-hbox0(0)+(nx0-1)*dx,-hbox0(1),-hbox0(1)+(ny0-1)*dy,-hbox0(2),-hbox0(2)+(nz0-1)*dz
        !close(200)

        !open(2001, file=trim(outname)//"occ.dat", access="stream")
        !write(2001) occ
        !close(2001)

        !open(200, file=trim(outname)//"occ_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) 2*(n_OH_ID+n_Os_ID)
        !close(200)

        !open(200, file=trim(outname)//"occ_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) 1,2*(n_OH_ID+n_Os_ID)
        !close(200)

        !open(2001, file=trim(outname)//"occkappa.dat", access="stream")
        !write(2001) occkappa
        !close(2001)

        !open(200, file=trim(outname)//"occkappa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) 2*(n_OH_ID+n_Os_ID)
        !close(200)

        !open(200, file=trim(outname)//"occkappa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) 1,2*(n_OH_ID+n_Os_ID)
        !close(200)

        !open(2001, file=trim(outname)//"zrdf.dat", access="stream")
        !write(2001) solv2drdfarr
        !close(2001)

        !open(200, file=trim(outname)//"zrdf_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) nbincrdf,lensolv2drdfarr
        !close(200)

        !open(200, file=trim(outname)//"zrdf_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) -hbox0(2),-hbox0(2)+(nbincrdf-1)*dz,0.0,(lensolv2drdfarr-1)*dr_rdf
        !close(200)

        !open(2001, file=trim(outname)//"3drdf.dat", access="stream")
        !write(2001) solv3drdfarr
        !close(2001)

        !open(200, file=trim(outname)//"3drdf_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) lensolv3drdfarr
        !close(200)

        !open(200, file=trim(outname)//"3drdf_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) 0.0,(lensolv3drdfarr-1)*dr_3drdf
        !close(200)

        open(2001, file=trim(outname)//"costhdistz.dat", access="stream")
        write(2001) costhdistz
        close(2001)

        open(200, file=trim(outname)//"costhdistz_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0,ldist
        close(200)

        open(200, file=trim(outname)//"costhdistz_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz,-1.0,1.0
        close(200)

        !npairs_z_wdwa, npairs_z_wdsOa, npairs_z_wdsOHa, npairs_z_sOHdsOHa, npairs_z_sOHdsOa, npairs_z_sOHdwa
        !print *, npairs_z_wdwa

        open(2001, file=trim(outname)//"wdbOa.dat", access="stream")
        write(2001) npairs_z_wdbOa
        close(2001)
        open(200, file=trim(outname)//"wdbOa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)
        open(200, file=trim(outname)//"wdbOa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        open(2001, file=trim(outname)//"wdbOHa.dat", access="stream")
        write(2001) npairs_z_wdbOHa
        close(2001)
        open(200, file=trim(outname)//"wdbOHa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)
        open(200, file=trim(outname)//"wdbOHa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        open(2001, file=trim(outname)//"sOHdbOa.dat", access="stream")
        write(2001) npairs_z_sOHdbOa
        close(2001)
        open(200, file=trim(outname)//"sOHdbOa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)
        open(200, file=trim(outname)//"sOHdbOa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        open(2001, file=trim(outname)//"sOHdbOHa.dat", access="stream")
        write(2001) npairs_z_sOHdbOHa
        close(2001)
        open(200, file=trim(outname)//"sOHdbOHa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)
        open(200, file=trim(outname)//"sOHdbOHa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        open(2001, file=trim(outname)//"bOHdwa.dat", access="stream")
        write(2001) npairs_z_bOHdwa
        close(2001)
        open(200, file=trim(outname)//"bOHdwa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)
        open(200, file=trim(outname)//"sOHdwa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        open(2001, file=trim(outname)//"bOHdsOa.dat", access="stream")
        write(2001) npairs_z_bOHdsOa
        close(2001)
        open(200, file=trim(outname)//"bOHdsOa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)
        open(200, file=trim(outname)//"bOHdsOa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        open(2001, file=trim(outname)//"bOHdsOHa.dat", access="stream")
        write(2001) npairs_z_bOHdsOHa
        close(2001)
        open(200, file=trim(outname)//"bOHdsOHa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)
        open(200, file=trim(outname)//"bOHdsOHa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)



        open(2001, file=trim(outname)//"wdwa.dat", access="stream")
        write(2001) npairs_z_wdwa
        close(2001)

        open(200, file=trim(outname)//"wdwa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)

        open(200, file=trim(outname)//"wdwa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        open(2001, file=trim(outname)//"wdsOa.dat", access="stream")
        write(2001) npairs_z_wdsOa
        close(2001)

        open(200, file=trim(outname)//"wdsOa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)

        open(200, file=trim(outname)//"wdsOa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        open(2001, file=trim(outname)//"wdsOHa.dat", access="stream")
        write(2001) npairs_z_wdsOHa
        close(2001)

        open(200, file=trim(outname)//"wdsOHa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)

        open(200, file=trim(outname)//"wdsOHa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        open(2001, file=trim(outname)//"sOHdsOHa.dat", access="stream")
        write(2001) npairs_z_sOHdsOHa
        close(2001)

        open(200, file=trim(outname)//"sOHdsOHa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)

        open(200, file=trim(outname)//"sOHdsOHa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        open(2001, file=trim(outname)//"sOHdsOa.dat", access="stream")
        write(2001) npairs_z_sOHdsOa
        close(2001)

        open(200, file=trim(outname)//"sOHdsOa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)

        open(200, file=trim(outname)//"sOHdsOa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        open(2001, file=trim(outname)//"sOHdwa.dat", access="stream")
        write(2001) npairs_z_sOHdwa
        close(2001)

        open(200, file=trim(outname)//"sOHdwa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)

        open(200, file=trim(outname)//"sOHdwa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        !do m = frame, tot_frame - maxdtdynsteps - 10
        !print *, "At frame = ",m," , Last frame",tot_frame-maxdtdynsteps-10
        !call alldyn_z(dr_dt_rc, n_dt_rc, dr_dt_2dmsd, n_dt_2dmsd, n_dt_RT, &
        !& f_filename, outname, O_ID, dbindynwidth, dbincdynwidth, &
        !& m, aid, tot_frame, tot_atom, maxdtdynsteps, nbinsdyn, n_O_ID)
        !end do 

        !dr_dt_rc,n_dt_rc,dr_dt_2dmsd,n_dt_2dmsd,n_dt_RT
        !do i = 0,nbinsdyn-1
        !do j = 0,maxdtdynsteps-1
        !dr_dt_rc(i,j) = dr_dt_rc(i,j)/max(1.0,n_dt_rc(i,j))
        !dr_dt_2dmsd(i,j) = dr_dt_2dmsd(i,j)/max(1.0,n_dt_2dmsd(i,j))
        !end do
        !dix = n_dt_RT(i,0)
        !n_dt_RT(i,:) = n_dt_RT(i,:)/dix
        !end do

        !open(2001, file=trim(outname)//"rc_dr.dat", access="stream")
        !write(2001) dr_dt_rc
        !close(2001)

        !open(200, file=trim(outname)//"rc_dr_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) nbinsdyn,maxdtdynsteps
        !close(200)

        !open(200, file=trim(outname)//"rc_dr_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz,dt,maxdtdynsteps*dt
        !close(200)

        !open(2001, file=trim(outname)//"rc_n.dat", access="stream")
        !write(2001) n_dt_rc
        !close(2001)

        !open(200, file=trim(outname)//"rc_n_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) nbinsdyn,maxdtdynsteps
        !close(200)

        !open(200, file=trim(outname)//"rc_n_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz,dt,maxdtdynsteps*dt
        !close(200)

        !open(2001, file=trim(outname)//"2dmsd_dr.dat", access="stream")
        !write(2001) dr_dt_2dmsd
        !close(2001)

        !open(200, file=trim(outname)//"2dmsd_dr_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) nbinsdyn,maxdtdynsteps
        !close(200)

        !open(200, file=trim(outname)//"2dmsd_dr_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz,dt,maxdtdynsteps*dt
        !close(200)

        !open(2001, file=trim(outname)//"2dmsd_n.dat", access="stream")
        !write(2001) n_dt_2dmsd
        !close(2001)

        !open(200, file=trim(outname)//"2dmsd_n_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) nbinsdyn,maxdtdynsteps
        !close(200)

        !open(200, file=trim(outname)//"2dmsd_n_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz,dt,maxdtdynsteps*dt
        !close(200)

        !open(2001, file=trim(outname)//"nRT.dat", access="stream")
        !write(2001) n_dt_RT
        !close(2001)

        !open(200, file=trim(outname)//"nRT_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) nbinsdyn,maxdtdynsteps
        !close(200)

        !open(200, file=trim(outname)//"nRT_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz,dt,maxdtdynsteps*dt
        !close(200)
end subroutine structure_z

subroutine structure(aid,t_t_outname,O_ID,HO_ID1,HO_ID2,OH_ID,HOH_ID1,HOH_ID2,HbOH_ID1,HbOH_ID2,Os_ID,subs_ID,bO_ID,bOH_ID,filenames,dt,nfiles,frame,n_O_ID,n_OH_ID,n_Os_ID,n_subs_ID,n_bO_ID,n_bOH_ID) 
        implicit none

        integer, intent(in) :: n_O_ID, n_OH_ID, n_Os_ID, n_subs_ID, n_bO_ID, n_bOH_ID

        integer, intent(in) :: aid, nfiles, frame

        character(LEN=1), intent(in), &
         &  dimension(1:200) :: filenames,t_t_outname


        character(LEN=1), dimension(1:200) :: t_filename,t_outname

        integer, intent(in), dimension(0:n_O_ID-1) :: O_ID

        integer, intent(in), &
        & dimension(0:n_O_ID-1) :: HO_ID1,HO_ID2

        integer, intent(in), dimension(0:n_OH_ID-1) :: OH_ID

        integer, intent(in), &
        & dimension(0:n_OH_ID-1) :: HOH_ID1,HOH_ID2

        integer, intent(in), dimension(0:n_bOH_ID-1) :: bOH_ID

        integer, intent(in), &
        & dimension(0:n_bOH_ID-1) :: HbOH_ID1,HbOH_ID2

        integer, intent(in), &
        & dimension(0:n_subs_ID-1) :: subs_ID

        integer, intent(in), &
        & dimension(0:n_Os_ID-1) :: Os_ID

        integer, intent(in), &
        & dimension(0:n_bO_ID-1) :: bO_ID

        real(4), intent(in) :: dt

        character(LEN=200) :: f_filename

        character(LEN=200) :: outname

        !real(4), &
        !& dimension(0:2*(n_OH_ID+n_Os_ID)-1) :: occ,occkappa

        integer :: tot_frame, &
        & tot_atom, &
        & nx, ny, nz , nx0, ny0, nz0, nzdd, &
        & nbincrdf, lensolv2drdfarr, &
        & lensolv3drdfarr

        real(4) :: dx, dy, dz, dxt, dyt, dzt

        real(4) :: M_H2O=29910.8567, &
        & costheta=0.0, &
        & qt=0.0, dr2=0.0, &
        & const1=0.0, const2=0.0

        real(8), dimension(0:5) :: box5

        real(8), dimension(0:2) :: box,hbox,binv,hbox0,box0,binv0

        real(4), dimension(0:2) :: dr,dr1,dra

        !real(4), dimension(0:n_O_ID-1,0:3) :: neighbor

        !integer, dimension(0:n_O_ID-1,0:3) :: neighborID

        integer :: i, j, k, n, n1, ja, m, o, m1, tnpairs, &
        & maxdtdynsteps, nbinsdyn, ldist=50

        real(4) :: dix=0.0, diy=0.0, diz=0.0, firstsolv2drdf, &
        & lastsolv2drdf, firstsolv3drdf, lastsolv3drdf, &
        & dbincrdf, dbinrdf, pi=4.0*atan(1.0), dr_3drdf, &
        & dr_rdf, dbindynwidth, dbincdynwidth

        real :: T1, T2
        t_filename(:) = filenames(:)
        print *, "HOH_ID1"
        print *, HOH_ID1
        print *, "HOH_ID2"
        print *, HOH_ID2
        print *, "O_ID"
        print *, O_ID
        print *, t_filename
        call convert_c_string_f_string(t_filename,200,f_filename)
        print *, "f_filename is ",f_filename

        call read_header(tot_atom, tot_frame, box, &
        & nx0, ny0, nz0, dx, dy, dz, &
        & f_filename)
        !nx0 = nx0+1
        !ny0 = ny0+1
        !nz0 = nz0+1

        !########### Control parameters ###
        nbincrdf = nz0
        dbincrdf = dz
        dbinrdf = 0.5 !Angstrom

        !########### Initialize ###########
        firstsolv2drdf = 0.0
        dr_rdf = 0.2
        lensolv2drdfarr = floor(0.5*min(dx*nx0,dy*ny0)/dr_rdf)
        lastsolv2drdf = lensolv2drdfarr*(dr_rdf)
        print *, nbincrdf*lensolv2drdfarr
        print *, nbincrdf,"x",lensolv2drdfarr

        !########### Initialize ###########
        firstsolv3drdf = 0.0
        dr_3drdf = 0.2
        lensolv3drdfarr = floor(0.5*min(dx*nx0,dy*ny0)/dr_3drdf)
        lastsolv3drdf = lensolv3drdfarr*(dr_3drdf)

        !########### Initialize ###########
        !nzdd = nz0/2
        !nz0 = 2*(int(nz0+int(6.0/dz))/2+1)+1 
        box0(:) = box(:)
        !box0(2) = (nz0-1)*dz
        hbox0(:) = 0.5*box0(:)
        binv0(:) = 1.0/box0(:)

        !########### Initialize ###########
        dbindynwidth = 3.2
        dbincdynwidth = dz
        maxdtdynsteps = int(1500.0/dt)
        nbinsdyn = nz0 !int(box0(2)/dbincdynwidth+1)

        print *, nbinsdyn," X ",maxdtdynsteps

        !########### Initialize ############
        print *, "Grid created of cell width : ",dx,dy,dz
        print *, "Length of grid box : ",nx0,ny0,nz0
        print *, "Box length of grid : ",box0
        print *, "Half box : ",hbox0
        print *, "Nested function started"
        print *, "Testing inputs"
        print *, SIZE(O_ID)," == ",n_O_ID
        print *, SIZE(OH_ID)," == ",n_OH_ID
        call structure_z(aid, t_t_outname, O_ID, HO_ID1, HO_ID2, OH_ID, &
           & HOH_ID1, HOH_ID2, HbOH_ID1, HbOH_ID2, &
           & Os_ID, subs_ID, bO_ID, bOH_ID, filenames, dt, frame, &
           & n_O_ID, n_OH_ID, n_Os_ID, n_subs_ID, n_bO_ID, n_bOH_ID, &
           & nfiles, tot_frame, tot_atom, &
           & nx0, ny0, nz0, nzdd, &
           & nbincrdf, lensolv2drdfarr, &
           & lensolv3drdfarr, maxdtdynsteps, nbinsdyn, &
           & firstsolv2drdf, lastsolv2drdf, firstsolv3drdf, &
           & lastsolv3drdf, dbincrdf, dbinrdf, dr_3drdf, &
           & dr_rdf, dbindynwidth, dbincdynwidth, ldist, box0, &
           & dx, dy, dz)
 
end subroutine structure
