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
        & HOH_ID1, HOH_ID2, &
        & Os_ID, subs_ID, filenames, dt, frame, &
        & n_O_ID, n_OH_ID, n_Os_ID, n_subs_ID, &
        & nfiles, tot_frame, tot_atom, &
        & nx0, ny0, nz0, nzdd, &
        & ldist, box0, &
        & dx, dy, dz)

        implicit none

        integer, intent(in) :: n_O_ID, n_OH_ID, n_Os_ID, &
            & n_subs_ID

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

        integer, intent(in), &
        & dimension(0:n_subs_ID-1) :: subs_ID

        integer, intent(in), &
        & dimension(0:n_Os_ID-1) :: Os_ID

        real(4), intent(in) :: dt

        character(LEN=200) :: f_filename

        character(LEN=200) :: outname

        integer, intent(in) :: tot_frame, &
        & tot_atom, &
        & nx0, ny0, nz0, nzdd, ldist

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
        & tnshell, nx, ny, nz

        real(4) :: dix=0.0, diy=0.0, diz=0.0        

        real :: T1, T2

        real(4), dimension(0:tot_atom-1,0:2) :: r

        real(4) :: nsolv3drdfarr,pi=4.0*atan(1.0)

        real(4), dimension(0:nz0-1) :: qz, costhz, tcosthz, Oz, tOz

        real(4), dimension(0:nz0-1,0:ldist-1) :: costhdistz, qdistz

        real(4), dimension(0:nz0-1) :: npairs_z_wdwa, nshell_z_wdwa !, &
        !& npairs_z_wdsOa, npairs_z_wdsOHa, &
        !& npairs_z_sOHdsOHa, npairs_z_sOHdsOa, &
        !& npairs_z_sOHdwa
        integer, dimension(0:n_O_ID-1) :: shell1
        integer, dimension(0:nz0-1) :: nwd !, nsd

        real(4), dimension(0:nx0-1,0:ny0-1,0:nz0-1) :: Oxyz, costhxyz, tO, Okappa

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
        
        !########### Initialize ###########
 
        !########### Initialize ###########
        hbox0(:) = 0.5*box0(:)
        binv0(:) = 1.0/box0(:)

        !########### Initialize ###########

        !########### Initialize ############
        print *, "Grid created of cell width : ",dx,dy,dz
        print *, "Length of grid box : ",nx0,ny0,nz0
        print *, "Box length of grid : ",box0
        print *, "Half box : ",hbox0

        print *, "Allocated"
        costhz = 0.0
        tcosthz = 0.0
        qz = 0.0
        Oz = 0.0
        tOz = 0.0
        Oxyz = 0.0
        costhxyz = 0.0
        tO = 0.0
        Okappa = 0.0
        r = 0.0      
        qdistz = 0.0
        costhdistz = 0.0
        npairs_z_wdwa = 0
        nshell_z_wdwa = 0
        !npairs_z_wdsOa = 0.0
        !npairs_z_wdsOHa = 0.0
        !npairs_z_sOHdsOa = 0.0
        !npairs_z_sOHdsOHa = 0.0
        !npairs_z_sOHdwa = 0.0
        !nwd = 0
        !nsd = 0
        !tnpairs = 0

        print *, "Initialized"
        
        neighborID = 0
        neighbor = 0.0
        box5 = 0.0
        box = 0.0
        hbox = 0.0
        binv = 0.0
        dr = 0.0
        dr1 = 0.0
        dra = 0.0
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
        do m = frame, tot_frame - 1
 
        r = 0.0 
        n = n+1     !!!! Keep count of frames read
        read(5200) box5
        read(5200) r(:,0)
        read(5200) r(:,1)
        read(5200) r(:,2)
        box(0) = box5(0)
        box(1) = box5(2)
        box(2) = box5(5)
        binv = 1.0/box
        hbox = box*0.5
  
        call pbcwrap(r, r(aid-1,:), tot_atom, box, binv)
 
        !!!!!!!!!!!!!!!!!!!!!! Call functions that are frame property

        !!!!!!!!!!!!!!!! Loop over water oxygens
        tO = 0.0
        do j = 0, n_O_ID-1
        !!!!!!!!!!!!!!!! Call functions that are atom property
        ja = O_ID(j)-1
        i = nint((r(ja,2))*1/dz) + (nz0+1)/2-1
        o = nint((r(ja,1)+hbox0(1))*1/dy)
        m1 = nint((r(ja,0)+hbox0(0))*1/dx)

        !!!!!!!!!!!!!!!! Hydrogen bond information !!!!!!!!!!!!!
        
        !-call nhbond_nshell(shell1, tnpairs, tnshell, &
        !-& O_ID(j),HO_ID1(j),HO_ID2(j), &
        !-& O_ID,HO_ID1,HO_ID2, &
        !-& r,box,binv,12.25, COS(30.0*pi/180.0), &
        !-& tot_atom,1,n_O_ID)
 
        !-nwd(i) = nwd(i)+1
        !-npairs_z_wdwa(i) = npairs_z_wdwa(i) + tnpairs
        !-nshell_z_wdwa(i) = nshell_z_wdwa(i) + tnshell
        !!!!!!!!!!!!!!!!  Dipole orientation wrt z

        costheta = 0.0
        call waterdipole(costheta, r(ja,:), &
        & r(ja+1,:), r(ja+2,:), box, binv)
        if (costheta < -1 .OR. costheta > 1) then
        print *, "Unusual dipole orientation &
        & wrt surface normal calculated for", &
        & O_ID(j), "=", costheta
        end if
        costhxyz(m1,o,i) = costhxyz(m1,o,i) + costheta
        !!!!!!!!!!!!!!!! costheta varies between -1 and 1
        costhz(i) = costhz(i) + costheta
        costhdistz(i, nint((costheta - (-1))*(ldist-1)/2.0)) = &
        & costhdistz(i,nint((costheta - (-1))*(ldist-1)/2.0)) + 1.0

        !!!!!!!!!!!!!!!! Number density

        Oz(i) = Oz(i)+1.0
        tO(m1,o,i) = tO(m1,o,i) + 1.0
        Oxyz(m1,o,i) = Oxyz(m1,o,i)+1.0
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

        Okappa = Okappa + tO*tO

        end do
        !!!!!!!!!!!!!!!! End frame loop

        close(5200)
        end do
        !!!!!!!!!!!!!!! End file loop
        print *, "Closing calculations; &
        & Starting normalization"

        do i = 0, nz0-1
        do j = 0, ny0-1
        do k = 0, nx0-1
        if (Oxyz(k,j,i) .eq. 0.0) then
        costhxyz(k,j,i) = -9999.0
        else
        costhxyz(k,j,i) = &
        & costhxyz(k,j,i)/Oxyz(k,j,i)
        end if
        end do
        end do

        !if (Oz(i) .eq. 0.0) then
        !qz(i) = -9999
        !costhz(i) = -9999
        !else
        !qz(i) = qz(i)/Oz(i)
        !costhz(i) = costhz(i)/Oz(i)
        !end if

        !-if (nwd(i) .eq. 0) then
        !-npairs_z_wdwa(i) = -9999
        !_nshell_z_wdwa(i) = -9999
        !-else
        !-npairs_z_wdwa(i) = npairs_z_wdwa(i)/nwd(i)
        !-nshell_z_wdwa(i) = nshell_z_wdwa(i)/nwd(i)
        !-end if
        end do

        do i = 6,nz0-7
        tOz(i) = Oz(i) !SUM(Oz(i-6:i+6))
        tcosthz(i) = costhz(i) !SUM(costhz(i-6:i+6))
        end do

        Oz = tOz
        costhz = tcosthz
        
        do i = 6,nz0-7
        if (Oz(i) .eq. 0.0) then
        costhz(i) = -9999
        else
        costhz(i) = costhz(i)/Oz(i)
        end if
        end do

        print *, "Normalization 1 done"

        const2 = 1.0/n
        Okappa = ( Okappa/Oxyz ) - (Oxyz*const2)

        const1 = M_H2O*const2*binv(0)*binv(1)*diz !/13
        Oz = Oz*const1
        
        const1 = M_H2O*const2*dix*diy*diz
        Oxyz = Oxyz*const1

        print *, "Normalization 2 done"

        open(2001, file=trim(outname)//"Oxyz.dat", access="stream")
        write(2001) Oxyz
        close(2001)

        open(200, file=trim(outname)//"Oxyz_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nx0,ny0,nz0
        close(200)

        open(200, file=trim(outname)//"Oxyz_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(0),-hbox0(0)+(nx0-1)*dx,-hbox0(1),-hbox0(1)+(ny0-1)*dy,-hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        open(2001, file=trim(outname)//"costhxyz.dat", access="stream")
        write(2001) costhxyz
        close(2001)

        open(200, file=trim(outname)//"costhxyz_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nx0,ny0,nz0
        close(200)

        open(200, file=trim(outname)//"costhxyz_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(0),-hbox0(0)+(nx0-1)*dx,-hbox0(1),-hbox0(1)+(ny0-1)*dy,-hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        open(2001, file=trim(outname)//"costhz.dat", access="stream")
        write(2001) costhz
        close(2001)

        open(200, file=trim(outname)//"costhz_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)

        open(200, file=trim(outname)//"costhz_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        open(2001, file=trim(outname)//"Oz.dat", access="stream")
        write(2001) Oz
        close(2001)

        open(200, file=trim(outname)//"Oz_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0
        close(200)

        open(200, file=trim(outname)//"Oz_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        open(2001, file=trim(outname)//"Okappa.dat", access="stream")
        write(2001) Okappa
        close(2001)

        open(200, file=trim(outname)//"Okappa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nx0,ny0,nz0
        close(200)

        open(200, file=trim(outname)//"Okappa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(0),-hbox0(0)+(nx0-1)*dx,-hbox0(1),-hbox0(1)+(ny0-1)*dy,-hbox0(2),-hbox0(2)+(nz0-1)*dz
        close(200)

        open(2001, file=trim(outname)//"costhdistz.dat", access="stream")
        write(2001) costhdistz
        close(2001)

        open(200, file=trim(outname)//"costhdistz_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nz0,ldist
        close(200)

        open(200, file=trim(outname)//"costhdistz_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz,-1.0,1.0
        close(200)

        !-open(2001, file=trim(outname)//"npairs_z_wdwa.dat", access="stream")
        !-write(2001) npairs_z_wdwa
        !-close(2001)
        !
        !-open(200, file=trim(outname)//"npairs_z_wdwa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !-write(200,FMT=*) nz0
        !-close(200)
        !
        !-open(200, file=trim(outname)//"npairs_z_wdwa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !-write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        !-close(200)
        !
        !-open(2001, file=trim(outname)//"nshell_z_wdwa.dat", access="stream")
        !-write(2001) nshell_z_wdwa
        !-close(2001)
        !
        !-open(200, file=trim(outname)//"nshell_z_wdwa_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !-write(200,FMT=*) nz0
        !-close(200)
        !
        !-open(200, file=trim(outname)//"nshell_z_wdwa_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !-write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz
        !-close(200)

end subroutine structure_z

subroutine structure(aid,t_t_outname,O_ID,HO_ID1,HO_ID2,OH_ID,HOH_ID1,HOH_ID2,Os_ID,subs_ID,filenames,dt,nfiles,frame,n_O_ID,n_OH_ID,n_Os_ID,n_subs_ID) 
        implicit none

        integer, intent(in) :: n_O_ID, n_OH_ID, n_Os_ID, n_subs_ID

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

        integer, intent(in), &
        & dimension(0:n_subs_ID-1) :: subs_ID

        integer, intent(in), &
        & dimension(0:n_Os_ID-1) :: Os_ID

        real(4), intent(in) :: dt

        character(LEN=200) :: f_filename

        character(LEN=200) :: outname

        !real(4), &
        !& dimension(0:2*(n_OH_ID+n_Os_ID)-1) :: occ,occkappa

        integer :: tot_frame, &
        & tot_atom, &
        & nx, ny, nz , nx0, ny0, nz0, nzdd

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

        integer :: i, j, k, n, n1, ja, m, o, m1, tnpairs, ldist=50

        real(4) :: dix=0.0, diy=0.0, diz=0.0, pi=4.0*atan(1.0)

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

        !########### Control parameters ###

        !########### Initialize ###########
        box0 = box
        hbox0 = 0.5*box0
        binv0 = 1.0/box0

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
           & HOH_ID1, HOH_ID2, &
           & Os_ID, subs_ID, filenames, dt, frame, &
           & n_O_ID, n_OH_ID, n_Os_ID, n_subs_ID, &
           & nfiles, tot_frame, tot_atom, &
           & nx0, ny0, nz0, nzdd, ldist, box0, &
           & dx, dy, dz)
 
end subroutine structure
