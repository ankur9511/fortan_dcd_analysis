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

subroutine dyncal_z(aid, t_t_outname, O_ID, HO_ID1, HO_ID2, OH_ID, &
        & HOH_ID1, HOH_ID2, &
        & Os_ID, subs_ID, filenames, dt, frame, &
        & n_O_ID, n_OH_ID, n_Os_ID, n_subs_ID, &
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

        real(4), &
        & dimension(0:2*(n_OH_ID+n_Os_ID)-1) :: occ,occkappa

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

        real(4), dimension(0:nz0-1) :: qz, costhz, Oz, &
           & tOz, Okappaz

        real(4), dimension(0:nz0-1,0:ldist-1) :: costhdistz, qdistz

        real(4), dimension(0:nbincrdf-1,0:lensolv2drdfarr-1) :: solv2drdfarr

        real(4), dimension(0:lensolv3drdfarr-1) :: solv3drdfarr

        real(4), dimension(0:nbinsdyn-1,0:maxdtdynsteps-1) :: dr_dt_rc, n_dt_rc, &
        & dr_dt_2dmsd, n_dt_2dmsd, n_dt_RT

        real(4), dimension(0:nz0-1) :: npairs_z_wdwa, &
        & npairs_z_wdsOa, npairs_z_wdsOHa, &
        & npairs_z_sOHdsOHa, npairs_z_sOHdsOa, &
        & npairs_z_sOHdwa

        integer, dimension(0:nz0-1) :: nwd, nsd

        real(4), dimension(0:nx0-1,0:ny0-1,0:nz0-1) :: Oxyz, costhxyz
        print *, "Testing inputs"
        print *, SIZE(O_ID)," == ",n_O_ID
        print *, SIZE(OH_ID)," == ",n_OH_ID
        print *, SIZE(HO_ID1)," == ",n_O_ID
        print *, SIZE(HO_ID2)," == ",n_O_ID
        print *, SIZE(HOH_ID1)," == ",n_OH_ID
        print *, SIZE(HOH_ID2)," == ",n_OH_ID
        print *, SIZE(Os_ID)," == ",n_Os_ID

        print *, filenames

        !t_outname(:) = t_t_outname(:)
        call convert_c_string_f_string(t_t_outname,200,outname)
        occ(:) = 0.0
        occkappa(:) = 0.0
        print *, filenames
        !t_filename(:) = filenames(:)
        !print *, t_filename
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
        Oxyz(:,:,:) = 0.0
        costhxyz(:,:,:) = 0.0
        tOz(:) = 0.0
        Okappaz(:) = 0.0
        r(:,:) = 0.0      
        qdistz(:,:) = 0.0
        costhdistz(:,:) = 0.0
        npairs_z_wdwa(:) = 0.0
        npairs_z_wdsOa(:) = 0.0
        npairs_z_wdsOHa(:) = 0.0
        npairs_z_sOHdsOa(:) = 0.0
        npairs_z_sOHdsOHa(:) = 0.0
        npairs_z_sOHdwa(:) = 0.0
        nwd(:) = 0
        nsd(:) = 0
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
        !print *, nfiles
        n = 0
        do n1 = 1, nfiles
        
        do m = frame, tot_frame - maxdtdynsteps - 10
        !print *, "At frame = ",m," , Last frame",tot_frame-maxdtdynsteps-10
        call alldyn_z(dr_dt_rc, n_dt_rc, dr_dt_2dmsd, n_dt_2dmsd, n_dt_RT, &
     & f_filename, outname, O_ID, dbindynwidth, dbincdynwidth, &
     & m, aid, tot_frame, tot_atom, maxdtdynsteps, nbinsdyn, n_O_ID)
        end do 
        
        end do
        !dr_dt_rc,n_dt_rc,dr_dt_2dmsd,n_dt_2dmsd,n_dt_RT
        do i = 0,nbinsdyn-1
        do j = 0,maxdtdynsteps-1
        dr_dt_rc(i,j) = dr_dt_rc(i,j)/max(1.0,n_dt_rc(i,j))
        dr_dt_2dmsd(i,j) = dr_dt_2dmsd(i,j)/max(1.0,n_dt_2dmsd(i,j))
        end do
        dix = n_dt_RT(i,0)
        n_dt_RT(i,:) = n_dt_RT(i,:)/dix
        end do

        open(2001, file=trim(outname)//"rc_dr.dat", access="stream")
        write(2001) dr_dt_rc
        close(2001)

        open(200, file=trim(outname)//"rc_dr_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nbinsdyn,maxdtdynsteps
        close(200)

        open(200, file=trim(outname)//"rc_dr_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz,dt,maxdtdynsteps*dt
        close(200)

        !open(2001, file=trim(outname)//"rc_n.dat", access="stream")
        !write(2001) n_dt_rc
        !close(2001)

        !open(200, file=trim(outname)//"rc_n_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) nbinsdyn,maxdtdynsteps
        !close(200)

        !open(200, file=trim(outname)//"rc_n_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz,dt,maxdtdynsteps*dt
        !close(200)

        open(2001, file=trim(outname)//"2dmsd_dr.dat", access="stream")
        write(2001) dr_dt_2dmsd
        close(2001)

        open(200, file=trim(outname)//"2dmsd_dr_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nbinsdyn,maxdtdynsteps
        close(200)

        open(200, file=trim(outname)//"2dmsd_dr_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz,dt,maxdtdynsteps*dt
        close(200)

        !open(2001, file=trim(outname)//"2dmsd_n.dat", access="stream")
        !write(2001) n_dt_2dmsd
        !close(2001)

        !open(200, file=trim(outname)//"2dmsd_n_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) nbinsdyn,maxdtdynsteps
        !close(200)

        !open(200, file=trim(outname)//"2dmsd_n_range.dat",STATUS="REPLACE",ACTION="WRITE")
        !write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz,dt,maxdtdynsteps*dt
        !close(200)

        open(2001, file=trim(outname)//"nRT.dat", access="stream")
        write(2001) n_dt_RT
        close(2001)

        open(200, file=trim(outname)//"nRT_dim.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) nbinsdyn,maxdtdynsteps
        close(200)

        open(200, file=trim(outname)//"nRT_range.dat",STATUS="REPLACE",ACTION="WRITE")
        write(200,FMT=*) -hbox0(2),-hbox0(2)+(nz0-1)*dz,dt,maxdtdynsteps*dt
        close(200)
end subroutine dyncal_z

subroutine dyncal(aid,t_t_outname,O_ID,HO_ID1,HO_ID2,OH_ID,HOH_ID1,HOH_ID2,Os_ID,subs_ID,filenames,dt,nfiles,frame,n_O_ID,n_OH_ID,n_Os_ID,n_subs_ID) 
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

        real(4), &
        & dimension(0:2*(n_OH_ID+n_Os_ID)-1) :: occ,occkappa

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

        real(4), dimension(0:n_O_ID-1,0:3) :: neighbor

        integer, dimension(0:n_O_ID-1,0:3) :: neighborID

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
        nx0 = nx0+1
        ny0 = ny0+1
        nz0 = nz0+1

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
        nzdd = nz0/2
        nz0 = 2*(int(nz0+int(6.0/dz))/2+1)+1 
        box0(:) = box(:)
        box0(2) = (nz0-1)*dz
        hbox0(:) = 0.5*box0(:)
        binv0(:) = 1.0/box0(:)

        !########### Initialize ###########
        dbindynwidth = 3.3
        dbincdynwidth = dz
        maxdtdynsteps = int(200.0/dt)
        nbinsdyn = int(box0(2)/dbincdynwidth+1)

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
        call dyncal_z(aid, t_t_outname, O_ID, HO_ID1, HO_ID2, OH_ID, &
           & HOH_ID1, HOH_ID2, &
           & Os_ID, subs_ID, filenames, dt, frame, &
           & n_O_ID, n_OH_ID, n_Os_ID, n_subs_ID, &
           & nfiles, tot_frame, tot_atom, &
           & nx0, ny0, nz0, nzdd, &
           & nbincrdf, lensolv2drdfarr, &
           & lensolv3drdfarr, maxdtdynsteps, nbinsdyn, &
           & firstsolv2drdf, lastsolv2drdf, firstsolv3drdf, &
           & lastsolv3drdf, dbincrdf, dbinrdf, dr_3drdf, &
           & dr_rdf, dbindynwidth, dbincdynwidth, ldist, box0, &
           & dx, dy, dz)
 
end subroutine dyncal
