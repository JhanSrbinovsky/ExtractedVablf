! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Subroutine: U_MODEL -----------------------------------------------
!
!    Purpose: High level control program for the Unified Model
!             (master routine).  Calls lower level control routines
!             according to top level switch settings. Called by
!             top level routine UMSHELL which provides dimension sizes
!             for dynamic allocation of data arrays.
!
!    Programming standard: UM Doc Paper 3, version 8.3
!
!    External documentation: On-line UM document C0 - The top-level
!                            control system
!
!    -------------------------------------------------------------------

!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Top Level

SUBROUTINE U_MODEL(                                           &
   NFT,NFTU,                                                  &
#include "argszsp.h"
#include "argszspa.h"
#include "argszspc.h"
    )

! D1 replacement module
  USE atm_fields_mod

! OASIS Modules 
  USE oasis_atm_data_mod


! For deleting superseded dumps
  USE MPPIO_job_control,        ONLY : jobCntrl     
  USE MPPIO_job_control_common, ONLY : JC_delete
  USE filenamelength_mod,       ONLY : filenamelength

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  USE io, ONLY             : io_timestep
  USE model_file, ONLY : storeAllLookups, SynchroniseAll
  USE um_types

  USE ereport_mod, ONLY : ereport
  USE printstatus_mod
  USE um_parvars
  USE control_max_sizes
  USE decomp_db

  ! ensure that module variables cannot go out of scope
  USE wet_to_dry_n_calc_mod, ONLY : wet_to_dry_n

  USE um_input_control_mod,  ONLY:                     &
  l_oasis, oasis_couple_freq
  USE river_inputs_mod, ONLY: l_rivers
  USE ppxlook_mod
  USE acp_namel_mod, ONLY: l_ac
  USE ukca_option_mod, ONLY: l_ukca

  USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim
  USE Submodel_Mod

  USE nlstcall_mod, ONLY : model_basis_time, &
                           LPP, &
                           ldump, &
                           lmean, &
                           lprint, &
                           linterface, &
                           lexit, &
                           ljobrelease, &
                           lancillary, &
                           lboundary, &
                           ltimer, &
                           ft_select

  USE chsunits_mod, ONLY : nunits

  use cable_data_mod, ONLY : set_endstep_umodel

  IMPLICIT NONE

!    Interface and arguments: ------------------------------------------
!       Sizes of super arrays
#include "typszsp.h"
#include "typszspa.h"
#include "typszspc.h"

!       Model sizes
#include "typsize.h"

!       Addresses of component arrays within super arrays
#include "spindex.h"

!*----------------------------------------------------------------------

!  Common blocks
#include "chistory.h"
#include "ccontrol.h"
#include "ctime.h"
#include "cintfa.h"
#include "c_global.h"

! Data kind parameters
!  DYNAMIC ALLOCATION OF SUPER ARRAYS:
!       Main D1 data array
#include "typspd1.h"
#if defined(NECSX6)
!dir$ cache_align spd1
#endif

! STASH related arrays
#include "typspst.h"

! Dump headers and lookups
#include "typspdua.h"

! Pointers (addresses) of model variables and constants
#include "typsppta.h"
! Maximum sizes of fields limited by User Interface
! CMAXSIZE now included earlier in routine

! Model derived constants arrays
#include "typspcoa.h"

! Generation of output interface fields
#include "typspina.h"

! Updating of model from ancillary files
#include "typspana.h"

! Boundary updating for Limited Area Models
#include "typspbo.h"
#include "typspboa.h"

! Coupled model arrays (atmosphere-river routing)
#include "typspcpl.h"

! Sizes for allocation of AC scheme arrays
#include "csizeobs.h"
  INTEGER :: obs_flag_len,obs_len
  INTEGER, ALLOCATABLE, DIMENSION(:) :: obs_flag
  REAL,    ALLOCATABLE, DIMENSION(:) :: obs

! Local variables

  INTEGER internal_model    ! Work - Internal model identifier
  INTEGER internal_model_prev!Work - Previous internal model ident
  INTEGER submodel          ! Work - Submodel id for dump partition
  INTEGER ngroup            ! Work - Number of steps in "group"
  INTEGER meanlev           ! Work - Mean level indicator
  INTEGER iabort            ! Work - Internal return code
  INTEGER i_step            ! Work - Loop counter over timesteps

! Number of dust bins in LBC generation, disabled if not using MakeBC
  INTEGER :: ndustbin_in, ndustbin_out

! River routing
  INTEGER g_river_field_size   ! Sizes for MPP dynamic allocation
  INTEGER g_theta_field_size   ! for river routing fields

  LOGICAL lexitnow          ! Work - Immediate exit indicator
  INTEGER nft           ! Unit no. for standard STASHmaster files
  INTEGER nftu          ! Do. user STASH files (for GET_FILE)
  INTEGER rownumber     ! Row no. counter for PPXI, PPXC arrays
  INTEGER i,j,k,idx     ! Loop counters
! CO2 array dimensions
  INTEGER co2_dima,                                            &
     co2_dimo,                                                 &
     co2_dimo2
! DMS array dimensions
  INTEGER dms_dima,                                            &
     dms_dimo,                                                 &
     dms_dimo2
  INTEGER info   ! Return code from GCom routines

  INTEGER len_runid       !  No of chars in RUNID

  CHARACTER(LEN=4) runtype_char  !  Run Type (ie. NRUN, CRUN)

! 3-D fields of species to be passed down to radiation
  INTEGER, PARAMETER :: ngrgas = 8
  INTEGER, SAVE :: grgas_addr(ngrgas)

  LOGICAL :: put_step ! True when we're going to do a put
                      ! for OASIS purposes, false if we use dummy data.
  LOGICAL :: get_step ! True if we need to get data from OASIS, false
                      ! if we perform a get but discard results.
  LOGICAL :: cpl_update_step ! True when we need to make sure
                             ! the D1 prognostic data is updated
                             ! in preparation for dump creation or
                             ! coupling actions.

#if defined(IBM)
  INTEGER                      :: pos               ! string position
  INTEGER                      :: icode_env         ! error code
  CHARACTER (LEN=256)          :: c_xlfrteopts      ! XLFRTEOPTS content
  CHARACTER (LEN=*), PARAMETER :: c_disable='buffering=disable_all'
#endif

! Error reporting
  INTEGER       icode       ! =0 normal exit; >0 error exit
  CHARACTER(LEN=256) cmessage    ! Error message
  CHARACTER(LEN=*) routinename

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  PARAMETER (   routinename='U_MODEL')

  IF (lhook) CALL dr_hook('U_MODEL',zhook_in,zhook_handle)


  ALLOCATE(ppxi(ppxrecs,ppxref_codelen))
  ALLOCATE(ppxptr(n_internal_model,0:ppxref_sections ,ppxref_items))


  icode=0
  cmessage=''

!----------------------------------------------------------------------
! 0. Start Timer call for U_MODEL (NB: not conditional on LTIMER)

! DEPENDS ON: timer
  IF (ltimer) CALL timer('U_MODEL ',5)

! DEPENDS ON: lbc_coup_setup
  CALL lbc_coup_setup(ltimer)

!  Routine GETPPX_PART reads those ppxref records which correspond to
!  entries in the stash list into the ppx look-up arrays PPXI, PPXC.
!  It also sets the ppx pointer array PPXPTR. The lengths of PPXI, PPXC
!  have been dynamically allocated to the value of ppxRecs.

! Initialise row number in PPXI, PPXC arrays
  rownumber = 1

! Initialise lookup and pointer array
  DO i=1,ppxrecs
    DO j=1,ppxref_codelen
      ppxi(i,j)=0
    END DO
    DO j=1,ppxref_charlen
      ppxc(i,j) = ' '
    END DO
  END DO
  DO i = 1,n_internal_model
    DO j   = 0,ppxref_sections
      DO k = 1,ppxref_items
        ppxptr(i,j,k)=0
      END DO
    END DO
  END DO

! Read in STASHmaster records
  IF (internal_model_index(a_im) >  0) THEN
! DEPENDS ON: getppx_part
    CALL getppx_part(nft,nftu,'STASHmaster_A',a_im,rownumber,         &
    icode,cmessage)

    IF (icode  /=  0) CALL ereport(routinename,icode,cmessage)
  END IF

!----------------------------------------------------------------------
! 1. General initialisation of control and physical data blocks
  icode=0
! DEPENDS ON: timer
      CALL timer('INITIAL ',5)
! DEPENDS ON: initial
  CALL initial(                                                   &
#include "argszsp.h"
#include "argszspa.h"
#include "argszspc.h"
     ,                                                            &
#include "argsp.h"
#include "argspa.h"
#include "argspc.h"

     ngrgas,grgas_addr,                                           &
     internal_model,submodel,ngroup,meanlev)

  IF (icode  /=  0) CALL ereport(routinename,icode,cmessage)
! DEPENDS ON: timer
      CALL timer('INITIAL ',6)

! DEPENDS ON: icopenoutput
  CALL icopenoutput(runtype_char)

!  Allocate AC scheme arrays using sizes from AC_INIT
  IF (l_ac) THEN
    obs_flag_len = a_max_no_obs
    ALLOCATE (obs_flag(obs_flag_len))
! 2048 gives enough space in WORK in RDOBS2 when no or very few obs.
    obs_len = a_max_obs_size+2048
    ALLOCATE (obs(obs_len))
    WRITE(6,'(a,i7,i7)')'U_MODEL - OBS arrays allocated with sizes ',   &
       a_max_no_obs,a_max_obs_size+2048
  ELSE
    a_max_no_obs = 1
    a_max_obs_size = 1
    obs_flag_len = a_max_no_obs
    ALLOCATE (obs_flag(obs_flag_len))
    obs_len = a_max_obs_size
    ALLOCATE (obs(obs_len))
    WRITE(6,*)'U_MODEL - OBS arrays allocated with length 1'
  END IF

!----------------------------------------------------------------------
! 2. Check for nothing-to-do

! DEPENDS ON: exitchek
  CALL exitchek( internal_model, lexitnow)
  IF (lexitnow) GO TO 9999

!----------------------------------------------------------------------
! 3. Start group of timesteps

  IF (l_oasis) THEN

    put_step=.FALSE.
    get_step=.FALSE.
    cpl_update_step=.FALSE.
  
  END IF

  WRITE(6,'(a,f10.2,a)') 'Model running with timestep ',           &
     secs_per_stepim(a_im),' seconds'

#if defined(IBM)
! On IBM turn off buffering of Fortran I/O once initialisation complete
! provided the disabling of buffering has been asked for. This is done
! by reading  the XLFRTEOPTS environment variable and checking is for the
! 'buffering=disable_all' substring (note case and space requirements)
  CALL fort_get_env('XLFRTEOPTS', 10, c_xlfrteopts, 256, icode_env)

! If ICODE_ENV is 0, the env var has been read
  IF (icode_env == 0) THEN
    pos = INDEX(c_xlfrteopts, c_disable)

    IF (pos /= 0) THEN    ! Substring found
      IF (printstatus  >=  prstatus_oper .AND. mype == 0) THEN
        WRITE (6,*) 'Disabling Fortran I/O Buffering'
      END IF

      CALL setrteopts(c_disable)
    END IF
  END IF
#endif
!----------------------------------------------------------------------
! 3. Start group of timesteps
1 CONTINUE
!----------------------------------------------------------------------
! 3.1. Start main timestep loop

! 3.1.1 Increment model time ..
! DEPENDS ON: incrtime
  CALL incrtime (                                                   &
#include "artduma.h"
     internal_model,icode,cmessage)

  IF (icode  /=  0) CALL ereport(routinename,icode,cmessage)

! Keep tabs on PRISM PUT/GET timesteps.
! At the moment we say a put and get timestep are one and the same.
! We're also hard coding to do operations every 48th timestep i.e.
! assuming a half hour atmos TS!.
! Non-ENDGAME code
  IF (l_oasis) THEN

! Is this a genuine exchange timestep.
    put_step = (MOD(stepim(atmos_im),oasis_couple_ts).eq.1)
    get_step = (MOD(stepim(atmos_im),oasis_couple_ts).eq.1)
    cpl_update_step=(MOD(stepim(atmos_im),oasis_couple_ts).eq.0)

! Perform coupling exchanges relating to TS-1

! DEPENDS ON: oasis3_geto2a
    CALL oasis3_geto2a(                                         &
#include "artd1.h"
#include "arg_atm_fields.h"
         get_step)

! DEPENDS ON: oasis3_puta2o
    CALL oasis3_puta2o(                                         &
#include "artd1.h"
#include "arg_atm_fields.h"
         put_step)

! Advance date ready for next timestep (if there is one)
! DEPENDS ON: OASIS3_ADVANCE_DATE
    CALL oasis3_advance_date()

  END IF ! l_oasis=true
! End of Non-ENDGAME code

! 3.1.2 .. set timestep control switches
! DEPENDS ON: settsctl
  CALL settsctl (                                                   &
#include "artduma.h"
#include "artsts.h"
#include "artinfa.h"
     internal_model,.FALSE.,meanlev,icode,cmessage)

  IF (icode  /=  0) CALL ereport(routinename,icode,cmessage)
! 3.1.3 If PPfile initialisation time call PP control routine
!          for instantaneous data (MEANLEV=0)
  IF (lpp) THEN
! DEPENDS ON: ppctl_reinit
    CALL ppctl_reinit(                                              &
#include "artduma.h"
#include "artinfa.h"
       internal_model,icode,cmessage)

    IF (icode  /=  0) CALL ereport(routinename,icode,cmessage)
  END IF

! Send an EOT message to IOS
! He may or may not act on this to purge outstanding items
  CALL io_timestep()


!       Integrate atmosphere or ocean by 1 timestep
  IF (internal_model == atmos_im) THEN
! Synchronize before the timestep starts
    CALL gc_gsync(nproc,info)

! River routing
    IF (l_rivers) THEN

! Get 'global' atmos horizontal domain sizes from database
! in DECOMPDB to set dynamic allocation in ATM_STEP for River routing
! on PE 0.

      g_theta_field_size=                                              &
         decompdb(decomp_standard_atmos)%glsize(1,fld_type_p) *        &
         decompdb(decomp_standard_atmos)%glsize(2,fld_type_p)
      g_river_field_size=                                              &
         decompdb(decomp_standard_atmos)%glsize(1,fld_type_r) *        &
         decompdb(decomp_standard_atmos)%glsize(2,fld_type_r)
    ELSE
      g_theta_field_size=1
      g_river_field_size=1
    END IF
!CABLE:
   call set_endstep_umodel( TARGET_END_STEPim(a_im) )

! DEPENDS ON: atm_step
    CALL atm_step (                                                &
#include "artd1.h"
#include "artduma.h"
#include "artlndm.h"
#include "arg_atm_fields.h"
#include "artbnd.h"
#include "artsts.h"
! River routing
#include "artatcpl.h"
       g_theta_field_size,                                              &
       g_river_field_size,                                              &
       obs_flag,obs,obs_flag_len,obs_len,                               &
       ngrgas,grgas_addr)

    IF (l_ukca) THEN
! DEPENDS ON: timer
      IF (ltimer) CALL timer('UKCA_MAIN1',5)
! DEPENDS ON: ukca_main1
      CALL ukca_main1(                                                &
#include "artd1.h"
#include "artduma.h"
#include "artsts.h"
#include "artlndm.h"
#include "artptra.h"
      i)                  ! dummy to terminate call
! DEPENDS ON: timer
      IF (ltimer) CALL timer('UKCA_MAIN1',6)
    END IF

! Generate Atmosphere lateral boundary values

    IF (linterface) THEN

! DEPENDS ON: timer
      IF (ltimer) CALL timer ('GEN_INTF',3)
      ndustbin_in  = 6
      ndustbin_out = 6
! DEPENDS ON: gen_intf
      CALL gen_intf (                                                 &
#include "artd1.h"
#include "artduma.h"
#include "artsts.h"
#include "artptra.h"
#include "artinfa.h"
         submodel, ndustbin_in, ndustbin_out, icode, cmessage)

! DEPENDS ON: timer
      IF (ltimer) CALL timer ('GEN_INTF',4)

      IF (icode /= 0) THEN

        CALL ereport(routinename, icode ,cmessage)
      END IF

    END IF  ! LINTERFACE

  END IF        ! internal_model = atmos_im

  IF (l_oasis) THEN

! Applicable to OASIS3 and OASIS3-MCT

    IF (cpl_update_step) THEN
!---------------------------------------------------------------------
! Ensure atmos coupling data in prognostic areas of D1 are up to date.
! The logic here is that cpl_update_step will be TRUE on the timestep
! BEFORE coupling is due to take place.

! Newly generated coupling data is intercepted after ATM_STEP and
! copied to  D1 prognostics.

! On the next timestep i.e. a coupling timestep, the prognostic
! contents will be sent to the other components in the coupling process.

! If there is no subsequent timestep (i.e. if this is the last model
! timestep then the D1 contents will be written to the dump
! ready for any future restart). There is no "end-of-model"
! coupling exchange.
!---------------------------------------------------------------------

! Non-ENDGAME code
! DEPENDS ON: oasis_updatecpl
      CALL oasis_updatecpl(                                           &
#include "artd1.h"
#include "arg_atm_fields.h"
         cmessage)

! End of Non-ENDGAME code
    END IF
  END IF

! 3.1.4 If dump time, call dump control routine

  IF (ldump) THEN
    IF (ltimer) THEN

! DEPENDS ON: timer
      CALL timer('DUMPCTL',5)
    END IF
! DEPENDS ON: dumpctl
    CALL dumpctl (                                             &
#include "artd1.h"
#include "artduma.h"
#include "artlndm.h"
#include "artptra.h"
#include "artsts.h"
       submodel,meanlev,.FALSE.,'           ',0,               &
       icode,cmessage)

    IF (icode  /=  0) CALL ereport(routinename,icode,cmessage)

    IF (ltimer) THEN
! DEPENDS ON: timer
      CALL timer('DUMPCTL ',6)
    END IF
! 3.1.4.2 If atmosphere timestep recalculate prognostic data and
!         wrap-around fields using rounded off values
    IF (submodel == atmos_sm) THEN
! Non-ENDGAME code

! DEPENDS ON: resetatm
      CALL resetatm(                                            &
#include "artd1.h"
#include "artptra.h"
         icode,cmessage)

      IF (icode >  0) GO TO 9999
! End of non-ENDGAME code
    END IF
  END IF ! ldump

! 3.1.5 If printed output time, call print control routine
  IF (lprint) THEN

    IF (printstatus >= prstatus_oper) THEN
      WRITE(6,*) routinename,':Warning, Printing of climate global ' &
         ,'and zonal diagnostics no longer supported'
    END IF  ! PrintStatus test

  END IF
! 3.1.6.1 Release job to process output created so far, if selected
  IF (ljobrelease) THEN
! Write all pp files' cached lookups to ensure that current data is
! is in the files output stream
    CALL storeAllLookups()

! Flush/Sync all pp files' units to ensure that data is committed from the
! application
    CALL SynchroniseAll()

! DEPENDS ON: jobctl
    CALL jobctl(internal_model,icode,cmessage)

    IF (icode  /=  0) CALL ereport(routinename,icode,cmessage)
  END IF
! 3.1.7 If partial sum/mean creation time, call means control routine
!       (calls mean PPfield and diagnostic print routines internally)

  IF (lmean) THEN
! DEPENDS ON: meanctl
    CALL meanctl (                                              &
#include "artd1.h"
#include "artduma.h"
#include "artptra.h"
#include "artsts.h"
#include "artlndm.h"
#include "artinfa.h"
       submodel,meanlev,icode,cmessage)
    IF (icode >  0) THEN
      CALL ereport(routinename,icode,cmessage)
    END IF
  END IF

! 3.1.8 Update history file once dumping and any meaning is complete
  IF (ldump) THEN
! checkpoint_dump_im now successfully written. Update history file to restart
! from this instead of superseded dump stored in save_dumpname. Only delete 
! superseded dump when this has been done.

! DEPENDS ON: set_history_values
    CALL set_history_values
    ! Write history to backup location before overwriting main file
    IF (mype  ==  0) THEN
! DEPENDS ON: temphist
      CALL temphist(thist_unit,icode,cmessage)
      IF (icode /= 0) THEN
        WRITE(6,*)routinename,':Failure writing temporary restart file'
        WRITE(6,*)'Check for problems and restart from main file'
        CALL ereport(routinename,icode,cmessage)
      END IF
      CALL temphist(xhist_unit,icode,cmessage)
      IF (icode /= 0) THEN
        WRITE(6,*)routinename,':Failure writing main restart file'
        WRITE(6,*)'Check for problems and restart from temporary file'
        WRITE(6,*)'by overwriting main file with temporary file'
        CALL ereport(routinename,icode,cmessage)
      ELSE
        ! Main restart file successfully written, so delete backup
! DEPENDS ON: del_hist
        CALL del_hist(thist_unit)
      END IF

      IF (ft_select(22) == "Y") THEN
        ! History files written so can delete superseded dump.
        ! Variable will be blank if last dump was timestep zero start dump
        IF (save_dumpname_im(a_im) /= blank_file_name) THEN
          CALL jobcntrl(jc_delete,save_dumpname_im(a_im))
        END IF
        save_dumpname_im(a_im)=checkpoint_dump_im(a_im)
      END IF ! ft_select(22) == "Y"
    END IF ! If (mype /= 0)
  END IF ! If (ldump)


! 3.1.9 If exit check time, check for immediate exit
  IF (lexit) THEN
! DEPENDS ON: exitchek
    CALL exitchek(internal_model, lexitnow)
    IF (lexitnow) THEN
      IF (.NOT. ldump) THEN
        WRITE(6,*)routinename,                                                 &
           ': Warning: exiting at a period that is not a dump period'
        WRITE(6,*)'Therefore continuing the run will rerun preceding timesteps'
        WRITE(6,*)'This is inefficient and can cause restart problems'
      END IF
      GO TO 9999
    END IF
  END IF ! IF (lexit)
! 3.1.10 Update ancillary fields if necessary
  IF (lancillary) THEN
! DEPENDS ON: up_ancil
    CALL up_ancil (                                   &
#include "artd1.h"
#include "artduma.h"
#include "artptra.h"
#include "artanc.h"
       submodel,                                      &
       icode,cmessage)

    IF (icode  /=  0) CALL ereport(routinename,icode,cmessage)
  END IF
! 3.1.11 Update boundary fields if necessary
  IF (lboundary) THEN
! DEPENDS ON: lbc_coup_update
    CALL lbc_coup_update(      &
#include "artd1.h"
#include "artduma.h"
#include "artptra.h"
#include "artbnd.h"
       submodel,icode,cmessage)
  END IF

!      End main timestep loop
!----------------------------------------------------------------------

  GO TO 1

9999 CONTINUE

!----------------------------------------------------------------------
! 4. Exit processing: Output error messages and perform tidy-up

  IF (l_oasis) THEN
! DEPENDS ON: oasis_tidy
    CALL oasis_tidy(                       &
#include "artd1.h"
#include "artsts.h"
#include "artduma.h"
#include "artptra.h"
       icode,cmessage)
  END IF

! DEPENDS ON: iccloseoutput
  CALL iccloseoutput()

! 4.1 Exit processing: If abnormal completion, output error message
  iabort = icode
  IF (icode /= 0) THEN

    CALL ereport(routinename,icode,cmessage)
  END IF
! 4.2 Exit processing: Perform tidy-up
! DEPENDS ON: exitproc
  CALL exitproc(icode,cmessage)
! 4.3 Exit processing: If error in exit processing, output error mess
  IF (icode /= 0) THEN

    CALL ereport(routinename,icode,cmessage)
  END IF


  DEALLOCATE(ppxi)
  DEALLOCATE(ppxptr)

!----------------------------------------------------------------------
! 5. Complete Timer call and return

  icode=iabort
! DEPENDS ON: timer
  IF (ltimer) CALL timer('U_MODEL ',6)
  IF (lhook) CALL dr_hook('U_MODEL',zhook_out,zhook_handle)

  RETURN
END SUBROUTINE u_model
