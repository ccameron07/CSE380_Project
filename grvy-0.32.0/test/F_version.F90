!! -*-f90-*-
!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!! 
!! libGRVY - a utility library for scientific computing.
!!
!! Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
!!
!! This library is free software; you can redistribute it and/or
!! modify it under the terms of the Version 2.1 GNU Lesser General
!! Public License as published by the Free Software Foundation.
!!
!! This library is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!! Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public
!! License along with this library; if not, write to the Free Software
!! Foundation, Inc. 51 Franklin Street, Fifth Floor, 
!! Boston, MA  02110-1301  USA
!!
!!-----------------------------------------------------------------------el-
!!
!! Regression/Test Suite for libGRVY.
!!
!! $Id: F_version.F90.in 29241 2012-04-05 18:44:56Z karl $
!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------

#define GRVY_MAJOR_VERSION  0
#define GRVY_MINOR_VERSION  32
#define GRVY_MICRO_VERSION  0

program main
  use grvy
  implicit none

  integer :: grvy_version
  integer :: expected_version

  expected_version = GRVY_MAJOR_VERSION*1000 + GRVY_MINOR_VERSION*100 + GRVY_MICRO_VERSION;
  grvy_version     = grvy_get_numeric_version()

  if(grvy_version .ne. expected_version)then
     call exit(1)
  endif

  ! Verify stdout call is cooperating

  call grvy_log_setlevel(GRVY_NOLOG)

  call grvy_version_stdout()

  stop
end program main





