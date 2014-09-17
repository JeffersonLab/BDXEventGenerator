      program pgs
      implicit none

c generated particle list
      include 'pgs.inc'

      nevhep = 12
      jmohep(1,2) = 11
      phep(1,1) = 15.123456789E10
      
      numtrk = 5
      ptrk(2,1) = 1.23456789e15

      hcal(2,1) = 2.34567890e15
      met_cal = pi

      pclu(2,1) = 3.456789012e15
      mulclu(2) = 123

      vectrg(9,2) = 4.567890123e15
      indtrg(2) = 234

      vecobj(10,3) = 5.678901234e15
      indobj(3) = 345

      call test_cpp(10)
      
      call pgs2root_ini
      call pgs2root_evt
      call pgs2root_end
      
      end
