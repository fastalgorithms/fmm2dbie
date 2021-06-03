      implicit real *8 (a-h,o-z)

      ntests =1 
      call test_dchunkints(i1)

      ntests = ntests + 1
      call test_zchunkints(i2)

      nsuc = i1 + i2

      open(unit=33,file='../../print_testres.txt',access='append')
      write(33,'(a,i2,a,i2,a)') 'Successfully completed ',nsuc,
     1  ' out of ',ntests,'tests in chunk_routs testing suite'
      close(33)


      stop
      end
