      implicit real *8 (a-h,o-z)
      ntests = 3
      call test_near_field_routs(i1)
      call test_zself_adap(i2)

      nsuc = i1+i2
      open(unit=33,file='../../print_testres.txt',access='append')
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ',nsuc,
     1  ' out of ',ntests,' tests in quadratures testing suite'
      close(33)


      return
      end
