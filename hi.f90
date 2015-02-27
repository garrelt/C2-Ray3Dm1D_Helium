    count = 1
    pass = .false.
    do while (count.le.length)

      a1 = real(triangular_a_array(count)) - parity_a*half !sign is + if parity a is +
      a2 = real(triangular_a_array(count)) + parity_a*half !sign is + if parity a is +
      c1 = a1*(cc/aa)
      c2 = a2*(cc/aa)
      !c-distance test
      if (abs(c1-triangular_c_array(count)).le.half .and. abs(c2-triangular_c_array(count)).le.half ) then
        b1 = a1*(bb/aa)
        b2 = a2*(bb/aa)
        pass_b = .true.
      elseif (abs(c1-triangular_c_array(count)).le.half .and. abs(c2-triangular_c_array(count)).gt.half ) then
        b1 = a1*(bb/aa)
        b2 = (triangular_c_array(count)+parity_c*half)*(bb/cc)
        pass_b = .true.
      elseif (abs(c1-triangular_c_array(count)).gt.half .and. abs(c2-triangular_c_array(count)).le.half ) then
        b1 = (triangular_c_array(count)-parity_c*half)*(bb/cc)
        b2 = a2*(bb/aa)
        pass_b = .true.
      else
        pass_b = .false.
      endif
			
      !b-distance test
      if (pass_b.eqv..true.) then
        if (abs(b1-triangular_b_array(count)).le.half .or. abs(b2-triangular_b_array(count)).le.half) then
          pass = .true.
          passed_cell = count
          triangular_array(count) = 2
          count = length					
        endif
      endif
			
      count = count + 1
			
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !no Ifront cells are passed through
    if (pass.eqv..false.) then 
      bb = bb + parity_b*1.0 

    !An Ifront cell is passed through
    else
      a1 = real(triangular_a_array(passed_cell)) - parity_a*half 
      b1 = a1*(bb/aa)
      c1 = a1*(cc/aa)
      !front test
      if (abs(b1-triangular_b_array(passed_cell)).le.half .and. &
          abs(c1-triangular_c_array(passed_cell)).le.half) then
        test_result = 0
      endif	
      c1 = real(triangular_c_array(passed_cell)) - parity_c*half 
      a1 = c1*(aa/cc)
      b1 = c1*(bb/cc)
      !right test
      if (abs(a1-triangular_a_array(passed_cell)).le.half .and. &
          abs(b1-triangular_b_array(passed_cell)).le.half) then
        test_result = 1	
      endif	
      b1 = real(triangular_b_array(passed_cell)) - parity_b*half 
      a1 = b1*(aa/bb)
      c1 = b1*(cc/bb)
      !bottom test
      if (abs(a1-triangular_a_array(passed_cell)).le.half .and. &
          abs(c1-triangular_c_array(passed_cell)).le.half) then
        a1 = real(triangular_a_array(passed_cell)) - parity_a*half 
        c1 = a1*(cc/aa)
      if (abs(c1-triangular_c_array(passed_cell)).le.half) test_result = 0
      if (abs(c1-triangular_c_array(passed_cell)).gt.half) test_result = 1

	      !update the ray parameter
      if (test_result.eq.0) then

        bb = aa*(real(triangular_b_array(passed_cell)) + parity_b*half)/ &
             (real(triangular_a_array(passed_cell)) - parity_a*half) + parity_b*0.01

      elseif (test_result.eq.1) then

        bb = cc*(real(triangular_b_array(passed_cell)) + parity_b*half)/ &
             (real(triangular_c_array(passed_cell)) - parity_c*half) + parity_b*0.01

      endif
    endif			
