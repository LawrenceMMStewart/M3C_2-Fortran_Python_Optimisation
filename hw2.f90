!Lawrence Stewart - ls3914@ic.ac.uk CID=00948972
!Module for solving n-d optimization problems with Newton's method and 2-d problems
!with bracket descent. Necessary cost function details are provided in separate cost
!module.
module hw2
  use cost
  implicit none
  integer :: itermax = 1000 !maximum number of iterations used by an optimizer
  real(kind=8) :: tol=1.0e-6 !stopping criteria for optimizer
  real(kind=8), allocatable :: jpath(:), xpath(:,:) !allocatable 1-d and 2-d arrays, should contain cost and location values used during optimization iterations

  contains


  subroutine newton(xguess,xf,jf) 
    !Use Newton's method to minimize cost function, costj
    !input: xguess -- initial guess for loaction of minimum
    !output: xf -- computed location of minimum, jf -- computed minimum
    !Should also set module variables xpath and jpath appropriately
    implicit none
    real(kind=8), dimension(:), intent(in) :: xguess !do not need to explicitly specify dimension of input variable when subroutine is within a module
    real(kind=8), intent(out) :: xf(size(xguess)),jf !location of minimum, minimum cost



    !variables that we have added:
    integer :: i1 ! for looping purposes
    logical :: stop_check !When stop check becomes false we stop newtons method

    real(kind=8) :: cur_hess(size(xguess),size(xguess)) ! current hessian matrix
    real(kind=8) :: cur_grad(size(xguess)) ! current gradient vector
    real(kind=8) :: tempx(size(xguess)) !  The position of the x that we are evaluating in the while looping
    real(kind=8) :: cur_cost, prev_cost

    !allocate the variables required for lapack matrix solving
    integer :: INFO
    integer, dimension(size(xguess)) :: IPIV
    real(kind=8) :: temp_hess(size(xguess),size(xguess)) ! This is the temporary hessian for lapack
    real(kind=8) :: temp_grad(size(xguess))   ! THis is the temporary gradient for lapack
    real(kind=8) ::  h(size(xguess))     ! this will be the solution h
    real(kind=8) :: temp_xpath(itermax+1,size(xguess)), temp_jpath(itermax+1)
    


    !start for i1=1
    i1=1

    !stop_check being false means the loop will continue (see later)
    stop_check=.False.

    !put the start point in the location vector 
    temp_xpath(1,:)=xguess

    !evaluate the cost at the start point 
    !-- we have prev_cost for the previous value and cur_cost for the new value being evaluated
    call costj(xguess,prev_cost)

    !put the starting cost in the cost store:
    temp_jpath(1)=prev_cost


    !begin Netwons method:
    do while(.not. stop_check) 
      
      !note here that i1 is still at the previous value last completed -

      !we are going to use lapack for linear equation solving -
      !do not forget to include lapack in your command when running from terminal

      tempx=temp_xpath(i1,:)

      !calculate the hessian at the current point:

      call costj_hess2d(tempx,cur_hess)
      call costj_grad2d(tempx,cur_grad)
     
      !solve H(x)h=-g(x)
      temp_hess=cur_hess
      temp_grad=cur_grad

      !solve the linear system
      call dgesv(size(xguess), 1, temp_hess, size(xguess), IPIV, temp_grad, size(xguess), INFO) 
      h = -temp_grad(1:size(xguess)) !extract the solution from lapack
     

      tempx=tempx+h !update the current x value

      !put the point into the xpath:
      temp_xpath(i1+1,:)=tempx

      !exaluate the cost at the new point
      call costj(tempx,cur_cost)

      !put the cost into the cost matrix:
      temp_jpath(i1+1)=cur_cost

      !Check constraints that cause Newtons method to stop
      
      !if itermax has been reached
      if(i1+1==itermax) then
        stop_check = .True.
      end if

      !see if the tolerance has been reached
      call convergence_check(prev_cost,cur_cost,stop_check)

      if(stop_check .eqv. .True.) then
        xf=tempx
        jf=cur_cost 
        ! print *, xf, jf, i1

      end if

      i1=i1+1
      prev_cost=cur_cost

    end do


    !deallocate if already in use
    if(allocated(xpath)) then
         deallocate(xpath, jpath)
    end if
    !allocate xpath and jpath
    allocate(xpath(i1,size(xguess)),jpath(i1)) 
    xpath=temp_xpath(1:i1,:)
    jpath=temp_jpath(1:i1)

  end subroutine newton















  subroutine bracket_descent(xguess,xf,jf)
    !Use bracket descent method to minimize cost function, costj
    !input: xguess -- initial guess for loaction of minimum
    !output: xf -- computed location of minimum, jf -- computed minimum
    !Should also set module variables xpath and jpath appropriately
    !Assumes size(xguess) = 2
    implicit none
    real(kind=8), dimension(2), intent(in) :: xguess
    real(kind=8), intent(out) :: xf(2),jf !location of minimum, minimum cost

    !Below are the variables we added to the function
    logical :: conv_check ! checks to see if we have converged
    integer :: i1
    real(kind=8), dimension(2) :: xinit_1,xinit_2,xinit_3
    real(kind=8) :: ja,jb,jc,jstar,jsstar, prev_sumj,sumj,tjf ! costs
    real(kind=8), dimension(2) :: va,vb,vc,xm,xstar,xsstar
    real(kind=8) :: temp_xpath(itermax+1,size(xguess)), temp_jpath(itermax+1) !this was changed in  newton as well!!!!
    real(kind=8), dimension(2) :: centroid

    ! set a boolean variable to notify convergance
    conv_check=.False.

    !set a counter for the number of iterations we have undergone:
    i1=1

    !set add xguess to the temp_xpath list:
    temp_xpath(1,:)=xguess

    !generate the three initial points from xguess
    xinit_3=(/ xguess(1) - sqrt(3.d0)/2.d0, xguess(2) - 1.d0/2.d0 /)
    xinit_2=(/ xguess(1) + sqrt(3.d0)/2.d0, xguess(2) - 1.d0/2.d0 /)
    xinit_1=(/ xguess(1), xguess(2) + 1.d0/)

    
    !calculate the cost of the xguess and put it in temp_jpath:
    call costj(xguess,tjf)
    temp_jpath(1)=tjf


    do while(i1<itermax .and. .not. conv_check)

    
      !order the points and evaluate the cost:
      !call order_then_cost( x1,x2,x3,)--> (va,vb,vc,ja,jb,jc)
      call order_then_cost(xinit_1,xinit_2,xinit_3,va,vb,vc,ja,jb,jc)


      !initialise |ja|+|jb|+|jc|
      prev_sumj=abs(ja)+abs(jb)+abs(jc)

      
      !Calculate xm, (the midpoint of line segment bc)
      xm = (/ (vb(1) + vc(1))/2.d0, (vb(2) + vc(2))/2.d0 /)
      
      !Calculate x star
      xstar = (/ 2*xm(1) - va(1), 2*xm(2) - va(2) /)
      
      !evaluate cost at x star
      call costj(xstar, jstar)


      if(jstar<jc) then 
        !step (1) of the question, and must extend the rotated line

        xsstar = (/3*xm(1) - 2*va(1), 3*xm(2) - 2*va(2)/)
        call costj(xsstar, jsstar)
        !If j** <j* set va =x**
        if(jsstar < jstar) then
          va = xsstar
        else
          !if not then replace va with x*
          va = xstar
        end if

        ! step (2) of the question if Jc<=J*<==Ja, replace va with x*
        else if(jc<=jstar .and. jstar <= ja) then
          va=xstar
          !if jstar >ja half the length of the line segment
        else if(jstar>ja) then
          xsstar = (/ (3*xm(1) - va(1))/2.d0, (3*xm(2) - va(2))/2.d0 /)
          !evaluate x** i.e the new cost at J**
          call costj(xsstar, jsstar)

          !step 3(i) of the question: if J**<Ja set va= x**
          if(jsstar<ja) then 
            va=xsstar

            !step 3(i) of the question: rotate shrunk line segment 180 deg and evaluate J** at the new x**
          else
            xsstar = (/ (xm(1) + va(1))/2, (xm(2) + va(2))/2 /)
            call costj(xsstar, jsstar)
            !step 3(i)(I)- if J**<Ja replace va with x** 
            if (jsstar<ja) then
              va=xsstar
              !step 3(i)(II) otherwise move va and vb to the midpoints of lac and lbc respectively
            else
             va=(/ (va(1)+vc(1))/2.d0, (va(2)+vc(2))/2.d0 /)
             vb=(/ (vb(1)+vc(1))/2.d0,(vb(2)+vc(2))/2.d0 /)

          end if 

        end if
    end if

    !now we do everything else
    xinit_1=va
    xinit_2=vb
    xinit_3=vc 
    !calculate new costs:
    call costj(va,ja)
    call costj(vb,jb)
    call costj(vc,jc)

    sumj=abs(ja)+abs(jb)+abs(jc)

    !check to see if iteration criterion has been reached
    if(i1+1==itermax) then
      conv_check=.True.
    end if

    !check to see if convergance criterion has been achieved:
    call convergence_check(prev_sumj,sumj,conv_check)

    !calculate the new centroid and the cost at the new centroid:
    centroid=(va+vb+vc)/3.d0 !check thi

    !put the centroid in the list
    temp_xpath(i1+1,:)=centroid

    !calculate cost at new centroid
    call costj(centroid,tjf)
    temp_jpath(i1+1)=tjf

    !add one to the iteration number - (note this works at i1 is always 1 behind the iteration we are on in the above)
    i1=i1+1
    !reset prev_sumj
    prev_sumj=sumj
    end do 


    !create results:
    if(allocated(xpath)) then
         deallocate(xpath)
         deallocate(jpath)
    end if

    !allocate jpath and xpath
    allocate(xpath(i1,size(xguess)),jpath(i1))
    xpath=temp_xpath(1:i1,:)
    jpath=temp_jpath(1:i1)

    !set xf and jf - final output values:
    xf=xpath(i1,:)
    jf=jpath(i1)
  
  end subroutine bracket_descent









  subroutine order_then_cost(x1,x2,x3,va,vb,vc,ja,jb,jc)
    !Given an input of x1, x2, x3 order_then_cost will order the variables x1,x2,x3 into 
    !va,vb,vc where Ja>=Jb>=Jc. The function returns v_i and j_i.
    implicit none
    real(kind=8), dimension(2),intent(in) :: x1,x2,x3! input variables 
    real(kind=8), dimension(2),intent(out) :: va,vb,vc !ordered variables
    real(kind=8) :: j1,j2,j3 ! costs of unordered variables
    real(kind=8) ,intent(out) :: ja,jb,jc !costs ordered variables

    !calculate the costs of each point
    call costj(x1,j1)
    call costj(x2,j2)
    call costj(x3,j3)

    if(j1<=j2) then
      
      if(j1<=j3) then
        vc=x1
        jc=j1
        if(j2<=j3) then
          vb=x2
          jb=j2

          va=x3
          ja=j3

        else
          vb=x3
          jb=j3

          va=x2
          ja=j2

        end if 
      
      else
        vc=x3
        jc=j3

        vb=x1
        jb=j1

        va=x2
        ja=j2


      end if 
    else
      if(j3<=j2) then
        vc=x3
        jc=j3

        vb=x2
        jb=j2

        va=x1
        ja=j1
      else
        vc=x2
        jc=j2
        if(j3<=j1) then
          vb=x3
          jb=j3

          va=x1
          ja=j1

        else
          vb=x1
          jb=j1

          va=x3
          ja=j3

        end if 

      end if
    end if
    call costj(va, ja)
    call costj(vb, jb)
    call costj(vc, jc)

  end subroutine












  subroutine bd_initialize(xguess,x3,j3)
    !given xguess, generates vertices (x3) and corresponding costs (j3) for initial
    !bracket descent step
    implicit none
    real(kind=8), intent(in) :: xguess(2)
    real(kind=8), intent(out) :: j3(3),x3(3,2) !location of minimum
    integer :: i1
    real(kind=8), parameter :: l=1.d0

    x3(1,1) = xguess(1)
    x3(2,1) = xguess(1)+l*sqrt(3.d0)/2
    x3(3,1) = xguess(1)-l*sqrt(3.d0)/2
    x3(1,2) = xguess(2)+l
    x3(2,2) = xguess(2)-l/2
    x3(3,2) = xguess(2)-l/2

    do i1=1,3
      call costj(x3(i1,:),j3(i1))
    end do
  end subroutine bd_initialize





  subroutine convergence_check(j1,j2,flag_converged)
    !check if costs j1 and j2 satisfy convergence criteria
    implicit none
    real(kind=8), intent(in) :: j1,j2
    real(kind=8) :: test
    logical, intent(out) :: flag_converged

    test = abs(j1-j2)/max(abs(j1),abs(j2),1.d0)
    if (test .le. tol) then
      flag_converged = .True.
    else
      flag_converged = .False.
    end if
  end subroutine convergence_check




  subroutine newtonstep(xguess,h)
    !subroutine used for calculating the single step h of Newtons method, used in Bracket descent to compare the algorithm 
    !to Netwons method (works in exactly the same way as newton!)
    implicit none 
    real(kind=8), dimension(:), intent(in) :: xguess
    integer, dimension(size(xguess)) :: IPIV
    real(kind=8) :: temp_hess(size(xguess),size(xguess)) ! This is the temporary hessian for lapack
    real(kind=8) :: temp_grad(size(xguess))   ! THis is the temporary gradient for lapack
    real(kind=8),intent(out) ::  h(size(xguess))     ! this will be the solution h
    real(kind=8) :: cur_hess(size(xguess),size(xguess)) ! current hessian matrix
    real(kind=8) :: cur_grad(size(xguess)) ! current gradient vector
    real(kind=8) ,dimension(size(xguess)) :: tempx
    integer :: INFO
    
    tempx=xguess
    !calculate the hessian at the current point:

    call costj_hess2d(tempx,cur_hess)
    call costj_grad2d(tempx,cur_grad)
     
    !solve H(x)h=-g(x)
    !set up the hessian and gradient vectors for lapack
    temp_hess=cur_hess
    temp_grad=cur_grad

    !solve the linear system
    call dgesv(size(xguess), 1, temp_hess, size(xguess), IPIV, temp_grad, size(xguess), INFO) 
    h = -temp_grad(1:size(xguess)) !extract the solution from lapack
  end subroutine newtonstep




end module hw2


